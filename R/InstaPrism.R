
#' Build reference matrix
#' @description Sum up expression values over cell groups for a given single expression matrix
#' @param Expr single cell expression with genes in rows and cells in columns
#' @param cell_type_labels a character vector indicating cell types/states of each cell
#' @keywords internal
#' @noRd
build_ref_matrix<-function(Expr,cell_type_labels){
  stopifnot(ncol(Expr)==length(cell_type_labels))
  group = list()
  for(i in unique(cell_type_labels)){
    group[[i]] <- which(cell_type_labels %in% i)
  }
  C = do.call(cbind, lapply(group,function(x) Matrix::rowSums(Expr[,x,drop=F])))
  C
}

#' Fast posterior estimation of cell-state fractions and expression
#' @description Fixed-point implementation of the initial Gibbs sampling step in BayesPrism.
#'      Returns a list of posterior estimation of cell-state fraction and expression.
#' @param bulk_Expr bulk expression (un-log transformed) matrix with genes in rows and samples in columns
#' @param ref single cell reference summarized at per cell-type or per cell-state level
#' @param n.iter number of iterations.
#' @param n.core number of cores to use for parallel programming. Default = 1
#' @export
fastPost.ini.cs.cpp<-function(bulk_Expr,ref,n.iter,n.core=1){
  cms=intersect(rownames(bulk_Expr),rownames(ref))
  ref=ref[cms,]
  cl<- parallel::makeCluster(n.core)
  res.list=parallel::parApply(cl,bulk_Expr[cms,],2,function(x) InstaPrism:::bpFixedPointCPP(bulk=matrix(x),ref=ref,n_iter = n.iter))
  parallel::stopCluster(cl)

  theta=mapply(`[[`, res.list, 1)
  rownames(theta)=colnames(ref)
  N <- ncol(bulk_Expr)
  G <- length(cms)
  K <- ncol(ref)
  Z=array(NA,
          dim = c(N,G,K),
          dimnames=list(colnames(bulk_Expr), cms, colnames(ref)))
  for (n in 1:N){
    Z[n,,]=res.list[[n]]$perCell
  }

  theta_pre=mapply(`[[`, res.list, 3)
  rownames(theta_pre)=colnames(ref)

  return(list(posterior=new('posterior',theta=theta,Z=Z),
              dif_cs = abs(theta-theta_pre)))
}


#' Merge posterior information over cell states
#' @describeIn Function to merge posterior information over cell states within each cell type
#' @param jointPost.obj an S4 posterior object returned by fastPost.ini.cs() function
#' @param map a list of the format list(cell.type1=c(cell.stateA, cell.stateB), ...)
#' @export
mergeK_adapted <- function(jointPost.obj,
                           map){

  bulkID <- dimnames(jointPost.obj@Z)[[1]]
  geneID <- dimnames(jointPost.obj@Z)[[2]]
  cellType <- dimnames(jointPost.obj@Z)[[3]]
  cellType.merged <- names(map)

  N <- length(bulkID)
  G <- length(geneID)
  K <- length(cellType)
  K_merged <- length(cellType.merged)

  stopifnot(length(unlist(map)) == K)

  Z <- array(NA,
             dim = c(N, G, K_merged),
             dimnames=list(bulkID, geneID, cellType.merged))

  theta.cs=jointPost.obj@theta
  theta.ct <- do.call(rbind,lapply(map,function(x)colSums(theta.cs[rownames(theta.cs) %in% x,,drop=F])))

  for (k in 1:K_merged){
    cellType.merged.k <- names(map)[k]
    cellTypes.k <- map[[k]]
    Z[,, cellType.merged.k] <- rowSums(jointPost.obj@Z[,,cellTypes.k, drop=F], dims=2)
  }
  return(new('posterior',theta=theta.ct,Z=Z))
}


#' functions to prepare reference matrices and map
#' @param sc_Expr single cell expression with genes in rows and cells in columns
#' @param bulk_Expr bulk expression matrix to deconvolute
#' @param cell.type.labels a character vector indicating cell types of each cell
#' @param cell.state.labels a character vector indicating cell state of each cell
#' @param outlier.cut & outlier.fraction: Filter genes in mixture whose expression fraction is greater than outlier.cut (Default=0.01)
#'		in more than outlier.fraction (Default=0.1) of bulk data. Removal of outlier genes will ensure that the inference will not be
#'		dominated by outliers. These two parameters denote the same thing as in new.prism() function from BayesPrism package.
#' @param pseudo.min the desired min values to replace zero after normalization
#' @keywords internal
#' @noRd
bpPrepare<-function(sc_Expr,bulk_Expr,
                    cell.type.labels,cell.state.labels,
                    outlier.cut=0.01,outlier.fraction=0.1,pseudo.min=1E-8){

  sc_Expr <- sc_Expr[rowSums(sc_Expr)>0,]
  bulk_Expr=t(filter.bulk.outlier(t(bulk_Expr),outlier.cut=outlier.cut,outlier.fraction=outlier.fraction))

  # collapse reference
  ref.cs=build_ref_matrix(sc_Expr,cell.state.labels)
  ref.ct=build_ref_matrix(sc_Expr,cell.type.labels)

  # align reference and mixture
  cms=intersect(rownames(sc_Expr),rownames(bulk_Expr))
  ref.cs=ref.cs[cms,]
  ref.ct=ref.ct[cms,]

  # normalize
  ref.cs=norm.to.one(t(ref.cs),pseudo.min) %>% t()
  ref.ct=norm.to.one(t(ref.ct),pseudo.min) %>% t()

  map=list()
  for (ct in unique(cell.type.labels)){
    i=which(cell.type.labels==ct)
    map[[ct]]=unique(cell.state.labels[i])
  }

  return(new("bpPrepare",
             phi.cs=ref.cs,
             phi.ct=ref.ct,
             bulk_mixture=bulk_Expr[cms,],
             map=map))
}

#' function that updates the reference matrix based on initial theta
#' this function is adapted from the BayesPrism for implementation of the fixed-point results
#' @param Z a Z_ngk array with cell states merged
#' @param phi cell type reference
#' @param map a list to store the correspondence between cell states and cell types
#' @param key a charater string to denote the word that corresponds to the malignant cell type, belongs to names(map)
#'		  set to NULL if there is no malignant cells in the problem.
#' @param optimizer a character string to denote which algorithm to use
#' @param opt.control a list containing the parameters to control optimization
#' @param pseudo.min the desired min values to replace zero after normalization
#' @keywords internal
#' @noRd
updateReference_adapated <- function(Z,
                                     phi,
                                     map,
                                     key,
                                     optimizer = c("MAP","MLE"),
                                     opt.control,
                                     pseudo.min=1E-8){

  cat("Update the reference matrix ... \n")

  sigma <- opt.control$sigma
  opt.control$sigma <- NULL

  optimizer <- opt.control$optimizer
  opt.control$optimizer <- NULL

  if(is.na(key)){
    #if no reference for maligant cells
    Z_gt <- colSums(Z, dims=1)

    if(optimizer == "MAP"){
      psi <- optimize.psi (phi = phi,
                           Z_gt = Z_gt,
                           prior.num = -1 / (2* sigma ^2),
                           opt.control = opt.control)$psi
    }
    if(optimizer == "MLE"){
      psi <- optimize.psi.oneGamma (phi = phi,
                                    Z_gt = Z_gt,
                                    opt.control = opt.control)$psi
    }
    return(new("refPhi", phi = psi))
  }
  else{
    # if reference for maligant cells is present
    # get MLE for psi tumor in each bulk
    Z_ng_mal <- Z[,,key]
    if(is.null(dim(Z_ng_mal)))
      Z_ng_mal <- matrix(Z_ng_mal, nrow=1, dimnames=dimnames(Z)[1:2])

    psi_mal <- get.MLE.psi_mal(Z_ng_mal = Z_ng_mal,
                               pseudo.min = pseudo.min)

    #get MAP for psi
    cellType.env <-  names(map)[names(map)!= key]
    Z_gt_env <- colSums(Z[,,cellType.env,drop=F], dims=1)
    phi_env <- phi[cellType.env,,drop=F]

    if(optimizer == "MAP"){
      psi_env <- optimize.psi (phi = phi_env,
                               Z_gt = Z_gt_env,
                               prior.num = -1 / (2* sigma ^2),
                               opt.control = opt.control)$psi
    }
    if(optimizer == "MLE"){
      psi_env <- optimize.psi.oneGamma (phi = phi_env,
                                        Z_gt = Z_gt_env,
                                        opt.control = opt.control)$psi
    }
    return(new("refTumor",
               psi_mal = psi_mal,
               psi_env = psi_env,
               key = key))}
}

#' Fast posterior estimation of cell-type fractions and expression with updated reference
#' @description Fixed-point implementation of the updated Gibbs sampling step with the updated reference.
#'      Returns a list of posterior estimation of cell-type fraction.
#' @param bulk_Expr bulk expression matrix to deconvolute
#' @param updated_phi updated reference
#' @param n.iter number of iteration
#' @param n.core number of cores to use for parallel programming. Default = 1
#' @export
fastPost.updated.ct.cpp<-function(bulk_Expr,updated_phi,n.iter, n.core=1){
  # for refTumor only
  stopifnot(nrow(bulk_Expr)==ncol(updated_phi@psi_mal))
  psi_mal=updated_phi@psi_mal
  psi_env=updated_phi@psi_env
  N = ncol(bulk_Expr)

  wrap=function(n){
    ref=t(rbind(psi_mal[n,], psi_env))
    colnames(ref)[1]=updated_phi@key
    pp = InstaPrism:::bpFixedPointCPP(matrix(bulk_Expr[,n]),ref,n_iter = n.iter)$pp
  }

  cl<- parallel::makeCluster(n.core)
  theta=do.call(cbind,parallel::parLapply(cl,sapply(1:N, list),wrap))
  parallel::stopCluster(cl)

  colnames(theta)=colnames(bulk_Expr)
  rownames(theta)=c(updated_phi@key,rownames(psi_env))


  return(new('theta',theta=theta))
}

#' Run InstaPrism deconvolution
#' @description A fast version of BayesPrism deconvolution, takes single cell profiles as prior information
#'      and returns cellular composition estimation for a bulk expression matrix.
#' @param input_type either 'raw' or 'prism'. With 'raw', need to specify the following input manually:
#'       sc_Expr, bulk_Expr, cell.type.labels, cell.state.labels
#' @param sc_Expr single cell expression matrix, with genes in rows and cells in columns
#' @param bulk_Expr bulk expression to deconvolute, with genes in rows and samples in columns
#' @param cell.type.labels a character vector indicating cell types of each cell in sc_Expr
#' @param cell.state.labels a character vector indicating cell states of each cell in sc_Expr
#' @param outlier.cut Filter genes in bulk mixture whose expression fraction is greater than outlier.cut (Default=0.01)
#'		  in more than outlier.fraction (Default=0.1) of bulk data. Removal of outlier genes will ensure that the inference will not be
#'		  dominated by outliers. This parameter denotes the same thing as in new.prism() function from BayesPrism package.
#' @param outlier.fraction Filter genes in bulk mixture whose expression fraction is greater than outlier.cut (Default=0.01)
#'		  in more than outlier.fraction (Default=0.1) of bulk data. Removal of outlier genes will ensure that the inference will not be
#'		  dominated by outliers. This parameter denotes the same thing as in new.prism() function from BayesPrism package.
#' @param pseudo.min the desired min values to replace zero after normalization. Default = 1E-8
#' @param prismObj a Prism object, required when input_type='prism'
#' @param n.iter number of iterations in InstaPrism algorithm. Default = max(20, number of cell.states * 2)
#' @param update a logical variable to denote whether return deconvolution result with the updated theta. When True,
#'      InstaPrism will implement the updateReference module from Bayesprism. Default=F.
#' @param key a charater string to denote the word that corresponds to the malignant cell type, need to be define when update = T.
#'		  Set to NULL if there is no malignant cells in the problem
#' @param optimizer a character string to denote which algorithm to use. Can be either "MAP" or "MLE", required when update=T.
#'      This parameter denotes the same thing as in the updateReference() function from BayesPrism.
#' @param opt.control a list containing the parameters to control optimization. Set opt.control = NULL for default settings as from package BayesPrism
#' @param return.Z.cs a logical variable determining whether to return cell.state specific gene expression, use default=FALSE to save memory
#' @param return.Z.ct a logical variable determining whether to return cell.type specific gene expression, use default=FALSE to save memory
#' @param verbose a logical variable determining whether to display convergence status of the model. Default = F
#' @param convergence.plot a logical variable determining whether to visualize convergence status for cell types
#' @param n.core number of cores to use for parallel programming. Default = 1
#' @param snowfall.ncore number of cores to use for reference update. Default = 1
#'
#' @return a list containing cellular fraction estimates, summarized at both cell type level and cell state level
#' @export

InstaPrism<-function(input_type=c('raw','prism'),
                 sc_Expr=NULL,bulk_Expr=NULL,
                 cell.type.labels=NULL,cell.state.labels=NULL,
                 outlier.cut=0.01,outlier.fraction=0.1,pseudo.min=1E-8,key=NA,
                 prismObj=NULL,
                 n.iter=NULL,
                 update=F,optimizer='MAP',opt.control=NULL,
                 return.Z.cs=F,
                 return.Z.ct=F,
                 verbose=F,
                 convergence.plot=F,
                 n.core=1,
                 snowfall.ncore=1){
  if(is.null(n.iter)){
    n.iter = max(20,length(unique(cell.state.labels)))
  }
  if(input_type=='raw'){
    bp=bpPrepare(sc_Expr,bulk_Expr,cell.type.labels,cell.state.labels,outlier.cut,outlier.fraction,pseudo.min=pseudo.min)

    rl = fastPost.ini.cs.cpp(bp@bulk_mixture,bp@phi.cs,n.iter,n.core)

    Post.ini.cs= rl$posterior
    dif_cs = rl$dif_cs

    map=bp@map
    bulk_mixture=bp@bulk_mixture
    phi.ct=t(bp@phi.ct)
  }else if(input_type=='prism'){

    rl = fastPost.ini.cs.cpp(t(prismObj@mixture),t(prismObj@phi_cellState@phi),n.iter,n.core)

    Post.ini.cs=rl$posterior
    dif_cs = rl$dif_cs

    map=prismObj@map
    bulk_mixture=t(prismObj@mixture)
    phi.ct=prismObj@phi_cellType@phi
    key=prismObj@key
  }

  dif_ct = do.call(rbind,lapply(map,function(x)colSums(dif_cs[rownames(dif_cs) %in% x,,drop=F])))

  if(any(dif_ct>0.01)){
    warning("fraction estimation didn't converge for some samples, enable convergence.plot for details and try larger n.iter")
  }

  if(verbose==T){
    conv_cs=t(dif_cs) %>% as.data.frame() %>% gather(key = 'cell.state',value = 'abs_diff') %>%
      group_by(cell.state) %>% summarise(min=min(abs_diff),median = median(abs_diff),max = max(abs_diff)) %>% as.data.frame()

    conv_ct=t(dif_ct) %>% as.data.frame() %>% gather(key = 'cell.type',value = 'abs_diff') %>%
      group_by(cell.type) %>% summarise(min=min(abs_diff),median = median(abs_diff),max = max(abs_diff)) %>% as.data.frame()

    cat('==================== convergence status summary ================== \n',
        'instructions: \n',
        'the absolute difference in fraction estimates between the last two iterations is utilized as an indicator of convergence, \n',
        'with smaller values indicating convergence (usually consider abs_diff < 0.01 as convergence), \n',
        'below is a summarized convergence status for cell.states/cell.types across all the samples')

    cat('=============== convergence status summary for cell states ============== \n')
    print(conv_cs)

    cat('=============== convergence status summary for cell types ============== \n')
    print(conv_ct)
  }
  if(convergence.plot==T){
    show_colnames=ifelse(ncol(dif_ct)>100,F,T)
    pheatmap::pheatmap(dif_ct,cluster_rows = F,cluster_cols = F,color = hcl.colors(50, "OrRd") %>% rev () ,
                       breaks = seq(0,0.02, length.out=50),show_colnames=show_colnames,main = 'Convergence status for cell types')
  }

  Post.ini.ct=mergeK_adapted(Post.ini.cs,map=map)

  if(return.Z.cs==T){
    Post.ini.cs = Post.ini.cs
    if (nrow(Post.ini.cs@theta)>10 & update==T){
      warning('R memory limit warning: too much data in memory, set return.Z.cs = F to save memory')
    }

  }else{
    Post.ini.cs = new('theta',theta=Post.ini.cs@theta)
  }
  gc()

  if(update==F){
    if(return.Z.ct==F){
      Post.ini.ct = new('theta',theta = Post.ini.ct@theta)
    }
    return(new('InstaPrism',Post.ini.cs=Post.ini.cs,Post.ini.ct=Post.ini.ct))

  }else if(update==T){
    if(is.null(opt.control)){
      opt.control <- valid.opt.control(list(),snowfall.ncore)
    }

    updated_phi=updateReference_adapated(Z = Post.ini.ct@Z,
                                         phi = phi.ct,
                                         map=map,
                                         key=key,
                                         optimizer=optimizer,
                                         opt.control=opt.control,
                                         pseudo.min=pseudo.min)
    gc()

    if(is(updated_phi,'refTumor')){
      message('deconvolution with the updated reference (using refTumor)')
      Post.updated.ct=fastPost.updated.ct.cpp(bulk_mixture,updated_phi,n.iter,n.core)
    }else if(is(updated_phi,'refPhi')){
      message('deconvolution with the updated reference (using refPhi)')
      Post.updated.ct=new('theta',theta=fastPost.ini.cs.cpp(bulk_mixture,t(updated_phi@phi),n.iter,n.core)$posterior@theta)
    }
    if(return.Z.ct==F){
      Post.ini.ct = new('theta',theta = Post.ini.ct@theta)
    }
    return(new('InstaPrismExtra',Post.ini.cs=Post.ini.cs,Post.ini.ct=Post.ini.ct,Post.updated.ct=Post.updated.ct))
  }
}


