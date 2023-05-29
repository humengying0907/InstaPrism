
# this document contains functions required for InstaPrism_legacy() function

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

  pboptions(type = "txt", style = 3, char = "=")
  res.list=pbapply(bulk_Expr[cms,,drop=F],2,function(x) InstaPrism:::bpFixedPointCPP(bulk=matrix(x),ref=ref,n_iter = n.iter),cl = n.core)

  return(res.list)
}

#' Merge posterior information over cell states
#' @description Function to merge posterior information over cell states within each cell type
#' @param Z a 3d array of cell.state specific expression
#' @param map a list with mapping information from cell.states to cell.types
#' @keywords internal
#' @noRd
merge_Z <- function(Z, map){

  bulkID <- dimnames(Z)[[1]]
  geneID <- dimnames(Z)[[2]]
  cellState <- dimnames(Z)[[3]]
  cellType.merged <- names(map)

  N <- length(bulkID)
  G <- length(geneID)
  K <- length(cellState)
  K_merged <- length(cellType.merged)

  stopifnot(length(unlist(map)) == K)

  Z.ct <- array(NA,
                dim = c(N, G, K_merged),
                dimnames=list(bulkID, geneID, cellType.merged))

  for (k in 1:K_merged){
    cellType.merged.k <- names(map)[k]
    cellTypes.k <- map[[k]]
    Z.ct[,, cellType.merged.k] <- rowSums(Z[,,cellTypes.k, drop=F], dims=2)
  }
  return(Z.ct)
}


#' Fast posterior estimation of cell-type fractions and expression with updated reference
#' @description Fixed-point implementation of the updated Gibbs sampling step with the updated reference.
#'      Returns a list of posterior estimation of cell-type fraction.
#' @param bulk_Expr bulk expression matrix to deconvolute
#' @param updated_phi updated reference
#' @param n.iter number of iteration
#' @param n.core number of cores to use for parallel programming. Default = 1
#' @keywords internal
#' @noRd
fastPost.updated.ct.cpp<-function(bulk_Expr,updated_phi,n.iter, n.core){
  # for bpRefTumor only
  stopifnot(nrow(bulk_Expr)==ncol(updated_phi@psi_mal))
  psi_mal=updated_phi@psi_mal
  psi_env=updated_phi@psi_env
  N = ncol(bulk_Expr)

  wrap=function(n){
    ref=t(rbind(psi_mal[n,], psi_env))
    colnames(ref)[1]=updated_phi@key
    pp = InstaPrism:::bpFixedPointCPP(matrix(bulk_Expr[,n]),ref,n_iter = n.iter)$pp
  }

  pboptions(type = "txt", style = 3, char = "=")
  theta=do.call(cbind,pblapply(sapply(1:N, list),wrap,cl = n.core))

  colnames(theta)=colnames(bulk_Expr)
  rownames(theta)=c(updated_phi@key,rownames(psi_env))

  return(new('theta',theta=theta))
}


#' function that updates the reference matrix based on initial theta
#' this function is adapted from the BayesPrism for implementation of the fixed-point results
#' @param Z a Z_ngk array with cell states merged
#' @param phi cell type reference
#' @param map a list to store the correspondence between cell states and cell types
#' @param key a character string to denote the malignant cell type
#' @param optimizer a character string to denote which algorithm to use
#' @param opt.control a list containing the parameters to control optimization
#' @param pseudo.min pseudo.min value to replace zero for normalization. Default = 1E-8
#' @keywords internal
#' @noRd
updateReference_adapated <- function(Z,
                                     phi,
                                     map,
                                     key,
                                     optimizer = c("MAP","MLE"),
                                     opt.control,
                                     pseudo.min=1E-8){

  cat('\n')
  cat("Update the reference matrix \n")

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
    return(new("bpRefPhi", phi = psi))
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
    return(new("bpRefTumor",
               psi_mal = psi_mal,
               psi_env = psi_env,
               key = key))}
}


#' Run InstaPrism deconvolution (legacy version)
#' @description A fast version of BayesPrism deconvolution, takes single cell profiles as prior information
#'      and returns cellular composition estimation for a bulk expression matrix.
#' @param input_type one of the following input_types are allowed: 'raw', 'prism', 'refPhi' and 'refPhi_cs'.
#'       With 'raw', need to specify the following input manually: sc_Expr, bulk_Expr, cell.type.labels, cell.state.labels.
#'       With 'prism', need to provide a prism object. With 'refPhi', need to provide refPhi and bulk_Expr.
#'       With 'refPhi_cs', need to provide refPhi_cs and bulk Expr
#' @param sc_Expr single cell expression matrix (without log transformation), with genes in rows and cells in columns, required when input_type = T
#' @param bulk_Expr bulk expression to deconvolute (without log transformation), with genes in rows and samples in columns
#' @param cell.type.labels a character vector indicating cell types of each cell in sc_Expr
#' @param cell.state.labels a character vector indicating cell states of each cell in sc_Expr
#' @param filter a logical vatiable to determine whether to filter genes in bulk_Expr or not
#' @param outlier.cut Filter genes in bulk mixture whose expression fraction is greater than outlier.cut (Default=0.01)
#'		  in more than outlier.fraction (Default=0.01) of bulk data. Removal of outlier genes will ensure that the inference will not be
#'		  dominated by outliers. This parameter denotes the same thing as in new.prism() function from BayesPrism package.
#' @param outlier.fraction Filter genes in bulk mixture whose expression fraction is greater than outlier.cut (Default=0.01)
#'		  in more than outlier.fraction (Default=0.1) of bulk data. Removal of outlier genes will ensure that the inference will not be
#'		  dominated by outliers. This parameter denotes the same thing as in new.prism() function from BayesPrism package.
#' @param pseudo.min pseudo.min value to replace zero for normalization. Default = 1E-8
#' @param prismObj a Prism object, required when input_type='prism'
#' @param refPhi a refPhi object with single cell reference phi information, required when input_type='refPhi'
#' @param refPhi_cs a refPhi_cs object with single cell reference phi information (cell state only), required when input_type='refPhi_cs'
#' @param n.iter number of iterations in InstaPrism algorithm. Default = max(20, number of cell.states * 2)
#' @param update a logical variable to denote whether return deconvolution result with the updated theta. When True,
#'      InstaPrism will implement the updateReference module from Bayesprism. Default = F
#' @param optimizer a character string to denote which algorithm to use. Can be either "MAP" or "MLE", required when update=T.
#'      This parameter denotes the same thing as in the updateReference() function from BayesPrism.
#' @param opt.control a list containing the parameters to control optimization. Set opt.control = NULL for default settings as from package BayesPrism
#' @param key name of the malignant cell type, need to be defined when update = T. The updated malignant reference will be unique for each individual.
#'		  Set to NA if there is no malignant cells in the problem, and the updated reference will be the same for all the individuals
#' @param return.Z.cs a logical variable determining whether to return cell.state specific gene expression, use default=FALSE to save memory
#' @param return.Z.ct a logical variable determining whether to return cell.type specific gene expression, use default=FALSE to save memory
#' @param verbose a logical variable determining whether to display convergence status of the model. Default = F
#' @param convergence.plot a logical variable determining whether to visualize convergence status for cell types
#' @param max_n_per_plot max number of samples to visualize in one convergence plot
#' @param n.core number of cores to use for parallel programming. Default = 1
#' @param snowfall.ncore number of cores to use for reference update. Default = 1
#'
#' @return an InstaPrism object containing posterior information for cell states and cell types,
#'      which includes theta (fraction estimates), and Z (cellular component specific expression, optional)
#' @export
#'
InstaPrism_legacy<-function(input_type=c('raw','prism','refPhi','refPhi_cs'),
                            sc_Expr=NULL,bulk_Expr=NULL,
                            cell.type.labels=NULL,cell.state.labels=NULL,
                            filter=TRUE,
                            outlier.cut=0.01,outlier.fraction=0.1,
                            pseudo.min=1E-8,
                            prismObj=NULL,
                            refPhi=NULL,
                            refPhi_cs=NULL,
                            n.iter=NULL,
                            update=F,optimizer='MAP',opt.control=NULL,key=NA,
                            return.Z.cs=F,
                            return.Z.ct=F,
                            verbose=F,
                            convergence.plot=F,max_n_per_plot=50,
                            n.core=1,
                            snowfall.ncore=1){

  if(input_type=='raw'){

    if(any(is.null(sc_Expr),is.null(cell.type.labels),is.null(cell.state.labels),is.null(bulk_Expr))){
      stop('Need to specify all raw input objects. One or more input objects from the following is missing: sc_Expr, cell.type.labels, cell.state.labels, bulk_Expr')
    }

    if(is.null(n.iter)){
      n.iter = max(20, 2*length(unique(cell.state.labels)))
    }

    bp=bpPrepare(input_type='raw',
                 sc_Expr,
                 cell.type.labels,cell.state.labels,
                 bulk_Expr,filter,
                 outlier.cut,outlier.fraction,
                 pseudo.min)

    cat('deconvolution with scRNA reference phi \n')
    res.list = fastPost.ini.cs.cpp(bp@bulk_mixture,bp@phi.cs,n.iter,n.core)

    map=bp@map
    bulk_mixture=bp@bulk_mixture
    initial.reference = new('initial_reference',
                            phi.cs = bp@phi.cs,
                            phi.ct = bp@phi.ct)
    cms = intersect(rownames(bp@bulk_mixture),rownames(bp@phi.cs))
    cell.states = colnames(bp@phi.cs)

  }else if(input_type=='prism'){

    if(is.null(prismObj)){
      stop('Need to specify a prismObj when input_type = "prism"')
    }

    if(is.null(n.iter)){
      n.iter = max(20, 2*nrow(prismObj@phi_cellState@phi))
    }

    cat('deconvolution with scRNA reference phi \n')
    res.list = fastPost.ini.cs.cpp(t(prismObj@mixture),t(prismObj@phi_cellState@phi),n.iter,n.core)

    map=prismObj@map
    bulk_mixture=t(prismObj@mixture)

    initial.reference = new('initial_reference',
                            phi.cs = t(prismObj@phi_cellState@phi),
                            phi.ct = t(prismObj@phi_cellType@phi))
    cms = intersect(colnames(prismObj@mixture),colnames(prismObj@phi_cellState@phi))
    cell.states = rownames(prismObj@phi_cellState@phi)

    # key=prismObj@key

  }else if(input_type=='refPhi'){
    if(any(is.null(refPhi),is.null(bulk_Expr))){
      stop('Need to specify refPhi and bulk_Expr when input_type = "refPhi"')
    }

    if(is.null(n.iter)){
      n.iter = max(20, 2* ncol(refPhi@phi.cs))
    }

    bp = bpPrepare(input_type = 'refPhi',
                   bulk_Expr = bulk_Expr,
                   filter = filter,
                   outlier.cut=outlier.cut,
                   outlier.fraction=outlier.fraction,
                   pseudo.min=pseudo.min,
                   refPhi = refPhi)

    cat('deconvolution with scRNA reference phi \n')
    res.list = fastPost.ini.cs.cpp(bp@bulk_mixture,bp@phi.cs,n.iter,n.core)

    map=bp@map
    bulk_mixture=bp@bulk_mixture
    initial.reference = new('initial_reference',
                            phi.cs = bp@phi.cs,
                            phi.ct = bp@phi.ct)
    cms = intersect(rownames(bp@bulk_mixture),rownames(bp@phi.cs))
    cell.states = colnames(bp@phi.cs)

  }else if(input_type=='refPhi_cs'){
    if(any(is.null(refPhi_cs),is.null(bulk_Expr))){
      stop('Need to specify refPhi_cs and bulk_Expr when input_type = "refPhi_cs"')
    }

    if(is.null(n.iter)){
      n.iter = max(20, 2* ncol(refPhi_cs@phi.cs))
    }

    bp = bpPrepare(input_type = 'refPhi_cs',
                   bulk_Expr = bulk_Expr,
                   filter = filter,
                   outlier.cut=outlier.cut,
                   outlier.fraction=outlier.fraction,
                   pseudo.min=pseudo.min,
                   refPhi_cs = refPhi_cs)

    cat('deconvolution with scRNA reference phi \n')
    res.list = fastPost.ini.cs.cpp(bp@bulk_mixture,bp@phi.cs,n.iter,n.core)

    map=bp@map
    bulk_mixture=bp@bulk_mixture

    initial.reference = new('initial_reference',
                            phi.cs = bp@phi.cs,
                            phi.ct = matrix(NA))


    phi.ct=NA
    cms = intersect(rownames(bp@bulk_mixture),rownames(bp@phi.cs))
    cell.states = colnames(bp@phi.cs)

  }

  theta=mapply(`[[`, res.list, 1)
  rownames(theta)=cell.states

  theta_pre=mapply(`[[`, res.list, 3)
  rownames(theta_pre)=cell.states

  dif_cs = abs(theta-theta_pre)
  dif_ct = do.call(rbind,lapply(map,function(x)colSums(dif_cs[rownames(dif_cs) %in% x,,drop=F])))

  # cat('merge information from cell states to cell types \n')
  theta.ct = do.call(rbind,lapply(map,function(x)colSums(theta[rownames(theta) %in% x,,drop=F])))


  if(any(dif_ct>0.01)){
    warning("fraction estimation didn't converge for some samples, enable convergence.plot for details and try larger n.iter")
  }

  if(verbose==T){
    conv_cs=t(dif_cs) %>% as.data.frame() %>% gather(key = 'cell.state',value = 'abs_diff') %>%
      group_by(cell.state) %>% summarise(min=min(abs_diff),median = median(abs_diff),max = max(abs_diff)) %>% as.data.frame()

    conv_ct=t(dif_ct) %>% as.data.frame() %>% gather(key = 'cell.type',value = 'abs_diff') %>%
      group_by(cell.type) %>% summarise(min=min(abs_diff),median = median(abs_diff),max = max(abs_diff)) %>% as.data.frame()

    cat('\n',
        '\n')
    cat(' ********************** convergence status summary ********************** \n',
        'instructions: \n',
        'the absolute difference in fraction estimates between the last two iterations is utilized as an indicator of convergence, \n',
        'with smaller values indicating convergence (usually consider abs_diff < 0.01 as convergence), \n',
        'below is a summarized convergence status for cell.states/cell.types across all the samples \n ',
        '\n')

    cat(' *************** convergence status summary for cell states *************** \n')
    print(conv_cs)

    cat('\n')
    cat(' *************** convergence status summary for cell types *************** \n')
    print(conv_ct)

    cat('\n',
        '******************* end of convergence status summary ******************** \n',
        '\n')

  }
  if(convergence.plot==T){

    if(input_type=='prism'){
      n_bulk = nrow(prismObj@mixture)
    }else{
      n_bulk = ncol(bulk_Expr)
    }

    bulk_index = seq(1,n_bulk)
    bulk_schedular = split(bulk_index, ceiling(seq_along(bulk_index) / max_n_per_plot))

    cat('\n')
    if(length(bulk_schedular)==1){
      cat(paste('display convergence quality check plot in',length(bulk_schedular), 'plot \n'))
    }else{
      cat(paste('display convergence quality check plot in',length(bulk_schedular), 'different plots \n'))
    }

    for(j in 1:length(bulk_schedular)){
      pheatmap::pheatmap(dif_ct[,bulk_schedular[[j]]],cluster_rows = F,cluster_cols = F,color = hcl.colors(50, "OrRd") %>% rev () ,
                         breaks = seq(0,0.02, length.out=50),show_colnames=T,main = 'Convergence status for cell types')

    }
  }

  if(return.Z.cs==T){
    cat('\n')
    cat('wrap up cell state specific expression \n')
    if(update == T & length(cell.states) > 20 & ncol(bulk_mixture) >= 100){
      warning('R memory limit warning: too much data in memory, set return.Z.cs = F to save memory')
    }

    N <- ncol(bulk_mixture)
    G <- length(cms)
    K <- length(cell.states)
    Z=array(NA,
            dim = c(N,G,K),
            dimnames=list(colnames(bulk_mixture), cms, cell.states))

    pb <- txtProgressBar(min = 0,
                         max = N,
                         style = 3,
                         width = 50,
                         char = "=")

    for (n in 1:N){
      Z[n,,]=res.list[[n]]$perCell
      setTxtProgressBar(pb, n)
    }
    close(pb)

    Post.ini.cs = new('posterior',theta=theta, Z=Z)
  }else{
    Post.ini.cs = new('theta',theta=theta)
  }

  if(return.Z.ct==T){
    cat('\n')
    cat('wrap up cell type specific expression \n')
    if(return.Z.cs==T){
      Z = Post.ini.cs@Z
    }else{
      N <- ncol(bulk_mixture)
      G <- length(cms)
      K <- length(cell.states)
      Z=array(NA,
              dim = c(N,G,K),
              dimnames=list(colnames(bulk_mixture), cms, cell.states))

      pb <- txtProgressBar(min = 0,
                           max = N,
                           style = 3,
                           width = 50,
                           char = "=")

      for (n in 1:N){
        Z[n,,]=res.list[[n]]$perCell
        setTxtProgressBar(pb, n)
      }
      close(pb)

    }

    Z.ct = merge_Z(Z,map)
    gc()

    Post.ini.ct =  new('posterior',theta = theta.ct, Z=Z.ct)
  }else{
    Post.ini.ct = new('theta', theta = theta.ct)
  }

  if(update==F){
    rm(res.list)
    return(new('InstaPrism_legacy',
               Post.ini.cs=Post.ini.cs,
               Post.ini.ct=Post.ini.ct,
               map = map,
               initial.reference = initial.reference))
  }else{
    if(input_type=='refPhi_cs'){
      stop('updatePhi is not available when input_type = "refPhi_cs"')
    }

    if(is.null(opt.control)){
      opt.control <- valid.opt.control(list(),snowfall.ncore)
    }

    if(return.Z.ct==T){
      rm(res.list)
      Z.ct = Post.ini.ct@Z
    }else{
      cat('\n')
      cat('wrap up cell type specific expression \n')
      if(return.Z.cs==T){
        rm(res.list)
        Z = Post.ini.cs@Z
      }else{
        N <- ncol(bulk_mixture)
        G <- length(cms)
        K <- length(cell.states)
        Z=array(NA,
                dim = c(N,G,K),
                dimnames=list(colnames(bulk_mixture), cms, cell.states))

        pb <- txtProgressBar(min = 0,
                             max = N,
                             style = 3,
                             width = 50,
                             char = "=")

        for (n in 1:N){
          Z[n,,]=res.list[[n]]$perCell
          setTxtProgressBar(pb, n)
        }
        close(pb)

      }
      Z.ct = merge_Z(Z,map)
      rm(res.list,Z)
    }

    updated_phi=updateReference_adapated(Z = Z.ct,
                                         phi = t(initial.reference@phi.ct),
                                         map=map,
                                         key=key,
                                         optimizer=optimizer,
                                         opt.control=opt.control,
                                         pseudo.min=pseudo.min)

    gc()

    if(is(updated_phi,'bpRefTumor')){
      cat('deconvolution with the updated reference \n')

      Post.updated.ct=fastPost.updated.ct.cpp(bulk_mixture,updated_phi,n.iter,n.core)

    }else if(is(updated_phi,'bpRefPhi')){
      cat('deconvolution with the updated reference \n')
      Post.updated.ct.res.list = fastPost.ini.cs.cpp(bulk_mixture,t(updated_phi@phi),n.iter,n.core)
      updated.theta.ct = mapply(`[[`, Post.updated.ct.res.list, 1)
      rownames(updated.theta.ct)=rownames(updated_phi@phi)
      Post.updated.ct = new('theta',theta = updated.theta.ct )
    }

    gc()
    return(new('InstaPrismExtra',Post.ini.cs=Post.ini.cs,
               Post.ini.ct=Post.ini.ct,
               Post.updated.ct=Post.updated.ct,
               map = map,
               initial.reference=initial.reference,
               updated.reference=updated_phi))
  }
}


