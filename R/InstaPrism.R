
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


#' function to make refPhi from single cell
#' @param sc_Expr single cell expression with genes in rows and cells in columns
#' @param cell.type.labels a character vector indicating cell types of each cell
#' @param cell.state.labels a character vector indicating cell state of each cell
#' @param pseudo.min pseudo.min value to replace zero for normalization. Default = 1E-8
#' @export
refPrepare <-function(sc_Expr,cell.type.labels,cell.state.labels,pseudo.min=1E-8){
  if(all(is.na(rownames(sc_Expr)))){
    stop('sc_Expr must have rownames')
  }

  sc_Expr <- sc_Expr[Matrix::rowSums(sc_Expr)>0,]

  cell.state.labels = as.vector(cell.state.labels)
  cell.type.labels = as.vector(cell.type.labels)

  # collapse reference
  ref.cs=build_ref_matrix(sc_Expr,cell.state.labels)

  # normalize
  ref.cs=norm.to.one(t(ref.cs),pseudo.min) %>% t()

  map=list()
  for (ct in unique(cell.type.labels)){
    i=which(cell.type.labels==ct)
    map[[ct]]=unique(cell.state.labels[i])
  }

  return(new('refPhi_cs',
             phi.cs=ref.cs,
             map=map))

}


#' function to prepare reference matrices and map for InstaPrism
#' @param input_type one of the following input_type is allowed: 'raw', 'refPhi' and 'refPhi_cs'.
#'       With 'raw', need to specify the following input manually: sc_Expr, bulk_Expr, cell.type.labels, cell.state.labels.
#'       With 'refPhi', need to specify a refPhi object. With 'refPhi_cs', need to specify a refPhi_cs object
#' @param bulk_Expr bulk expression matrix to deconvolute, with genes in rows and samples in columns
#' @param sc_Expr single cell expression with genes in rows and cells in columns
#' @param cell.type.labels a character vector indicating cell types of each cell
#' @param cell.state.labels a character vector indicating cell state of each cell
#' @param filter a logical variable to determine whether to filter the bulk expression or not
#' @param outlier.cut & outlier.fraction: Filter genes in mixture whose expression fraction is greater than outlier.cut (Default=0.01)
#'		in more than outlier.fraction (Default=0.1) of bulk data. Removal of outlier genes will ensure that the inference will not be
#'		dominated by outliers. These two parameters denote the same thing as in new.prism() function from BayesPrism package.
#' @param pseudo.min pseudo.min value to replace zero for normalization. Default = 1E-8
#' @param refPhi a refPhi object with pre-defined scRNA-seq reference phi and map
#' @keywords internal
#' @noRd

bpPrepare<-function(input_type=c('raw','refPhi','refPhi_cs'),
                    sc_Expr = NULL,
                    cell.type.labels = NULL,cell.state.labels = NULL,
                    bulk_Expr, filter = TRUE,
                    outlier.cut=0.01,outlier.fraction=0.1,pseudo.min=1E-8,
                    refPhi = NULL, refPhi_cs = NULL){

  if(filter==TRUE){
    bulk_Expr=t(filter.bulk.outlier(t(bulk_Expr),outlier.cut=outlier.cut,outlier.fraction=outlier.fraction))
  }

  if(input_type=='raw'){
    sc_Expr <- sc_Expr[Matrix::rowSums(sc_Expr)>0,]

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

    out = new("bpPrepareObj",
              phi.cs=ref.cs,
              phi.ct=ref.ct,
              bulk_mixture=bulk_Expr[cms,,drop=F],
              map=map)

  }else if(input_type == 'refPhi'){
    cms = intersect(rownames(refPhi@phi.cs),rownames(bulk_Expr))

    ref.cs = refPhi@phi.cs[cms,]
    ref.ct = refPhi@phi.ct[cms,]

    # normalize
    ref.cs=norm.to.one(t(ref.cs),pseudo.min) %>% t()
    ref.ct=norm.to.one(t(ref.ct),pseudo.min) %>% t()

    out = new("bpPrepareObj",
              phi.cs=ref.cs,
              phi.ct=ref.ct,
              bulk_mixture=bulk_Expr[cms,,drop=F],
              map=refPhi@map)
  }else if(input_type == 'refPhi_cs'){
    cms = intersect(rownames(refPhi_cs@phi.cs),rownames(bulk_Expr))
    ref.cs = refPhi_cs@phi.cs[cms,]

    # normalize
    ref.cs=norm.to.one(t(ref.cs),pseudo.min) %>% t()

    out = new("bpPrepareObj",
              phi.cs=ref.cs,
              phi.ct=matrix(NA),
              bulk_mixture=bulk_Expr[cms,,drop=F],
              map=refPhi_cs@map)
  }
  return(out)
}


# Fast posterior estimation of cell-state fractions and expression
fastPost.ini.cs.elite.cpp<-function(bulk_Expr,ref,n.iter,n.core=1){
  cms=intersect(rownames(bulk_Expr),rownames(ref))
  ref=ref[cms,]

  sublist = function(list){
    return(list[c(1,3)]) # use only pp and pp_pre
  }

  pboptions(type = "txt", style = 3, char = "=")
  res.list=pbapply(bulk_Expr[cms,,drop=F],2,function(x) InstaPrism:::bpFixedPointCPP(bulk=matrix(x),ref=ref,n_iter = n.iter) %>% sublist(),cl = n.core)

  return(res.list)
}

#' Run InstaPrism deconvolution
#' @description A fast version of BayesPrism deconvolution, takes single cell profiles as prior information
#'      and returns cellular composition estimation for a bulk expression matrix.
#' @param input_type one of the following input_types are allowed: 'raw', 'prism' and 'refPhi_cs'.
#'       With 'raw', need to specify the following input manually: sc_Expr, bulk_Expr, cell.type.labels, cell.state.labels.
#'       With 'prism', need to provide a prism object. With 'refPhi_cs', need to provide refPhi_cs and bulk Expr
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
#' @param refPhi_cs a refPhi_cs object with single cell reference phi information (cell state only), required when input_type='refPhi_cs'
#' @param n.iter number of iterations in InstaPrism algorithm. Default = max(100, number of cell.states * 2)
#' @param verbose a logical variable determining whether to display convergence status of the model. Default = F
#' @param convergence.plot a logical variable determining whether to visualize convergence status for cell types
#' @param max_n_per_plot max number of samples to visualize in one convergence plot
#' @param n.core number of cores to use for parallel programming. Default = 1
#'
#' @return an InstaPrism object containing posterior information for cell states and cell types
#' @export
#'
InstaPrism <-function(input_type=c('raw','prism','refPhi_cs'),
                      sc_Expr=NULL,bulk_Expr=NULL,
                      cell.type.labels=NULL,cell.state.labels=NULL,
                      filter=TRUE,
                      outlier.cut=0.01,outlier.fraction=0.1,
                      pseudo.min=1E-8,
                      prismObj=NULL,
                      refPhi_cs=NULL,
                      n.iter=NULL,
                      verbose=F,
                      convergence.plot=F,max_n_per_plot=50,
                      n.core=1){
  if(input_type == 'refPhi'){
    stop('input type "refPhi" is no longer supported since version v0.1.6, please specify a refPhi_cs object instead')
  }

  if(input_type=='raw'){

    if(any(is.null(sc_Expr),is.null(cell.type.labels),is.null(cell.state.labels),is.null(bulk_Expr))){
      stop('Need to specify all raw input objects. One or more input objects from the following is missing: sc_Expr, cell.type.labels, cell.state.labels, bulk_Expr')
    }

    if(length(commonRows(sc_Expr,bulk_Expr)) < 10){
      stop('few gene overlap detected between sc_Expr and bulk_Expr, please ensure consistent gene symbol formats')
      }

    if(is.null(n.iter)){
      n.iter = max(100, 2*length(unique(cell.state.labels)))
      }

    bp=bpPrepare(input_type='raw',
                 sc_Expr,
                 cell.type.labels,cell.state.labels,
                 bulk_Expr,filter,
                 outlier.cut,outlier.fraction,
                 pseudo.min)

    cat('deconvolution with scRNA reference phi \n')
    res.list = fastPost.ini.cs.elite.cpp(bp@bulk_mixture,bp@phi.cs,n.iter,n.core)

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
      n.iter = max(100, 2*nrow(prismObj@phi_cellState@phi))
    }

    cat('deconvolution with scRNA reference phi \n')
    res.list = fastPost.ini.cs.elite.cpp(t(prismObj@mixture),t(prismObj@phi_cellState@phi),n.iter,n.core)

    map=prismObj@map
    bulk_mixture=t(prismObj@mixture)

    initial.reference = new('initial_reference',
                            phi.cs = t(prismObj@phi_cellState@phi),
                            phi.ct = t(prismObj@phi_cellType@phi))
    cms = intersect(colnames(prismObj@mixture),colnames(prismObj@phi_cellState@phi))
    cell.states = rownames(prismObj@phi_cellState@phi)

    # key=prismObj@key

  }else if(input_type=='refPhi_cs'){
    if(any(is.null(refPhi_cs),is.null(bulk_Expr))){
      stop('Need to specify refPhi_cs and bulk_Expr when input_type = "refPhi_cs"')
    }

    if(length(commonRows(refPhi_cs@phi.cs,bulk_Expr)) < 10){
      stop('few gene overlap detected between reference and bulk_Expr, please ensure consistent gene symbol formats')
    }

    if(is.null(colnames(refPhi_cs@phi.cs))|is.null(rownames(refPhi_cs@phi.cs))){
      stop('please specify dimnames for refPhi_cs@phi.cs')
    }

    if(!setequal(colnames(refPhi_cs@phi.cs), do.call(c,refPhi_cs@map))){
      stop('please ensure the same cell state information in refPhi_cs@phi.cs and refPhi_cs@map')
    }


    if(is.null(n.iter)){
      n.iter = max(100, 2* ncol(refPhi_cs@phi.cs))
      }

    bp = bpPrepare(input_type = 'refPhi_cs',
                   bulk_Expr = bulk_Expr,
                   filter = filter,
                   outlier.cut=outlier.cut,
                   outlier.fraction=outlier.fraction,
                   pseudo.min=pseudo.min,
                   refPhi_cs = refPhi_cs)

    cat('deconvolution with scRNA reference phi \n')
    res.list = fastPost.ini.cs.elite.cpp(bp@bulk_mixture,bp@phi.cs,n.iter,n.core)

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

  theta_pre=mapply(`[[`, res.list, 2)
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
      pheatmap::pheatmap(dif_ct[,bulk_schedular[[j]],drop=F],cluster_rows = F,cluster_cols = F,color = hcl.colors(50, "OrRd") %>% rev () ,
                         breaks = seq(0,0.02, length.out=50),show_colnames=T,main = 'Convergence status for cell types')

    }
  }

  # calculate scaler
  scaler_initial_cs = function(index){
    X = matrix(bulk_mixture[,index],ncol=1)
    pp = matrix(theta[,index],ncol = 1)

    phi = initial.reference@phi.cs
    norm_factor = rowSums(phi)
    norm_factor = ifelse(norm_factor==0,pseudo.min,norm_factor)
    phi_rownormalized = sweep(phi,1,norm_factor,'/')

    intermediate = phi_rownormalized %*% pp
    intermediate = ifelse(intermediate==0,pseudo.min,intermediate)

    scaler = X/intermediate
    return(scaler)
  }

  N = ncol(bulk_mixture)

  cat('scaler calculation \n')
  pboptions(type = "txt", style = 3, char = "=")

  scaler = do.call(cbind, pblapply(1:N,scaler_initial_cs,cl = n.core))
  colnames(scaler) = colnames(bulk_mixture)

  Post.ini.cs = new('theta',theta=theta)
  Post.ini.ct = new('theta', theta = theta.ct)

  return(new('InstaPrism',
             Post.ini.cs=Post.ini.cs,
             Post.ini.ct=Post.ini.ct,
             map = map,
             initial.reference = initial.reference,
             initial.scaler = scaler))

}


#' Reconstruct cell.state specific expression from initial.reference
#' @description build cell.state specific expression, which is equivalent to Z of Post.ini.ct object from InstaPrism() function
#' @param InstaPrism_obj an InstaPrism object from InstaPrism() function
#' @param cell.type.of.interest cell.type to reconstruct
#' @param n.core number of cores to use for parallel programming. Default = 1
#' @param pseudo.min pseudo.min value to replace zero for normalization. Default = 1E-8
#'
#' @return a gene * sample matrix with cell.state specific expression
#' @export
#'
reconstruct_Z_cs_initial = function(InstaPrism_obj,
                                    cell.state.of.interest,
                                    n.core = 1,
                                    pseudo.min=1E-8){
  scaler = InstaPrism_obj@initial.scaler
  phi = InstaPrism_obj@initial.reference@phi.cs
  theta = InstaPrism_obj@Post.ini.cs@theta

  norm_factor = rowSums(phi)
  norm_factor = ifelse(norm_factor==0,pseudo.min,norm_factor)

  phi_rownormalized = sweep(phi,1,norm_factor,'/') # gene * cell.states

  z_cs_from_scaler <-function(index){
    pp = theta[cell.state.of.interest,index] # a numeric value
    intermediate = phi_rownormalized[,cell.state.of.interest,drop=F] * pp # gene * 1
    z = sweep(intermediate,1,scaler[,index],'*') # gene * 1
    return(z)
  }

  N = ncol(scaler)

  pboptions(type = "txt", style = 3, char = "=")
  reconstructed_Z = do.call(cbind,pblapply(seq(1,N),z_cs_from_scaler,cl = n.core))
  colnames(reconstructed_Z) = colnames(scaler)

  return(reconstructed_Z)
}


#' Reconstruct cell.type specific expression from initial.reference
#' @description build cell.type specific expression, which is equivalent to Z of Post.ini.ct object from InstaPrism() function
#' @param InstaPrism_obj an InstaPrism object from InstaPrism() function
#' @param cell.type.of.interest cell.type to reconstruct
#' @param n.core number of cores to use for parallel programming. Default = 1
#' @param pseudo.min pseudo.min value to replace zero for normalization. Default = 1E-8
#'
#' @return a gene * sample matrix with cell.type specific expression
#' @export
#'
reconstruct_Z_ct_initial = function(InstaPrism_obj,
                                    cell.type.of.interest,
                                    n.core = 1,
                                    pseudo.min=1E-8){

  scaler = InstaPrism_obj@initial.scaler
  phi = InstaPrism_obj@initial.reference@phi.cs
  theta = InstaPrism_obj@Post.ini.cs@theta
  map = InstaPrism_obj@map

  norm_factor = rowSums(phi)
  norm_factor = ifelse(norm_factor==0,pseudo.min,norm_factor)

  phi_rownormalized = sweep(phi,1,norm_factor,'/') # gene * cell.states

  cs = map[[cell.type.of.interest]]

  z_ct_from_scaler <-function(index){
    pp = theta[cs,index,drop=F]
    intermediate = sweep(phi_rownormalized[,cs,drop=F],2,pp,'*')
    z = sweep(intermediate,1,scaler[,index],'*')
    z = rowSums(z)
    return(z)
  }

  N = ncol(scaler)

  pboptions(type = "txt", style = 3, char = "=")
  reconstructed_Z = do.call(cbind,pblapply(1:N, z_ct_from_scaler, cl = n.core))

  colnames(reconstructed_Z) = colnames(scaler)
  return(reconstructed_Z)
}


# Summarize cell.type specific expression for reference update
Z_env_normalize = function(cell.type.of.interest,InstaPrism_obj,resolution = 'ct',n.core=1,pseudo.min=1E-8,map = NULL){
  if(resolution == 'ct'){
    Z_env = reconstruct_Z_ct_initial(InstaPrism_obj,
                                     cell.type.of.interest,
                                     n.core,
                                     pseudo.min) %>% suppressMessages()
    # Z_env: gene * samples
    Z_env_sum = matrix(rowSums(Z_env),ncol = 1)
    rownames(Z_env_sum) = rownames(Z_env)
    psi = Z_env_sum/colSums(Z_env_sum) # gene * 1
    colnames(psi) = cell.type.of.interest
  }else if(resolution == 'cs'){
    # note: using updated reference on cs level is not recommended
    cell_states = map[[cell.type.of.interest]]

    Z_wrap = function(cs){
      Z_env = reconstruct_Z_cs_initial(InstaPrism_obj,
                                       cs,
                                       n.core,
                                       pseudo.min) %>% suppressMessages()
      Z_env_sum = matrix(rowSums(Z_env),ncol = 1)
      rownames(Z_env_sum) = rownames(Z_env)
      Z_env_sum = Z_env_sum/colSums(Z_env_sum) # gene * 1
      colnames(Z_env_sum) = cs
      return(Z_env_sum)
    }

    psi=do.call(cbind,lapply(cell_states,Z_wrap))
  }

  return(psi)
}


#' Run InstaPrism update on InstaPrism object
#' @description fast and memory-saving ways to update reference and get the updated deconvolution results
#'
#' @param InstaPrism_obj an InstaPrism object from InstaPrism() function
#' @param bulk_Expr bulk expression to deconvolute (without log transformation), with genes in rows and samples in columns
#' @param n.iter number of iterations in InstaPrism algorithm. Default = 100
#' @param n.core number of cores to use for parallel programming. Default = 1
#' @param cell.types.to.update non-malignant cell.types to update. Default = 'all'. Set to NULL if only need to update the 'key' cell.type.
#' @param key name of the malignant cell type. Upon setting the key parameter, the updated malignant reference will be unique for each individual.
#'		  Set to NA if there is no malignant cells in the problem, and the updated reference will be the same for all the individuals
#' @param keep.phi either 'phi.ct' or 'phi.cs', denoting whether using phi.ct or phi.cs for cell.types that are not updated.
#'     Only applicable when there are cell-types don't need to be updated
#' @param pseudo.min pseudo.min value to replace zero for normalization. Default = 1E-8
#'
#' @return an InstaPrism_update object
#' @export
#'
InstaPrism_update = function(InstaPrism_obj,
                             bulk_Expr,
                             n.iter = 100,
                             n.core = 1,
                             cell.types.to.update = 'all',
                             key = NA,
                             keep.phi = 'phi.cs',
                             pseudo.min = 1E-08){
  stopifnot(all(rownames(InstaPrism_obj@initial.reference@phi.cs) %in% rownames(bulk_Expr)))

  map = InstaPrism_obj@map
  N = ncol(bulk_Expr)
  cell.types.env = names(map)


  if(is.na(key)){

    message('the "key" parameter is set to NA, the updated reference will be the same for all the individuals')

    wrap = function(n){
      pp = InstaPrism:::bpFixedPointCPP(matrix(bulk_Expr[,n]),psi_env,n_iter = n.iter)$pp
      return(pp)
    }

    get_scaler = function(index){
      X = matrix(bulk_Expr[,index],ncol=1)
      pp = matrix(theta[,index],ncol = 1)

      stopifnot(all.equal(rownames(theta),c(colnames(psi_env))))

      phi = psi_env
      norm_factor = rowSums(phi)
      norm_factor = ifelse(norm_factor==0,pseudo.min,norm_factor)
      phi_rownormalized = sweep(phi,1,norm_factor,'/')

      intermediate = phi_rownormalized %*% pp
      intermediate = ifelse(intermediate==0,pseudo.min,intermediate)

      scaler = X/intermediate
      return(scaler)
    }

    if(is.null(cell.types.to.update)){
      stop('please specify the cell.types to update when key = NA')
    }
    else{
      if(all(length(cell.types.to.update) == 1 & cell.types.to.update == 'all')){
        cell.types.to.update = cell.types.env
        cat('update reference for all cell types \n')
      }else{

        stopifnot(all(cell.types.to.update %in% cell.types.env))

        cat('update reference for user defined cell types \n')
      }

      # update cell.types.to.update with Z_env
      Z_gt_env_normalized = do.call(cbind,lapply(cell.types.to.update, Z_env_normalize,InstaPrism_obj,'ct',n.core,pseudo.min))
      colnames(Z_gt_env_normalized) = cell.types.to.update

      cms = intersect(rownames(bulk_Expr),rownames(Z_gt_env_normalized))
      bulk_Expr = bulk_Expr[cms,,drop=F]

      if(keep.phi =='phi.ct'){

        if(nrow(InstaPrism_obj@initial.reference@phi.ct)!=nrow(InstaPrism_obj@initial.reference@phi.cs)){
          # this suggest that InstaPrism_obj@initial.reference@phi.ct is matrix(NA)
          stop('cell.type specific reference not available')
        }

        # using scRNA based Phi ct for the remaining non-malignant cells
        cell.types.env.remaining = cell.types.env[!cell.types.env %in% cell.types.to.update]
        psi_env = cbind(Z_gt_env_normalized[cms,,drop=F],
                        InstaPrism_obj@initial.reference@phi.ct[cms,cell.types.env.remaining,drop=F])

        cat('deconvolution with the updated reference \n')
        pboptions(type = "txt", style = 3, char = "=")
        theta=do.call(cbind,pblapply(sapply(1:N, list),wrap,cl = n.core))

        colnames(theta)=colnames(bulk_Expr)
        rownames(theta)=colnames(psi_env)

      }else if(keep.phi =='phi.cs'){
        # using scRNA based phi cs for the remaining cell types
        cell.types.env.remaining = cell.types.env[!cell.types.env %in% cell.types.to.update]
        map.env.remaining = map[cell.types.env.remaining]
        cell.states.env.remaining = do.call(c,map.env.remaining) %>% unname()

        psi_env = cbind(Z_gt_env_normalized[cms,,drop=F],
                        InstaPrism_obj@initial.reference@phi.cs[cms,cell.states.env.remaining,drop=F])

        cat('deconvolution with the updated reference \n')
        pboptions(type = "txt", style = 3, char = "=")
        theta=do.call(cbind,pblapply(sapply(1:N, list),wrap,cl = n.core))

        colnames(theta)=colnames(bulk_Expr)
        rownames(theta)=colnames(psi_env)
      }

      cat('scaler calculation \n')
      pboptions(type = "txt", style = 3, char = "=")

      scaler = do.call(cbind, pblapply(1:N,get_scaler,cl = n.core))
      colnames(scaler) = colnames(bulk_Expr)
    }
    return(new('InstaPrism_update',
               theta = theta,
               psi_mal = matrix(NA),
               psi_env = psi_env,
               scaler = scaler,
               map = map,
               updated.cell.types = cell.types.to.update,
               key = key,
               keep.phi = keep.phi))
  }else{

    cell.types.env = cell.types.env[cell.types.env!=key]



    wrap=function(n){
      ref = cbind(psi_mal[,n],psi_env)
      colnames(ref)[1]=key
      pp = InstaPrism:::bpFixedPointCPP(matrix(bulk_Expr[,n]),ref,n_iter = n.iter)$pp
      return(pp)
    }

    get_scaler = function(index){
      X = matrix(bulk_Expr[,index],ncol=1)
      pp = matrix(theta[,index],ncol = 1)

      stopifnot(all.equal(rownames(theta),c(key,colnames(psi_env))))

      phi = cbind(psi_mal[,index,drop=F],psi_env)
      norm_factor = rowSums(phi)
      norm_factor = ifelse(norm_factor==0,pseudo.min,norm_factor)
      phi_rownormalized = sweep(phi,1,norm_factor,'/')

      intermediate = phi_rownormalized %*% pp
      intermediate = ifelse(intermediate==0,pseudo.min,intermediate)

      scaler = X/intermediate
      return(scaler)
    }

    cat('update reference for malignant cells \n')

    # update malignant reference, which is unique for each individual
    Z_mal = reconstruct_Z_ct_initial(InstaPrism_obj,key,n.core,pseudo.min)
    Z_mal_normalized = sweep(Z_mal,2,colSums(Z_mal),'/')

    cms = intersect(rownames(bulk_Expr),rownames(Z_mal_normalized))
    bulk_Expr = bulk_Expr[cms,,drop = F]
    psi_mal = Z_mal_normalized[cms,,drop = F]

    # non-malignant reference
    if(is.null(cell.types.to.update)){

      if(keep.phi =='phi.ct'){
        # using scRNA based Phi ct for all non-malignant cells
        psi_env = InstaPrism_obj@initial.reference@phi.ct[cms,cell.types.env,drop=F]

        cat('deconvolution with the updated reference \n')
        pboptions(type = "txt", style = 3, char = "=")
        theta=do.call(cbind,pblapply(sapply(1:N, list),wrap,cl = n.core))

        colnames(theta)=colnames(bulk_Expr)
        rownames(theta)=c(key,colnames(psi_env))

      }else if(keep.phi =='phi.cs'){
        # using scRNA based phi cs for all non-malignant cells
        map.env = map[cell.types.env]
        cell.states.env = do.call(c,map.env) %>% unname()

        psi_env = InstaPrism_obj@initial.reference@phi.cs[cms,cell.states.env,drop=F]

        cat('deconvolution with the updated reference \n')
        pboptions(type = "txt", style = 3, char = "=")
        theta=do.call(cbind,pblapply(sapply(1:N, list),wrap,cl = n.core))

        colnames(theta)=colnames(bulk_Expr)
        rownames(theta)=c(key,colnames(psi_env))
      }

    }else{

      if(length(cell.types.to.update) == 1 & cell.types.to.update == 'all'){
        cell.types.to.update = cell.types.env
        cat('update reference for all cell types \n')
      }else{

        stopifnot(all(cell.types.to.update %in% cell.types.env))

        cat('update reference for-user defined cell types \n')
      }


      # update cell.types.to.update with Z_env
      Z_gt_env_normalized = do.call(cbind,lapply(cell.types.to.update,Z_env_normalize,InstaPrism_obj,'ct',n.core,pseudo.min))
      colnames(Z_gt_env_normalized) = cell.types.to.update

      if(keep.phi =='phi.ct'){
        # using scRNA based Phi ct for the remaining non-malignant cells
        cell.types.env.remaining = cell.types.env[!cell.types.env %in% cell.types.to.update]
        psi_env = cbind(Z_gt_env_normalized[cms,,drop=F],
                        InstaPrism_obj@initial.reference@phi.ct[cms,cell.types.env.remaining,drop=F])

        cat('deconvolution with the updated reference \n')
        pboptions(type = "txt", style = 3, char = "=")
        theta=do.call(cbind,pblapply(sapply(1:N, list),wrap,cl = n.core))

        colnames(theta)=colnames(bulk_Expr)
        rownames(theta)=c(key,colnames(psi_env))

      }else if(keep.phi =='phi.cs'){
        # using scRNA based phi cs for the remaining cell types
        cell.types.env.remaining = cell.types.env[!cell.types.env %in% cell.types.to.update]
        map.env.remaining = map[cell.types.env.remaining]
        cell.states.env.remaining = do.call(c,map.env.remaining) %>% unname()

        psi_env = cbind(Z_gt_env_normalized[cms,,drop=F],
                        InstaPrism_obj@initial.reference@phi.cs[cms,cell.states.env.remaining,drop=F])

        cat('deconvolution with the updated reference \n')
        pboptions(type = "txt", style = 3, char = "=")
        theta=do.call(cbind,pblapply(sapply(1:N, list),wrap,cl = n.core))

        colnames(theta)=colnames(bulk_Expr)
        rownames(theta)=c(key,colnames(psi_env))
      }
    }

    cat('scaler calculation \n')
    pboptions(type = "txt", style = 3, char = "=")

    scaler = do.call(cbind, pblapply(1:N,get_scaler,cl = n.core))
    colnames(scaler) = colnames(bulk_Expr)

    return(new('InstaPrism_update',
               theta = theta,
               psi_mal = psi_mal,
               psi_env = psi_env,
               scaler = scaler,
               map = map,
               updated.cell.types = cell.types.to.update,
               key = key,
               keep.phi = keep.phi))
  }
}


# update map
update_map = function(map, updated.cell.types){
  if(is.null(updated.cell.types)){
    updated_map = map
  }else{
    for(ct in updated.cell.types){
      map[[ct]] = ct
    }
  }
  return(map)
}

#' merge cell.state level information to cell.type level
#' @description  InstaPrism_updated_obj contains cell.state level fraction estimation when keep.phi = 'phi.cs',
#'     this function merges cell.state level theta to cell.type level.
#' @param InstaPrism_updated_obj an InstaPrism_updated_obj from InstaPrism_update() function
#'
#' @return fraction estimation at cell.type level
#' @export
#'
merge_updated_theta = function(InstaPrism_updated_obj){

  theta = InstaPrism_updated_obj@theta
  updated.cell.types = InstaPrism_updated_obj@updated.cell.types
  key = InstaPrism_updated_obj@key
  map = InstaPrism_updated_obj@map
  keep.phi = InstaPrism_updated_obj@keep.phi

  if(keep.phi == 'phi.ct'){
    stop('merge_updated_theta() is not applied when keep.phi = "phi.ct" in InstaPrism_updated_obj')
  }


  if(is.na(key)){
    map.updated = update_map(map,updated.cell.types)
  }else{
    map.updated = update_map(map, c(key,updated.cell.types))
  }
  theta.ct = do.call(rbind,lapply(map.updated,function(x)colSums(theta[rownames(theta) %in% x,,drop=F])))
  return(theta.ct)
}

#' Reconstruct cell.type specific expression using the updated reference
#'
#' @param InstaPrism_updated_obj an InstaPrism_update object from InstaPrism_update() function
#' @param cell.type.of.interest cell.type to reconstruct
#' @param n.core number of cores to use for parallel programming. Default = 1
#' @param pseudo.min pseudo.min value to replace zero for normalization. Default = 1E-8
#'
#' @return a gene * sample matrix with cell.type specific expression
#' @export
#'
reconstruct_Z_ct_updated = function(InstaPrism_updated_obj,
                                    cell.type.of.interest,
                                    n.core = 1,
                                    pseudo.min=1E-8){
  map = InstaPrism_updated_obj@map
  scaler = InstaPrism_updated_obj@scaler
  theta = InstaPrism_updated_obj@theta
  psi_mal = InstaPrism_updated_obj@psi_mal
  psi_env = InstaPrism_updated_obj@psi_env
  key = InstaPrism_updated_obj@key
  keep.phi = InstaPrism_updated_obj@keep.phi

  N = ncol(scaler)

  if(keep.phi == 'phi.cs'){

    if(is.na(key)){
      stopifnot(all.equal(matrix(NA),psi_mal))
      map.updated = update_map(map, InstaPrism_updated_obj@updated.cell.types)
      cs = map.updated[[cell.type.of.interest]]

      Z_ct_from_updated_scaler<-function(index){
        pp = theta[cs,index,drop=F]
        phi = psi_env

        norm_factor = rowSums(phi)
        norm_factor = ifelse(norm_factor==0,pseudo.min,norm_factor)

        phi_rownormalized = sweep(phi,1,norm_factor,'/')

        intermediate = sweep(phi_rownormalized[,cs,drop=F],2,pp,'*')
        z = sweep(intermediate,1,scaler[,index],'*')
        z = rowSums(z)
        return(z)
      }
      pboptions(type = "txt", style = 3, char = "=")
      reconstructed_Z = do.call(cbind,pblapply(1:N, Z_ct_from_updated_scaler, cl = n.core))
      colnames(reconstructed_Z) = colnames(scaler)
      return(reconstructed_Z)

    }else{

      map.updated = update_map(map, c(key,InstaPrism_updated_obj@updated.cell.types))
      cs = map.updated[[cell.type.of.interest]]

      Z_ct_from_updated_scaler<-function(index){
        pp = theta[cs,index,drop=F]
        phi = cbind(psi_mal[,index,drop=F],psi_env)
        colnames(phi)[1] = key

        norm_factor = rowSums(phi)
        norm_factor = ifelse(norm_factor==0,pseudo.min,norm_factor)

        phi_rownormalized = sweep(phi,1,norm_factor,'/')

        intermediate = sweep(phi_rownormalized[,cs,drop=F],2,pp,'*')
        z = sweep(intermediate,1,scaler[,index],'*')
        z = rowSums(z)
        return(z)
      }
      pboptions(type = "txt", style = 3, char = "=")
      reconstructed_Z = do.call(cbind,pblapply(1:N, Z_ct_from_updated_scaler, cl = n.core))
      colnames(reconstructed_Z) = colnames(scaler)
      return(reconstructed_Z)
    }

  }else if(keep.phi == 'phi.ct'){
    if(is.na(key)){

      norm_factor = rowSums(psi_env)
      norm_factor = ifelse(norm_factor==0,pseudo.min,norm_factor)

      phi_rownormalized = sweep(psi_env,1,norm_factor,'/')

      Z_ct_from_updated_scaler <- function(index){
        pp = theta[cell.type.of.interest,index] # a numeric value
        intermediate = phi_rownormalized[,cell.type.of.interest,drop=F] * pp # gene * 1
        z = sweep(intermediate,1,scaler[,index],'*') # gene * 1
        return(z)
      }

      reconstructed_Z = do.call(cbind,lapply(seq(1,N),Z_ct_from_updated_scaler))
      colnames(reconstructed_Z) = colnames(scaler)
      return(reconstructed_Z)

    }else{

      Z_ct_from_updated_scaler <- function(index){
        pp = theta[cell.type.of.interest,index]  # a numeric value
        phi = cbind(psi_mal[,index,drop=F],psi_env)
        colnames(phi)[1] = key

        norm_factor = rowSums(phi)
        norm_factor = ifelse(norm_factor==0,pseudo.min,norm_factor)

        phi_rownormalized = sweep(phi,1,norm_factor,'/')
        intermediate = phi_rownormalized[,cell.type.of.interest,drop=F] * pp # gene * 1
        z = sweep(intermediate,1,scaler[,index],'*') # gene * 1

        return(z)
      }

      reconstructed_Z = do.call(cbind,lapply(seq(1,N),Z_ct_from_updated_scaler))
      colnames(reconstructed_Z) = colnames(scaler)
      return(reconstructed_Z)
    }
  }
}



#' Get 3d Z tensor
#' @description obtain sampleID by gene by cell.type/cell.state 3d Z array
#'
#' @param InstaPrism_obj An InstaPrism_obj obtained from InstaPrism() function or an InstaPrism_update_obj obtained from InstaPrism_update() function
#' @param resolution a character indicating whether to get cell.type specific gene expression or cell.state specific gene expression. Available options include 'ct' or 'cs'. Default = 'ct'
#' @param n.core number of cores to use for parallel programming. Default = 1
#'
#' @return a gene by cell.type(state) 3d Z array
#' @export
#'
get_Z_array = function(InstaPrism_obj,resolution = 'ct',n.core = 1){
  if(!resolution %in% c('ct','cs')){
    stop('invalid resolution parameter provided')
  }

  if(class(InstaPrism_obj) == 'InstaPrism'){
    if(resolution == 'cs'){
      cs = rownames(InstaPrism_obj@Post.ini.cs@theta)
      Z <- array(NA, dim = c(ncol(InstaPrism_obj@Post.ini.cs@theta), nrow(InstaPrism_obj@initial.reference@phi.cs), length(cs)))

      dimnames(Z)[[1]] = colnames(InstaPrism_obj@Post.ini.cs@theta)
      dimnames(Z)[[2]] = rownames(InstaPrism_obj@initial.reference@phi.cs)
      dimnames(Z)[[3]] = cs

      for(i in 1:length(cs)){
        Z[,,i] = t(reconstruct_Z_cs_initial(InstaPrism_obj = InstaPrism_obj, cell.state.of.interest  = cs[i]))
      }

    }else if(resolution == 'ct'){
      ct = rownames(InstaPrism_obj@Post.ini.ct@theta)

      Z <- array(NA, dim = c(ncol(InstaPrism_obj@Post.ini.cs@theta), nrow(InstaPrism_obj@initial.reference@phi.cs), length(ct)))

      dimnames(Z)[[1]] = colnames(InstaPrism_obj@Post.ini.cs@theta)
      dimnames(Z)[[2]] = rownames(InstaPrism_obj@initial.reference@phi.cs)
      dimnames(Z)[[3]] = ct

      for(i in 1:length(ct)){
        Z[,,i] = t(reconstruct_Z_ct_initial(InstaPrism_obj = InstaPrism_obj,
                                            cell.type.of.interest = ct[i],n.core = n.core))
      }

    }
  }else if(class(InstaPrism_obj) == 'InstaPrism_update'){
    if(resolution == 'cs'){
      stop('resolution cs is not supported with InstaPrism_update obj')
    }else{
      ct = names(InstaPrism_obj@map)
      Z <- array(NA, dim = c(ncol(InstaPrism_obj@theta), nrow(InstaPrism_obj@psi_env), length(ct)))
      dimnames(Z)[[1]] = colnames(InstaPrism_obj@theta)
      dimnames(Z)[[2]] = rownames(InstaPrism_obj@psi_env)
      dimnames(Z)[[3]] = ct

      for(i in 1:length(ct)){
        Z[,,i] = t(reconstruct_Z_ct_updated(InstaPrism_obj,cell.type.of.interest = ct[i],n.core = n.core))
      }
    }
  }
  return(Z)
}


#' Get InstaPrism built-in reference
#' @description Load InstaPrism built-in reference and/or download the reference locally
#'
#' @param refName Name of the built-in reference. Currently supported refName includes 'BRCA', 'CRC', 'GBM', 'LUAD','OV','RCC','SKCM'
#' @param download_locally a logical variable to determine whether to download the reference to the local working directory. Default = F
#'
#' @return an refPhi_cs object that can be pass to the InstaPrism() function with `input_type = 'refPhi_cs'`
#' @export
#'
#' @examples
#' brca_ref = InstaPrism_reference('BRCA')
InstaPrism_reference = function(refName,download_locally = F){
  url = paste0('https://github.com/humengying0907/InstaPrismSourceCode/raw/main/refPhi/',refName,'_refPhi.RDS')

  if(download_locally){
    tryCatch({
      # download reference data to working directory
      download.file(url,destfile = paste0(refName,'_refPhi.RDS'),mode = 'wb',method = 'auto')
      ref <- readRDS(paste0(refName,'_refPhi.RDS'))
      return(ref)
    }, error = function(e) {
      stop("Failed to download or read the RDS file, please provide a valid refName; if the problem persists, please visit https://github.com/humengying0907/InstaPrism to download reference file manually. Error details: ", e$message)
    })
  }else{
    tryCatch({
      # download reference data to a temp file
      temp_file = tempfile()
      download.file(url,destfile = temp_file,mode = 'wb',method = 'auto')
      ref <- readRDS(temp_file)
    }, error = function(e) {
      stop("Failed to download or read the RDS file, please provide a valid refName; if the problem persists, please visit https://github.com/humengying0907/InstaPrism to download reference file manually. Error details: ", e$message)
    })
  }
}




