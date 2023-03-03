#' function to validate opt.control
#' @param control a named list of parameters required to control optimization
#' @keywords internal
#' @noRd
valid.opt.control <- function(control,snowfall.ncore = 1){

  ctrl <- list(maxit = 100000, maximize = FALSE, trace = 0, eps = 1e-07,
               dowarn = TRUE, tol=0, maxNA=500, n.cores=snowfall.ncore, optimizer="MAP", sigma=2)
  namc <- names(control)

  if (!all(namc %in% names(ctrl)))
    stop("unknown names in opt.control: ", namc[!(namc %in% names(ctrl))])
  ctrl[namc] <- control

  if(! ctrl$optimizer %in% c("MAP","MLE"))
    stop("unknown names of optimizer: ", ctrl$optimizer)

  stopifnot(length(ctrl$optimizer)==1)

  if(ctrl$optimizer=="MAP"){
    if(!is.numeric(ctrl$sigma))
      stop("sigma needs to be a numeric variable")
    else{
      if(ctrl$sigma<0){
        stop("sigma needs to be positive")
      }
    }
  }

  return(ctrl)
}

#' function to normalize expression matrix, s.t. it sum up to one for each row, with the zero entries = pseudo.min
#' if no zero entries, return direct normalization (no pseudo.min returned)
#' @param ref a unnormalized matrix of dimension K*G (with rownames and colnames supplied)
#' @param pseudo.min the desired min values to replace zero after normalization
#' return a normalized matrix of the same dimension
#' @keywords internal
#' @noRd
norm.to.one <- function(ref,
                        pseudo.min){

  G <- ncol(ref)

  phi <- ref/rowSums(ref) * (1-pseudo.min*G) + pseudo.min

  #if the minimum value is greater than zero. simply normalize by total depth
  min.value <- apply(ref,1,min)
  which.row <- min.value>0
  if(any(which.row)){
    #cat("One or more cell types have all genes with non-zero expression. pseudo.min is not applied to these cell types. \n")
    phi[which.row,] <- ref[which.row,,drop=F]/rowSums(ref[which.row,,drop=F])
  }


  return(phi)

}

#' function to filter bulk outliers (integrated into new.prism)
#' @param mixture the bulk RNA-seq matrix (#of samples * # of genes).
#' @param outlier.cut & outlier.fraction: Filter genes in mixture whose expression fraction is greater than outlier.cut
#'		in more than outlier.fraction of bulk data. Removal of outlier genes will ensure that the inference will not be dominated by outliers.
#' @keywords internal
#' @noRd
filter.bulk.outlier <- function(mixture,
                                outlier.cut,
                                outlier.fraction){

  mixture.norm <- mixture / rowSums(mixture)

  outlier.idx <- colSums(mixture.norm > outlier.cut) / nrow(mixture.norm) > outlier.fraction
  mixture <- mixture[, !outlier.idx, drop=F]
  cat("Number of outlier genes filtered from mixture =", sum(outlier.idx),"\n")

  return(mixture)
}


