# function to normalize expression matrix, s.t. it sum up to one for each row, with the zero entries = pseudo.min ((adapted from BayesPrism))
norm.to.one <- function(ref,
                        pseudo.min){

  G <- ncol(ref)

  phi <- ref/rowSums(ref) * (1-pseudo.min*G) + pseudo.min

  #if the minimum value is greater than zero. simply normalize by total depth
  min.value <- apply(ref,1,min)
  which.row <- min.value>0
  if(any(which.row)){
    phi[which.row,] <- ref[which.row,,drop=F]/rowSums(ref[which.row,,drop=F])
  }


  return(phi)

}

#' function to filter bulk outliers (adapted from BayesPrism)
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
  # cat("Number of outlier genes filtered from mixture =", sum(outlier.idx),"\n")

  return(mixture)
}


