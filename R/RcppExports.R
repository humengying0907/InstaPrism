# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' @title bpFixedPoint function using Rcpp and Armadillo
#' @param bulk Numeric matrix
#' @param ref Numeric matrix
#' @param n_iter Number of iterations (default is 20)
#' @param ppguess Optional initial guess for proportion values
#' @return List with pp, perCell and pp_pre values
#' @export
NULL

bpFixedPointCPP <- function(bulk, ref, n_iter = 20L) {
    .Call(`_InstaPrism_bpFixedPointCPP`, bulk, ref, n_iter)
}

