#------------------------------------------------
#' @title Flexible Markov Chain Monte Carlo via Reparameterization
#'
#' @description Flexible Markov chain monte carlo via reparameterization using
#'   the Jacobean matrix.
#'
#' @docType package
#' @name drjacoby
NULL

#------------------------------------------------
#' @useDynLib drjacoby, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom magrittr %>%
NULL

#------------------------------------------------
# unload DLL when package is unloaded
#' @noRd
.onUnload <- function(libpath) {
  library.dynam.unload("drjacoby", libpath)  # nocov
}
