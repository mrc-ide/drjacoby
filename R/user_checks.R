#------------------------------------------------
#' Check Cpp loglikelihood function
#' 
#' Peforms compilation, return type and argument type checks
#'
#' @param cpp_loglike Cpp loglikelihood function
#'
#' @export
check_likelihood_compilation <- function(cpp_loglike){
  cat(sprintf("Checking compilation..."))
  temp <- RcppXPtrUtils::cppXPtr(cpp_loglike)
  cat(sprintf("\rChecking compilation...success\n"))
  cat(sprintf("Checking types..."))
  RcppXPtrUtils::checkXPtr(temp, "SEXP", c("std::vector<double>",
                                           "std::vector<double>"))
  cat(sprintf("\rChecking types...success\n"))
}

#------------------------------------------------
#' Check Cpp logprior function
#' 
#' Peforms compilation, return type and argument type checks
#'
#' @param cpp_logprior Cpp logprior function
#'
#' @export
check_prior_compilation <- function(cpp_logprior){
  cat(sprintf("Checking compilation..."))
  temp <- RcppXPtrUtils::cppXPtr(cpp_logprior)
  cat(sprintf("\rChecking compilation...success\n"))
  cat(sprintf("Checking types..."))
  RcppXPtrUtils::checkXPtr(temp, "SEXP", "std::vector<double>")
  cat(sprintf("\rChecking types...success\n"))
}
