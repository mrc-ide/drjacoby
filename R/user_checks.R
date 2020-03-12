#------------------------------------------------
#' Check Cpp loglikelihood function
#' 
#' Peforms compilation, return type and argument type checks
#'
#' @param cpp_loglike Cpp loglikelihood function
#'
#' @export
check_likelihood_compilation <- function(cpp_loglike){
  message(sprintf("Checking compilation..."))
  temp <- RcppXPtrUtils::cppXPtr(cpp_loglike)
  message(sprintf("\rChecking compilation...success\n"))
  message(sprintf("Checking types..."))
  RcppXPtrUtils::checkXPtr(temp, "SEXP", c("std::vector<double>",
                                           "std::vector<double>"))
  message(sprintf("\rChecking types...success\n"))
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
  message(sprintf("Checking compilation..."))
  temp <- RcppXPtrUtils::cppXPtr(cpp_logprior)
  message(sprintf("\rChecking compilation...success\n"))
  message(sprintf("Checking types..."))
  RcppXPtrUtils::checkXPtr(temp, "SEXP", "std::vector<double>")
  message(sprintf("\rChecking types...success\n"))
}
