
#------------------------------------------------
# custom print output
#' @method print drjacoby_output
#' @export
print.drjacoby_output <- function(x, ...) {
  
  # print summary
  summary(x)
  
  # return invisibly
  invisible(x)
}

#------------------------------------------------
# custom summary output
#' @method summary drjacoby_output
#' @export
summary.drjacoby_output <- function(x, ...) {
  
  # print summary
  message("drjacoby output:")
  message(sprintf("%s chains", length(x)))
  message(sprintf("%s rungs", length(x[[1]]$loglike_burnin)))
  message(sprintf("%s burn-in iterations", length(x[[1]]$loglike_burnin[[1]])))
  message(sprintf("%s sampling iterations", length(x[[1]]$loglike_sampling[[1]])))
  message(sprintf("%s parameters", ncol(x[[1]]$theta_sampling[[1]])))
  
  # return invisibly
  invisible(x)
}
