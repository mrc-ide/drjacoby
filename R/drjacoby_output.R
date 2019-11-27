
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
summary.drjacoby_output <- function(object, ...) {
  
  # print summary
  message("drjacoby output:")
  message(sprintf("%s chains", object$parameters$chains))
  message(sprintf("%s rungs", object$parameters$rungs))
  message(sprintf("%s burn-in iterations", object$parameters$burnin))
  message(sprintf("%s sampling iterations", object$parameters$samples))
  message(sprintf("%s parameters", length(object$parameters$df_params$name)))
  
  # return invisibly
  invisible(object)
}
