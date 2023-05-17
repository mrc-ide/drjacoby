#' Append output from the cpp mcmc to R6 object
#'
#' @param current Any current output (from previous steps)
#' @param new NEw output from cpp mcmc call
#' @param phase Phase: tune, burn or sample
#' @param theta_names Vector of parameter names
#' @param chain Chain
append_output <- function(current, new, phase, theta_names, chain){
  new <- as.data.frame(new)
  colnames(new) <- c("iteration", theta_names, "logprior", "loglikelihood")
  new$phase <- phase
  new$chain <- chain
  
  if(!is.null(current)){
    # Append iteration count
    new$iteration <- new$iteration + max(current$iteration)
    new <- rbind(current, new)
  }
  return(new)
}
