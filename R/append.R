append_output <- function(current, new, phase, theta_names){
  new <- as.data.frame(new)
  colnames(new) <- c("iteration", theta_names, "loglikelihood", "logprior")
  new$phase <- phase
  
  if(!is.null(current)){
    # Append iteration count
    new$iteration <- new$iteration + max(current$iteration)
    new <- rbind(current, new)
  }
  return(new)
}
