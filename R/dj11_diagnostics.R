dic <- function(output){
  deviance <- -2 * output[output$phase == "sampling", "loglikelihood"]
  dic  <- mean(deviance) + 0.5 * var(deviance)
  return(dic)
}

ess <- function(output, parameter_names){
  input <- output[output$phase == "sampling", parameter_names]
  ess <- apply(input, 2, coda::effectiveSize)
  return(ess)
}

rhat <- function(output, parameter_names, n_chains, samples){
  rhat_est <- c()
  for (p in seq_along(parameter_names)) {
    rhat_est[p] <- output[output$phase == "sampling", c("chain", parameter_names)] |>
      drjacoby:::gelman_rubin(chains = n_chains, samples = samples)
  }
  names(rhat_est) <- parameter_names
  return(rhat_est)
}

#' @title Estimate autocorrelation
#'
#' @export
acf_data <- function(x, lag){
  stats::acf(x, plot = FALSE, lag.max = lag)$acf
}

