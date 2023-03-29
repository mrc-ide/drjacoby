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

#' @title MC coupling rates
#'
#' @export
mc_acceptance <- function(mcmc_runs, chains, rungs, burnin, samples){
  mc_accept_burnin <- (dplyr::bind_rows(
    sapply(mcmc_runs, '[', 'swap_acceptance_burnin')
  ) / burnin) |>
    dplyr::rename(value = swap_acceptance_burnin) |>
    dplyr::mutate(
      phase = "burnin",
      chain = rep(1:chains, each = (rungs - 1))
    )
  mc_accept_sampling <- (dplyr::bind_rows(
    sapply(mcmc_runs, '[', 'swap_acceptance_sampling')
  ) / samples) |>
    dplyr::rename(value = swap_acceptance_sampling) |>
    dplyr::mutate(
      phase = "sampling",
      chain = rep(1:chains, each = (rungs - 1))
    )
  dplyr::bind_rows(mc_accept_burnin, mc_accept_sampling)
}

