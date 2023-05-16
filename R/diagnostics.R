estimate_dic <- function(output){
  deviance <- -2 * output[output$phase == "sample", "loglikelihood"]
  dic  <- mean(deviance) + 0.5 * var(deviance)
  return(dic)
}

estimate_ess <- function(output, parameter_names){
  input <- output[output$phase == "sample", parameter_names, drop = FALSE]
  ess <- apply(input, 2, coda::effectiveSize)
  return(ess)
}

rhat <- function(output, parameter_names, n_chains, samples){
  rhat_est <- c()
  for (p in seq_along(parameter_names)) {
    rhat_est[p] <- output[output$phase == "sample", c("chain", parameter_names[p])] |>
      gelman_rubin(chains = n_chains, samples = samples)
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

#' Gelman-Rubin statistic
#' 
#' Estimate sthe Gelman-Rubin (rhat) convergence statistic for a single parameter
#' across multiple chains. Basic method, assuming all chains are of equal length
#'  
#' @references Gelman, A., and D. B. Rubin. 1992. 
#' Inference from Iterative Simulation Using Multiple Sequences. 
#' Statistical Science 7: 457–511.
#' @references \url{https://astrostatistics.psu.edu/RLectures/diagnosticsMCMC.pdf}
#'
#' @param par_matrix Matrix (interations x chains)
#' @param chains number of chains
#' @param samples number of samples
#'
#' @return Gelman-Rubin statistic
gelman_rubin <- function(par_matrix, chains, samples){
  
  # Check that >1 chains and >1 samples
  assert_gr(chains, 1)
  assert_gr(samples, 1)
  
  # Coerce to matrix
  par_matrix <- as.data.frame(par_matrix)
  
  # Mean over all samples
  all_mean <- mean(par_matrix[,2])
  
  # Mean of each chain
  chain_mean <- tapply(par_matrix[,2], par_matrix[,1], mean)
  
  # Variance of each chain
  chain_var <- tapply(par_matrix[,2], par_matrix[,1], stats::var)
  W <- (1 / chains) * sum(chain_var)
  B <- samples / (chains - 1) * sum((chain_mean - all_mean)^2)
  V <- (1 - 1 / samples) * W + (1 / samples) * B
  round(sqrt(V / W), 4)
}
