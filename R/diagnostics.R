estimate_dic <- function(output){
  deviance <- -2 * output[output$phase == "sample", "loglikelihood"]
  dic  <- mean(deviance) + 0.5 * stats::var(deviance)
  return(dic)
}

estimate_ess <- function(output, parameter_names){
  input <- output[output$phase == "sample", parameter_names, drop = FALSE]
  ess <- apply(input, 2, coda::effectiveSize)
  return(ess)
}

estimate_rhat <- function(output, parameter_names, n_chains, samples){
  rhat_est <- c()
  for (p in seq_along(parameter_names)) {
    rhat_est[p] <- output[output$phase == "sample", c("chain", parameter_names[p])] |>
      gelman_rubin(chains = n_chains, samples = samples[3,1])
  }
  names(rhat_est) <- parameter_names
  return(rhat_est)
}

estimate_acceptance_rate <- function(acceptance_counter, iteration_counter, chains, theta_names){
  # Extract only cold chain rates
  for(i in 1:chains){
    for(j in 1:3){
      acceptance_counter[[i]][[j]] <- acceptance_counter[[i]][[j]][1,]
    }
  }
  ar <- lapply(acceptance_counter, list_r_bind)
  for(i in 1:chains){
    its <- iteration_counter[,i]
    ar[[i]] <- apply(ar[[i]], 2, function(x, y){
      x / y
    }, y = its)
    colnames(ar[[i]]) <- theta_names
  }
  names(ar) <- paste0("Chain_", 1:chains)
  return(ar)
}

estimate_timing  <- function(seconds, iterations){
  seconds <- rbind(seconds, Total = colSums(seconds))
  iterations <- rbind(iterations, All = colSums(iterations))
  iterations_per_second <- round(iterations / seconds)
  return(
    list(
      seconds = seconds,
      iterations_per_second = iterations_per_second)
  )
}

estimate_mc_acceptance_rate <- function(swap_acceptance_counter, iteration_counter){
  swap_acceptance_rate <- apply(swap_acceptance_counter, 2, function(x){
    x / iteration_counter
  })
  rownames(swap_acceptance_rate) <- rownames(swap_acceptance_counter)
  swap_acceptance_rate  
}


#' @title Estimate autocorrelation
#'
#' @param x samples
#' @param lag lag
#' @export
acf_data <- function(x, lag){
  stats::acf(x, plot = FALSE, lag.max = lag)$acf
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
  stopifnot(chains > 1)
  stopifnot(samples > 1)
  
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
