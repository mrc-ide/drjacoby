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
#'
#' @return Gelman-Rubin statistic
gelman_rubin <- function(par_matrix){
  chains <- ncol(par_matrix)
  n <- nrow(par_matrix)
  # Mean over all samples
  all_mean <- mean(par_matrix)
  # Mean of each chain
  chain_mean <- apply(par_matrix, 2, mean)
  # Variance of each chain
  chain_var <- apply(par_matrix, 2, stats::var)
  W <- (1 / chains) * sum(chain_var)
  B <- n / (chains - 1) * sum((chain_mean - all_mean)^2)
  V <- (1 - 1 / n) * W + (1 / n) * B
  round(sqrt(V / W), 4)
}

#' Gelman-Rubin statistic: all parameters
#' 
#' Estimates the Gelman-Rubin (rhat) convergence statistic for a all parameters
#'
#' @param output MCMC output
#' @inheritParams run_mcmc
#'
#' @return Vector of Gelman-Rubin statistics
set_rhat <- function(output, chains){
  # Extract theta matrix for each chain
  out <- c()
  for(i in 1:chains){
    out[[i]] <- output[[i]]$theta_sampling$rung1
  }
  # Number of parameters
  n_par <- ncol(out[[1]])
  # Empty vector to store rhat estimates
  rhat <- vector(mode = "numeric", length = n_par)
  # Estimate Rhat for each parameter
  for(j in 1:n_par){
    par_matrix <- sapply(out, function(x, n){x[,n]}, n = j)
    rhat[j] <- gelman_rubin(par_matrix)
  }
  # Add parameter names
  names(rhat) <- names(out[[1]])
  return(rhat)
}

