#------------------------------------------------
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

#------------------------------------------------
#' @title Estimate autocorrelation
#'
#' @inheritParams plot_autocorrelation
#' @param x Single chain.
acf_data <- function(x, lag){
  stats::acf(x, plot = FALSE, lag.max = lag)$acf
}
