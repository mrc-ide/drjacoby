# return 95% quantile
#' @importFrom stats quantile
#' @noRd
quantile_95 <- function(x) {
  ret <- quantile(x, probs = c(0.025, 0.5, 0.975))
  names(ret) <- c("Q2.5", "Q50", "Q97.5")
  return(ret)
}

#' Sample N draws from all available chains
#'
#' @param x an object of class \code{drjacoby_output}
#' @param sample_n An integer number of samples
#'
#' @return A data.frame of posterior samples
#' @export
sample_chains <- function(x, sample_n) {
  
  # check inputs
  assert_class(x, "drjacoby_output")
  assert_int(sample_n, "sample_n")
  assert_gr(sample_n, 0)
  
  # Join chains
  all_chains <- dplyr::filter(x$output, .data$phase == "sampling") %>%
    dplyr::select(-.data$chain, -.data$iteration, -.data$phase, -.data$logprior, -.data$loglikelihood)
  assert_leq(sample_n, nrow(all_chains))
  
  # Sample chains
  sampled_chains <- all_chains[seq.int(1, nrow(all_chains), length.out = sample_n),, drop = FALSE]
  sampled_chains$sample <- 1:nrow(sampled_chains)
  
  # ESS
  ess_est_sampled <- round(apply(sampled_chains[,1:(ncol(sampled_chains) - 1), drop = FALSE], 2, coda::effectiveSize))
  message("Effective sample size of sample has range: ", min(ess_est_sampled),
          " to ", max(ess_est_sampled), ". See function ess to estimate.")
  
  return(sampled_chains)
}

