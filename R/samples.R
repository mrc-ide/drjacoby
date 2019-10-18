#' Sample N draws from all available chains
#'
#' @param x an object of class \code{drjacoby_output}
#' @param sample_n An integer number of samples
#'
#' @return A dataframe of posterior samples
#' @export
sample_chains <- function(x, sample_n) {
  
  # check inputs
  assert_custom_class(x, "drjacoby_output")
  assert_int(sample_n, "sample_n")
  assert_gr(sample_n, 0)
  assert_leq(sample_n, nrow(x[[1]]$theta_sampling$rung1) * length(x))
  
  # Join chains
  all_chains <- dplyr::bind_rows(lapply(x, function(x){
    x$theta_sampling$rung1
  }))
  
  # Sample chains
  sampled_chains <- all_chains[seq.int(1, nrow(all_chains), length.out = sample_n),]
  rownames(sampled_chains) <- 1:nrow(sampled_chains)
  
  # Ess
  ess_est_sampled <- apply(sampled_chains, 2, ess)
  message("Effective sample size of sample has range: ", min(ess_est_sampled),
          " to ", max(ess_est_sampled), ". See function ess to estimate.")
  
  return(sampled_chains)
}

