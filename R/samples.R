#------------------------------------------------
# return 95% quantile
#' @importFrom stats quantile
#' @noRd
quantile_95 <- function(x) {
  ret <- quantile(x, probs = c(0.025, 0.5, 0.975))
  names(ret) <- c("Q2.5", "Q50", "Q97.5")
  return(ret)
}

#------------------------------------------------
#' Sample posterior draws from all available chains
#'
#' @param x an object of class \code{drjacoby_output}.
#' @param sample_n An integer number of samples.
#' @param keep_chain_index if \code{TRUE} then the column giving the chain is
#'   retained.
#'
#' @return A data.frame of posterior samples
#' @export
sample_chains <- function(x, sample_n, keep_chain_index = FALSE) {
  
  # avoid "no visible binding" note
  phase <- chain <- iteration <- logprior <- loglikelihood <- NULL
  
  # check inputs
  assert_class(x, "drjacoby_output")
  assert_pos_int(sample_n, zero_allowed = FALSE)
  assert_single_logical(keep_chain_index)
  
  # join chains
  all_chains <- x$output |>
    filter(phase == "sampling") |>
    select(-c(iteration, phase, logprior, loglikelihood))
  if (!keep_chain_index) {
    all_chains <- all_chains |>
      select(-chain)
  }
  assert_leq(sample_n, nrow(all_chains), message = sprintf("sample_n cannot exceed the total number of samples over all chains (%s)", nrow(all_chains)))
  
  # sample chains
  sampled_chains <- all_chains[seq.int(1, nrow(all_chains), length.out = sample_n),, drop = FALSE]
  sampled_chains$sample <- 1:nrow(sampled_chains)
  
  return(sampled_chains)
}

