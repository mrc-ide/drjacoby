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
  assert_gr(nrow(x[[1]]$theta_sampling$rung1), sample_n / length(x))
  
  # Number of samples from each chain
  parts <- length(x)
  d <- as.integer(sample_n / parts)
  r <- as.integer(sample_n %% parts)
  sample_n_chain <- c(rep(d + 1, r), rep(d, parts - r))
  
  # Take sample_n samples from each chain and combine
  samples <- list()
  for(i in seq_along(x)){
    samp <- seq.int(1, nrow(x[[i]]$theta_sampling$rung1), length.out = sample_n_chain[i])
    samples[[i]] <- x[[i]]$theta_sampling$rung1[samp,]
  }
  samples <- dplyr::bind_rows(samples)
  return(samples)
}

