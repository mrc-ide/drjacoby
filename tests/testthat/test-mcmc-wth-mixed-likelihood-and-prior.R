context("test-mcmc-with-mixed-likelihood-and-prior")

test_that("Cpp likelihood and R prior, and vice versa", {
  # define true parameter values
  mu_true <- 3
  sigma_true <- 2
  
  # draw example data
  x <- rnorm(1000, mean = mu_true, sd = sigma_true)
  data_list <- list(x = x)
  
  # define parameters dataframe
  df_params <- define_params(name = "mu", min = -10, max = 10, init = 5,
                             name = "sigma", min = 0, max = 10, init = 1)
  # Source Rcpp likehood and prior functions
  fpath <- system.file("extdata", "cpp_likelihood_and_prior_functions.cpp", package = "drjacoby")
  cpp11::cpp_source(fpath)
  
  # expect no warnings when run MCMC with both combinations of prior and
  # likelihood
  mcmc1 <- dj$new(
    data = data_list,
    df_params = df_params,
    loglike = loglike_null_cpp11,
    logprior = logprior_null,
    chains = 1L
  )
  
  expect_silent(
    mcmc1$burn(iterations = 100L, silent = TRUE)
  )
  expect_silent(
    mcmc1$sample(iterations = 100L, silent = TRUE)
  )
  
  mcmc2 <- dj$new(
    data = data_list,
    df_params = df_params,
    loglike = loglike_null,
    logprior = logprior_null_cpp11,
    chains = 1L
  )
  
  expect_silent(
    mcmc2$burn(iterations = 100L, silent = TRUE)
  )
  expect_silent(
    mcmc2$sample(iterations = 100L, silent = TRUE)
  )
})
