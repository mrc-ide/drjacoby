context("test-mcmc-with-mixed-likelihood-and-prior")

test_that("Cpp likelihood and R prior, and vice versa", {
  set.seed(1)
  
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
  Rcpp::sourceCpp("test_input_files/loglike_logprior.cpp")
  
  # R likelihood and prior
  r_loglike_null <- function(params, data, misc) {
    return(0)
  }
  r_logprior_null <- function(params, misc) {
    return(0)
  }
  
  # expect no warnings when run MCMC with both combinations of prior and
  # likelihood
  expect_silent(run_mcmc(data = data_list,
                         df_params = df_params,
                         loglike = "loglikenull",
                         logprior = r_logprior_null,
                         burnin = 1e2,
                         samples = 1e2,
                         chains = 1,
                         silent = TRUE))
  
  expect_silent(run_mcmc(data = data_list,
                         df_params = df_params,
                         loglike = r_loglike_null,
                         logprior = "logpriornull",
                         burnin = 1e2,
                         samples = 1e2,
                         chains = 1,
                         silent = TRUE))
  
})
