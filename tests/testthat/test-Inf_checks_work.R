test_that("Checks on likelihood errors work", {
  
  # define true parameter values
  mu_true <- 3
  sigma_true <- 2
  
  # draw example data
  data_list <- list(x = 1)
  
  # define parameters dataframe
  df_params <- define_params(name = "mu", min = -10, max = 10, init = 1)
  
  
  # define log-likelihood functions
  r_loglike1 <- function(params, data, misc) {
    Inf
  }
  
  r_loglike2 <- function(params, data, misc) {
    NA
  }
  
  r_loglike3 <- function(params, data, misc) {
    NaN
  }
  
  # define log-prior function
  r_logprior <- function(params, misc) {
    0
  }
  
  expect_error(run_mcmc(data = data_list,
                   df_params = df_params,
                   loglike = r_loglike1,
                   logprior = r_logprior,
                   burnin = 1e3,
                   samples = 1e3,
                   chains = 1))
  expect_error(run_mcmc(data = data_list,
                        df_params = df_params,
                        loglike = r_loglike2,
                        logprior = r_logprior,
                        burnin = 1e3,
                        samples = 1e3,
                        chains = 1))
  expect_error(run_mcmc(data = data_list,
                        df_params = df_params,
                        loglike = r_loglike3,
                        logprior = r_logprior,
                        burnin = 1e3,
                        samples = 1e3,
                        chains = 1))
})
