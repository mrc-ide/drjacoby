context("test-mcmc-with-r-likelihood-and-prior")

test_that("R likelihood and prior", {
  set.seed(1)
  # define true parameter values
  mu_true <- 3
  sigma_true <- 2
  # draw example data
  x <- rnorm(1000, mean = mu_true, sd = sigma_true)
  # define parameters dataframe
  df_params <- data.frame(name = c("mu", "sigma"),
                          min = c(-10, 0),
                          max = c(10, Inf),
                          init = c(5, 1))
  # Null lL
  r_loglike_null <- function(params, x) {
    return(0)
  }
  # lL
  r_loglike <- function(params, x) {
    sum(dnorm(x, mean = params[1], sd = params[2], log = TRUE))
  }
  # lPrior
  r_logprior_strong <- function(params) {
    dnorm(params[1], 6, 0.1, log = TRUE) +
      dnorm(params[2], 1, 0.1, log = TRUE)
  }
  # LPrior
  r_logprior_null <- function(params) {
    return(0)
  }
  
  r_mcmc_null <- run_mcmc(data = x,
                         df_params = df_params,
                         loglike = r_loglike_null,
                         logprior = r_logprior_strong,
                         burnin = 1e3,
                         samples = 1e3,
                         silent = TRUE)
  
  posterior_estimate <- apply(r_mcmc_null$chain1$theta_sampling$rung1, 2, median)
  expect_lt(posterior_estimate[1] - 6, 0.1)
  expect_lt(posterior_estimate[2] - 1, 0.1)
  
  r_mcmc_data <- run_mcmc(data = x,
                          df_params = df_params,
                          loglike = r_loglike,
                          logprior = r_logprior_null,
                          burnin = 1e3,
                          samples = 1e3,
                          silent = TRUE)
  posterior_estimate2 <- apply(r_mcmc_data$chain1$theta_sampling$rung1, 2, median)
  expect_lt(posterior_estimate2[1] - 3, 0.25)
  expect_lt(posterior_estimate2[2] - 2, 0.25)
  
  ## Sample chains
  sampled <- sample_chains(r_mcmc_data, 100)
  expect_type(sampled, "list")
  expect_equal(nrow(sampled), 100)  
  expect_equal(ncol(sampled), 2)
  expect_error(sample_chains(r_mcmc_data, 1000000))
  expect_error(sample_chains(r_mcmc_data, -1))
  expect_error(sample_chains(r_mcmc_data, "Q"))
  expect_error(sample_chains(1, 100))
})
