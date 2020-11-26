context("test-mcmc-with-r-likelihood-and-prior")

test_that("R likelihood and prior", {
  set.seed(1)
  
  # define true parameter values
  mu_true <- 3
  sigma_true <- 2
  
  # draw example data
  x <- rnorm(1000, mean = mu_true, sd = sigma_true)
  data_list <- list(x = x)
  
  # define parameters dataframe
  df_params <- data.frame(name = c("mu", "sigma"),
                          min = c(-10, 0),
                          max = c(10, Inf),
                          init = c(5, 1))
  
  # Null log likelihood
  r_loglike_null <- function(params, param_i, data, misc) {
    return(0)
  }
  
  # Log likelihood
  r_loglike <- function(params, param_i, data, misc) {
    sum(dnorm(data$x, mean = params["mu"], sd = params["sigma"], log = TRUE))
  }
  
  # Log prior
  r_logprior_strong <- function(params, param_i, misc) {
    dnorm(params["mu"], 6, 0.1, log = TRUE) +
      dnorm(params["sigma"], 1, 0.1, log = TRUE)
  }
  
  # Null log prior
  r_logprior_null <- function(params, param_i, misc) {
    return(0)
  }
  
  # run MCMC
  r_mcmc_null <- run_mcmc(data = data_list,
                          df_params = df_params,
                          loglike = r_loglike_null,
                          logprior = r_logprior_strong,
                          burnin = 1e3,
                          samples = 1e3,
                          silent = TRUE)
  
  # subset output
  pe <- dplyr::filter(r_mcmc_null$output, stage == "sampling", chain == "chain1") %>%
    dplyr::select(mu, sigma)
  
  # check posterior estimates
  posterior_estimate <- apply(pe, 2, median)
  expect_lt(posterior_estimate["mu"] - 6, 0.1)
  expect_lt(posterior_estimate["sigma"] - 1, 0.1)
  
  # Test parameter names are ordered correctly
  namei <- match(df_params$name, names(r_mcmc_null$output))
  namei <- 1 + namei - min(namei)
  expect_equal(namei, 1:length(namei))
  
  # run MCMC with null prior
  r_mcmc_data <- run_mcmc(data = data_list,
                          df_params = df_params,
                          loglike = r_loglike,
                          logprior = r_logprior_null,
                          burnin = 1e3,
                          samples = 1e3,
                          silent = TRUE)
  
  # subset output
  pe <- dplyr::filter(r_mcmc_data$output, stage == "sampling", chain == "chain1") %>%
    dplyr::select(mu, sigma)
  
  # check posterior estimates
  posterior_estimate2 <- apply(pe, 2, median)
  expect_lt(posterior_estimate2[1] - 3, 0.25)
  expect_lt(posterior_estimate2[2] - 2, 0.25)
  
  # Test parameter names are ordered correctly
  namei <- match(df_params$name, names(r_mcmc_data$output))
  namei <- 1 + namei - min(namei)
  expect_equal(namei, 1:length(namei))
  
  ## Sample chains
  sampled <- sample_chains(r_mcmc_data, 100)
  expect_type(sampled, "list")
  expect_equal(nrow(sampled), 100)
  expect_equal(ncol(sampled), 3)
  expect_error(sample_chains(r_mcmc_data, 1000000))
  expect_error(sample_chains(r_mcmc_data, -1))
  expect_error(sample_chains(r_mcmc_data, "Q"))
  expect_error(sample_chains(1, 100))
  
  expect_type(r_mcmc_data, "list")
  expect_type(summary(r_mcmc_data), "list")
  
  ## Check behaviour with multiple chains
  # set variable init values
  df_params$init <- list(c(-5, 5), 1)
  
  # check for error when init does not match chain number
  expect_error(run_mcmc(data = data_list,
                        df_params = df_params,
                        loglike = r_loglike,
                        logprior = r_logprior_null,
                        burnin = 1e2,
                        samples = 1e2,
                        chains = 3,
                        silent = TRUE))
  
  # run with correct number of chains
  r_mcmc_chains <- run_mcmc(data = data_list,
                            df_params = df_params,
                            loglike = r_loglike,
                            logprior = r_logprior_null,
                            burnin = 1e2,
                            samples = 1e2,
                            chains = 2,
                            silent = TRUE)
  
  # check that first values in output match initial values over all chains
  first_it <- subset(r_mcmc_chains$output, iteration == 1)
  expect_equal(first_it$mu, df_params$init[[1]])
  expect_equal(first_it$sigma, rep(df_params$init[[2]], 2))
  
})
