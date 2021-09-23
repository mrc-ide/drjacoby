context("test-mcmc-with-cpp-likelihood-and-prior")

test_that("Cpp likelihood and prior", {
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
  # run MCMC
  cpp_mcmc_null <- run_mcmc(data = data_list,
                            df_params = df_params,
                            loglike = "loglikenull",
                            logprior = "logprior",
                            burnin = 1e3,
                            samples = 1e3,
                            chains = 1,
                            silent = TRUE)
  
  # subset output
  pe <- dplyr::filter(cpp_mcmc_null$output, phase == "sampling") %>%
    dplyr::select(mu, sigma)
  
  # check posterior estimates
  posterior_estimate <- apply(pe, 2, median)
  expect_lt(posterior_estimate["mu"] - 6, 0.1) # should approximately follow the values in logprior_strong
  expect_lt(posterior_estimate["sigma"] - 1, 0.1) # should approximately follow the values in logprior_strong
  
  # run MCMC with null prior
  cpp_mcmc_data <- run_mcmc(data = data_list,
                            df_params = df_params,
                            loglike = "loglike",
                            logprior = "logpriornull",
                            burnin = 1e3,
                            samples = 1e4,
                            silent = TRUE)
  
  # subset output
  pe <- dplyr::filter(cpp_mcmc_data$output, phase == "sampling") %>%
    dplyr::select(mu, sigma)
  
  # check posterior estimates
  posterior_estimate2 <- apply(pe, 2, median)
  expect_lt(posterior_estimate2["mu"] - mu_true, 0.1)
  expect_lt(posterior_estimate2["sigma"] - sigma_true, 0.1)
  
  ## Multiple chains
  cpp_mcmc_chains <- run_mcmc(data = data_list,
                              df_params = df_params,
                              loglike = "loglike",
                              logprior = "logpriornull",
                              chains = 2,
                              burnin = 1e3,
                              samples = 1e4,
                              silent = TRUE)
  expect_length(cpp_mcmc_chains, 4)
  
  # subset output to chain1 and check posterior estimates
  pe <- dplyr::filter(cpp_mcmc_chains$output, phase == "sampling", chain == 1) %>%
    dplyr::select(mu, sigma)
  posterior_estimate3a <- apply(pe, 2, median)
  expect_lt(posterior_estimate3a["mu"] - mu_true, 0.1)
  expect_lt(posterior_estimate3a["sigma"] - sigma_true, 0.1)
  
  # subset output to chain2 and check posterior estimates
  pe <- dplyr::filter(cpp_mcmc_chains$output, phase == "sampling", chain == 2) %>%
    dplyr::select(mu, sigma)
  posterior_estimate3b <- apply(pe, 2, median)
  expect_lt(posterior_estimate3b["mu"] - mu_true, 0.1)
  expect_lt(posterior_estimate3b["sigma"] - sigma_true, 0.1)
  
  ## Metropolis coupling
  mcmc_out_MC <- run_mcmc(data = data_list,
                          df_params = df_params,
                          loglike = "loglike",
                          logprior = "logpriornull",
                          burnin = 1e3,
                          samples = 1e4,
                          rungs = 4,
                          silent = TRUE)
  
  # subset output
  pe <- dplyr::filter(mcmc_out_MC$output, phase == "sampling") %>%
    dplyr::select(mu, sigma)
  
  # check posterior estimates
  posterior_estimate4 <- apply(pe, 2, median)
  expect_lt(posterior_estimate4["mu"] - mu_true, 0.1)
  expect_lt(posterior_estimate4["sigma"] - sigma_true, 0.1)
})
