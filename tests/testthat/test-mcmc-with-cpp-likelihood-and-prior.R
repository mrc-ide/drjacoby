test_that("Cpp likelihood and prior", {
  # define true parameter values
  mu_true <- 3
  sigma_true <- 2
  
  # draw example data
  x <- rnorm(1000, mean = mu_true, sd = sigma_true)
  data_list <- list(x = x)
  
  # define parameters dataframe
  df_params <- define_params(name = "mu", min = -10, max = 10, init = 5,
                             name = "sigma", min = 0, max = 10, init = 1)
  
  # Source Rcpp likehood
  fpath <- system.file("extdata", "cpp_likelihood_and_prior_functions.cpp", package = "drjacoby")
  cpp11::cpp_source(fpath)
  # Source Rcpp likehood and prior functions
  cpp11::cpp_source("test_input_files/loglike_logprior.cpp")
  
  # run MCMC
  cpp_mcmc_null <- dj$new(
    data = data_list,
    df_params = df_params,
    loglike = loglike_null_cpp11,
    logprior = logprior_normal_cpp11,
    chains = 1L,
    seed = 1)
  cpp_mcmc_null$burn(iterations = 100L, silent = TRUE)
  cpp_mcmc_null$sample(iterations = 1000L, silent = TRUE)
  
  # subset output
  pe <- cpp_mcmc_null$output()[,c("mu", "sigma")]
  
  # check posterior estimates
  posterior_estimate <- apply(pe, 2, median)
  expect_lt(posterior_estimate["mu"] - 6, 0.1)
  expect_lt(posterior_estimate["sigma"] - 1, 0.1)
  
  # run MCMC with null prior
  cpp_mcmc_data <- dj$new(
    data = data_list,
    df_params = df_params,
    loglike = loglike_normal_cpp11,
    logprior = logprior_null_cpp11,
    seed = 1)
  cpp_mcmc_data$burn(iterations = 100L, silent = TRUE)
  cpp_mcmc_data$sample(iterations = 1000L, silent = TRUE)
  
  # subset output
  pe <- cpp_mcmc_data$output()[,c("mu", "sigma")]
  
  # check posterior estimates
  posterior_estimate2 <- apply(pe, 2, median)
  expect_lt(posterior_estimate2["mu"] - mu_true, 0.1)
  expect_lt(posterior_estimate2["sigma"] - sigma_true, 0.1)
  
  ## Multiple chains
  cpp_mcmc_chains <- dj$new(
    data = data_list,
    df_params = df_params,
    loglike = loglike_normal_cpp11,
    logprior = logprior_null_cpp11,
    chains = 2L,
    seed = 1)
  cpp_mcmc_chains$burn(iterations = 100L, silent = TRUE)
  cpp_mcmc_chains$sample(iterations = 1000L, silent = TRUE)
  
  # subset output to chain1 and check posterior estimates
  pe <- cpp_mcmc_chains$output(chain = 1)[,c("mu", "sigma")]
  posterior_estimate3a <- apply(pe, 2, median)
  expect_lt(posterior_estimate3a["mu"] - mu_true, 0.1)
  expect_lt(posterior_estimate3a["sigma"] - sigma_true, 0.1)
  
  # subset output to chain2 and check posterior estimates
  pe <- cpp_mcmc_chains$output(chain = 2)[,c("mu", "sigma")]
  posterior_estimate3b <- apply(pe, 2, median)
  expect_lt(posterior_estimate3b["mu"] - mu_true, 0.1)
  expect_lt(posterior_estimate3b["sigma"] - sigma_true, 0.1)
  
  ## Metropolis coupling
  mcmc_out_MC <- dj$new(
    data = data_list,
    df_params = df_params,
    loglike = loglike_normal_cpp11,
    logprior = logprior_null_cpp11,
    seed = 1)
  mcmc_out_MC$tune(beta = seq(1, 0, -0.1), silent = TRUE)
  mcmc_out_MC$burn(iterations = 100L, silent = TRUE)
  mcmc_out_MC$sample(iterations = 1000L, silent = TRUE)
  
  # subset output
  pe <- mcmc_out_MC$output()[,c("mu", "sigma")]
  
  # check posterior estimates
  posterior_estimate4 <- apply(pe, 2, median)
  expect_lt(posterior_estimate4["mu"] - mu_true, 0.1)
  expect_lt(posterior_estimate4["sigma"] - sigma_true, 0.1)
})
