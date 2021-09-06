test_that("Multi function cpp likelihood and prior", {
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
  Rcpp::sourceCpp("test_input_files/multi_loglike_logprior.cpp")
  
  mcmc <- run_mcmc(data = data_list,
                   df_params = df_params,
                   loglike = "loglike",
                   logprior = "logprior",
                   burnin = 1e3,
                   samples = 1e3,
                   chains = 1,
                   silent = TRUE)
  
  # subset output
  pe <- dplyr::filter(mcmc$output, phase == "sampling") %>%
    dplyr::select(mu, sigma)
  
  # check posterior estimates
  posterior_estimate <- apply(pe, 2, median)
  expect_lt(posterior_estimate["mu"] - mu_true, 0.1)
  expect_lt(posterior_estimate["sigma"] - sigma_true, 0.1)
})
