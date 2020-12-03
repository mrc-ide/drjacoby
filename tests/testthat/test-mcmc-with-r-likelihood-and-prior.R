context("test-mcmc-with-r-likelihood-and-prior")

#------------------------------------------------
test_that("define_params() working as expected", {
  
  # expect error if bad name
  expect_error(define_params(name = "mu", min = -10, max = 10, foo = 5))
  
  # expect error if duplicated names
  expect_error(define_params(name = "mu", min = -10, max = 10,
                             name = "mu", min = -10, max = 10))
  
  # expect error if missing some args
  expect_error(define_params(name = "mu1", min = -10, max = 10, init = 5,
                             name = "mu2", min = -10, max = 10))
  
  # expect blocks in correct format
  df_params <- define_params(name = "mu", min = -10, max = 10, block = 1,
                             name = "sigma", min = 0, max = Inf, block = 1)
  expect_equal(df_params$block, list(block = 1, block = 1))
  
  # correctly defined dataframe
  df_params <- define_params(name = "mu", min = -10, max = 10, init = 5,
                             name = "sigma", min = 0, max = Inf, init = 1)
  
  # same parameters dataframe using base R method
  df_params_base <- data.frame(name = c("mu", "sigma"),
                               min = c(-10, 0),
                               max = c(10, Inf))
  df_params_base$init <- list(init = 5, init = 1)
  
  # check identical
  expect_identical(df_params, df_params_base)
  
})

#------------------------------------------------
test_that("R likelihood and prior", {
  set.seed(1)
  
  # define true parameter values
  mu_true <- 3
  sigma_true <- 2
  
  # draw example data
  data_list <- list(x = rnorm(1000, mean = mu_true, sd = sigma_true))
  
  # define parameters dataframe
  df_params <- define_params(name = "mu", min = -10, max = 10, init = 5,
                             name = "sigma", min = 0, max = Inf, init = 1)
  
  # Log likelihood and log prior
  r_loglike <- function(params, data, misc) {
    sum(dnorm(data$x, mean = params["mu"], sd = params["sigma"], log = TRUE))
  }
  r_logprior_strong <- function(params, misc) {
    dnorm(params["mu"], 6, 0.1, log = TRUE) +
      dnorm(params["sigma"], 1, 0.1, log = TRUE)
  }
  
  # Null log likelihood and log prior
  r_loglike_null <- function(params, data, misc) {
    return(0)
  }
  r_logprior_null <- function(params, misc) {
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
  pe <- dplyr::filter(r_mcmc_null$output, phase == "sampling") %>%
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
  pe <- dplyr::filter(r_mcmc_data$output, phase == "sampling") %>%
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
  
})

#------------------------------------------------
test_that("Run with multiple chains", {
  set.seed(1)
  
  # define data
  data_list <- list(x = rnorm(10))
  
  # define parameters dataframe with multiple init values
  df_params <- define_params(name = "mu", min = -10, max = 10, init = c(-5, 5),
                             name = "sigma", min = 0, max = Inf, init = 1)
  
  # log likelihood and log prior
  r_loglike <- function(params, data, misc) {
    return(0)
  }
  r_logprior <- function(params, misc) {
    dnorm(params["mu"], 0, sd = 1, log = TRUE) +
      dgamma(params["sigma"], shape = 1, rate = 1, log = TRUE)
  }
  
  # check for error when init does not match chain number
  expect_error(run_mcmc(data = data_list,
                        df_params = df_params,
                        loglike = r_loglike,
                        logprior = r_logprior,
                        burnin = 1e2,
                        samples = 1e2,
                        chains = 3,
                        silent = TRUE))
  
  # run with correct number of chains
  mcmc <- run_mcmc(data = data_list,
                   df_params = df_params,
                   loglike = r_loglike,
                   logprior = r_logprior,
                   burnin = 1e2,
                   samples = 1e2,
                   chains = 2,
                   silent = TRUE)
  
  # check that first values in output match initial values over all chains
  first_it <- subset(mcmc$output, iteration == 1)
  expect_equal(first_it$mu, df_params$init[[1]])
  expect_equal(first_it$sigma, rep(df_params$init[[2]], 2))
  
})

#------------------------------------------------
test_that("All parameter transformation types", {
  set.seed(1)
  
  # define data
  data_list <- list(x = rnorm(10))
  
  # define parameters with all transformation types
  df_params <- define_params(name = "p1", min = -Inf, max = Inf,
                             name = "p2", min = -Inf, max = 0,
                             name = "p3", min = 0, max = Inf,
                             name = "p4", min = 0, max = 1)
  
  # log likelihood and log prior
  r_loglike <- function(params, data, misc) {
    return(0)
  }
  r_logprior <- function(params, misc) {
    dnorm(params["p1"], 0, sd = 1, log = TRUE) +
      dgamma(-params["p2"], shape = 1, rate = 1, log = TRUE) +
      dgamma(params["p3"], shape = 1, rate = 1, log = TRUE) +
      dbeta(params["p4"], shape1 = 1, shape2 = 1, log = TRUE)
  }
  
  # should run without error
  mcmc <- run_mcmc(data = data_list,
                   df_params = df_params,
                   loglike = r_loglike,
                   logprior = r_logprior,
                   burnin = 1e2,
                   samples = 1e2,
                   chains = 3,
                   silent = TRUE)
  
  # should run without error
  expect_silent(run_mcmc(data = data_list,
                         df_params = df_params,
                         loglike = r_loglike,
                         logprior = r_logprior,
                         burnin = 1e2,
                         samples = 1e2,
                         chains = 3,
                         silent = TRUE))
  
})

#------------------------------------------------
test_that("Checks on misc input", {
  set.seed(1)
  
  # define data
  data_list <- list(x = rnorm(10))
  
  # define parameters
  df_params <- define_params(name = "mu", min = -Inf, max = Inf)
  
  # log likelihood and log prior
  r_loglike <- function(params, data, misc) {
    return(0)
  }
  r_logprior <- function(params, misc) {
    dnorm(params["mu"], 0, sd = 1, log = TRUE)
  }
  
  # expect error if misc not a list
  expect_error(run_mcmc(data = data_list,
                        df_params = df_params,
                        loglike = r_loglike,
                        logprior = r_logprior,
                        misc = 5,
                        burnin = 1e2,
                        samples = 1e2))
  
  # expect error if misc contains element "block"
  expect_error(run_mcmc(data = data_list,
                        df_params = df_params,
                        loglike = r_loglike,
                        logprior = r_logprior,
                        misc = list("block" = 5),
                        burnin = 1e2,
                        samples = 1e2))
  
})
