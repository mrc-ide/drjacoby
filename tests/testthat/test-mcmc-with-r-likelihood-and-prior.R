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
  r_mcmc_null <- dj$new(
    data = data_list,
    df_params = df_params,
    loglike = r_loglike_null,
    logprior = r_logprior_strong)
  r_mcmc_null$burn(iterations = 1000L, silent = TRUE)
  r_mcmc_null$sample(iterations = 1000L, silent = TRUE)
  
  # subset output
  pe <- r_mcmc_null$output()[,c("mu", "sigma")]
  
  # check posterior estimates
  posterior_estimate <- apply(pe, 2, median)
  expect_lt(posterior_estimate["mu"] - 6, 0.1)
  expect_lt(posterior_estimate["sigma"] - 1, 0.1)
  
  expect_equal(
    names(r_mcmc_null$output()),
    c("iteration", "mu", "sigma", "logprior", "loglikelihood", "phase", "chain")
  )
  
  # run MCMC with null prior
  r_mcmc_data <- dj$new(
    data = data_list,
    df_params = df_params,
    loglike = r_loglike,
    logprior = r_logprior_null)
  r_mcmc_data$burn(iterations = 1000L, silent = TRUE)
  r_mcmc_data$sample(iterations = 1000L, silent = TRUE)
  
  # subset output
  pe <- r_mcmc_data$output()[,c("mu", "sigma")]
  
  # check posterior estimates
  posterior_estimate2 <- apply(pe, 2, median)
  expect_lt(posterior_estimate2[1] - 3, 0.25)
  expect_lt(posterior_estimate2[2] - 2, 0.25)
  
  expect_equal(
    names(r_mcmc_data$output()),
    c("iteration", "mu", "sigma", "logprior", "loglikelihood", "phase", "chain")
  )
})

#------------------------------------------------
test_that("Run with multiple chains", {
  # define data
  data_list <- list(x = rnorm(10))
  
  # define parameters dataframe with multiple init values
  df_params <- define_params(name = "mu", min = -10, max = 10, init = c(-5, 5),
                             name = "sigma", min = 0, max = Inf, init = c(1, 2))
  
  # log likelihood and log prior
  r_loglike <- function(params, data, misc) {
    return(0)
  }
  r_logprior <- function(params, misc) {
    dnorm(params["mu"], 0, sd = 1, log = TRUE) +
      dgamma(params["sigma"], shape = 1, rate = 1, log = TRUE)
  }
  
  # check for error when init does not match chain number
  expect_error(dj$new(
    data = data_list,
    df_params = df_params,
    loglike = r_loglike,
    logprior = r_logprior,
    chains = 3L),
    "length(df_params$init[[i]]) == chains is not TRUE",
    fixed = TRUE
  )
  
  # run with correct number of chains
  mcmc <- dj$new(
    data = data_list,
    df_params = df_params,
    loglike = r_loglike,
    logprior = r_logprior,
    chains = 2L)
  mcmc$burn(iterations = 100L, silent = TRUE)
  mcmc$sample(iterations = 100L, silent = TRUE)
  
  # check that first values in output match initial values over all chains
  output <- mcmc$output()
  first_it <- output[output$iteration == 1, ]
  expect_equal(first_it$mu, df_params$init[[1]])
  expect_equal(first_it$sigma, df_params$init[[2]])
  
})

test_that("All parameter transformation types", {
  # define data
  data_list <- list(x = rnorm(10))
  
  # define parameters with all transformation types, including a fixed parameter
  df_params <- define_params(name = "p1", min = -Inf, max = Inf,
                             name = "p2", min = -Inf, max = 0,
                             name = "p3", min = 0, max = Inf,
                             name = "p4", min = 0, max = 1,
                             name = "p5", min = 0, max = 0)
  
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
  expect_silent({
    mcmc <- dj$new(
      data = data_list,
      df_params = df_params,
      loglike = r_loglike,
      logprior = r_logprior,
      chains = 3L)
    mcmc$burn(iterations = 100L, silent = TRUE)
    mcmc$sample(iterations = 100L, silent = TRUE)
  })
})

#------------------------------------------------
test_that("Checks on misc input", {
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
  expect_error(dj$new(
    data = data_list,
    df_params = df_params,
    loglike = r_loglike,
    logprior = r_logprior,
    misc = 5),
    "is.list(misc) is not TRUE",
    fixed = TRUE)
  
  # expect error if misc contains element "block"
  expect_error(dj$new(
    data = data_list,
    df_params = df_params,
    loglike = r_loglike,
    logprior = r_logprior,
    misc = list("block" = 5),
    "!block %in% names(misc) is not TRUE",
    fixed = TRUE
  ))
  
})
