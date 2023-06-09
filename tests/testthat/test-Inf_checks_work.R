test_that("Checks on likelihood and prior errors work", {
  
  # define true parameter values
  mu_true <- 3
  sigma_true <- 2
  
  # draw example data
  data_list <- list(x = 1)
  
  # define parameters dataframe
  df_params <- data.frame() |>
    add_parameter(name = "mu", min = -10, max = 10, init = 1)
  
  
  # define failing log-likelihood functions
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
  
  mcmc <- dj$new(
    data = data_list,
    df_params = df_params,
    loglike = r_loglike1,
    logprior = r_logprior
  )
  
  expect_error(
    mcmc$burn(iterations = 10L, silent = TRUE),
    "Error in mcmc, check $error_debug",
    fixed = TRUE
  )
  
  mcmc <- dj$new(
    data = data_list,
    df_params = df_params,
    loglike = r_loglike2,
    logprior = r_logprior
  )
  
  expect_error(
    mcmc$burn(iterations = 10L, silent = TRUE),
    "Error in mcmc, check $error_debug",
    fixed = TRUE
  )
  mcmc <- dj$new(
    data = data_list,
    df_params = df_params,
    loglike = r_loglike3,
    logprior = r_logprior
  )
  
  expect_error(
    mcmc$burn(iterations = 10L, silent = TRUE),
    "Error in mcmc, check $error_debug",
    fixed = TRUE
  )
  
  
  # define failing log-prior functions
  r_loglike <- function(params, data, misc) {
    0
  }
  
  r_logprior1 <- function(params, data, misc) {
    Inf
  }
  
  r_logprior2 <- function(params, data, misc) {
    NA
  }
  
  r_logprior3 <- function(params, data, misc) {
    NaN
  }
  
  mcmc <- dj$new(
    data = data_list,
    df_params = df_params,
    loglike = r_loglike,
    logprior = r_logprior1
  )
  
  expect_error(
    mcmc$burn(iterations = 10L, silent = TRUE),
    "Error in mcmc, check $error_debug",
    fixed = TRUE
  )
  
  mcmc <- dj$new(
    data = data_list,
    df_params = df_params,
    loglike = r_loglike,
    logprior = r_logprior2
  )
  
  expect_error(
    mcmc$burn(iterations = 10L, silent = TRUE),
    "Error in mcmc, check $error_debug",
    fixed = TRUE
  )
  mcmc <- dj$new(
    data = data_list,
    df_params = df_params,
    loglike = r_loglike,
    logprior = r_logprior3
  )
  
  expect_error(
    mcmc$burn(iterations = 10L, silent = TRUE),
    "Error in mcmc, check $error_debug",
    fixed = TRUE
  )
  
})
