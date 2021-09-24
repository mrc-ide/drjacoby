#------------------------------------------------
test_that("Checks on likelihood and prior errors work", {
  
  # define true parameter values
  mu_true <- 3
  sigma_true <- 2
  
  # draw example data
  data_list <- list(x = 1)
  
  # define parameters dataframe
  df_params <- define_params(name = "mu", min = -10, max = 10, init = 1)
  
  
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
  
  expect_error(run_mcmc(data = data_list,
                        df_params = df_params,
                        loglike = r_loglike,
                        logprior = r_logprior1,
                        burnin = 1e3,
                        samples = 1e3,
                        chains = 1))
  expect_error(run_mcmc(data = data_list,
                        df_params = df_params,
                        loglike = r_loglike,
                        logprior = r_logprior2,
                        burnin = 1e3,
                        samples = 1e3,
                        chains = 1))
  expect_error(run_mcmc(data = data_list,
                        df_params = df_params,
                        loglike = r_loglike,
                        logprior = r_logprior3,
                        burnin = 1e3,
                        samples = 1e3,
                        chains = 1))
  
})

#------------------------------------------------
test_that("Negative Inf initialisation checks", {
  # define true parameter values
  mu_true <- 3
  sigma_true <- 2
  
  # draw example data
  data_list <- list(x = 1)
  
  # define parameters dataframe
  df_params <- define_params(name = "mu", min = -10, max = 10, init = 1)
  
  
  # Initialisation fails with initial params
  r_loglike <- function(params, data, misc) {
    if(params[[1]] == 1){
      -Inf
    } else {
      -1
    }
  }
  
  # define log-prior function
  r_logprior <- function(params, misc) {
    0
  }
  
  expect_error(run_mcmc(data = data_list,
                        df_params = df_params,
                        loglike = r_loglike,
                        logprior = r_logprior,
                        burnin = 1e3,
                        samples = 1e3,
                        chains = 1))
  
  # Define log likelihood
  r_loglike <- function(params, data, misc) {
    0
  }
  
  # Initialisation fails with initial params
  r_logprior <- function(params, misc) {
    if(params[[1]] == 1){
      -Inf
    } else {
      -1
    }
  }
  
  expect_error(run_mcmc(data = data_list,
                        df_params = df_params,
                        loglike = r_loglike,
                        logprior = r_logprior,
                        burnin = 1e3,
                        samples = 1e3,
                        chains = 1))
  
})

#------------------------------------------------
test_that("-Inf tolerated in likelihood and prior", {
  
  # define parameter mu with [0,1] range
  df_params <- define_params(name = "mu", min = 0, max = 1, init = 0.1)
  
  # likelihood is -Inf for half the domain
  r_loglike <- function(params, data, misc) {
    mu <- params["mu"]
    if (mu > 0.5) {
      return(-Inf)
    } else {
      return(2)
    }
  }
  
  # prior is finite everywhere
  r_logprior <- function(params, misc) {
    0
  }
  
  # run mcmc with two rungs at 0 and 1
  mcmc <- run_mcmc(data = list(x = 1),
                   df_params = df_params,
                   loglike = r_loglike,
                   logprior = r_logprior,
                   burnin = 1e3,
                   samples = 1e3,
                   chains = 1,
                   beta_manual = c(0, 1))
  
  # get maximum value of mu for each rung. Should be >0.5 for rung1 and <0.5 for
  # rung2
  mu_max_df <- mcmc$output %>%
    dplyr::group_by(rung) %>%
    dplyr::summarise(mu = max(mu))
  
  expect_equal((mu_max_df$mu > 0.5), c(TRUE, FALSE))
  
  # rung1 should contain -Inf loglike values whenever mu > 0.5, and finite values
  # otherwise
  rung1 <- mcmc$output %>%
    dplyr::filter(rung == 1)
  
  expect_true(all(rung1$loglikelihood[rung1$mu > 0.5] == -Inf))
  expect_true(all(is.finite(rung1$loglikelihood[rung1$mu < 0.5])))
  
  # expect no errors when plotting loglikelihoods
  expect_error(plot_rung_loglike(mcmc, y_axis_type = 1), NA)
  expect_error(plot_rung_loglike(mcmc, y_axis_type = 2), NA)
  expect_error(plot_rung_loglike(mcmc, y_axis_type = 3), NA)
  
  
  ################################
  
  
  # same as above, but with -Inf implemented in prior not likelihood
  r_loglike <- function(params, data, misc) {
    0
  }
  
  # prior is finite everywhere
  r_logprior <- function(params, misc) {
    mu <- params["mu"]
    if (mu > 0.5) {
      return(-Inf)
    } else {
      return(2)
    }
  }
  
  # run mcmc with two rungs at 0 and 1
  mcmc <- run_mcmc(data = list(x = 1),
                   df_params = df_params,
                   loglike = r_loglike,
                   logprior = r_logprior,
                   burnin = 1e3,
                   samples = 1e3,
                   chains = 1,
                   beta_manual = c(0, 1))
  
  # now all mu values should be in the range [0,0.5], irrespective of rung
  expect_true(all(mcmc$output$mu < 0.5))
  
})
