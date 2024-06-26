
#------------------------------------------------
test_that("plots do not produce errors", {
  
  # define example data
  # (use small series so that loglikelihood comes back +ve in some cases, which
  # forces down a particular path in some plotting functions)
  data_list <- list(x = (1:10)*1e-2)
  
  # define parameters dataframe
  df_params <- rbind.data.frame(list("mu", -10, 10, 5),
                                list("sigma", 0, 20, 1))
  names(df_params) <- c("name", "min", "max", "init")
  
  # log likelihood
  r_loglike <- function(params, data, misc) {
    sum(dnorm(data$x, mean = params["mu"], sd = params["sigma"], log = TRUE))
  }
  
  # log prior
  r_logprior <- function(params, misc) {
    dnorm(params["mu"], log = TRUE) + dlnorm(params["sigma"])
  }
  
  # run MCMC
  mcmc_out <- run_mcmc(data = data_list,
                       df_params = df_params,
                       loglike = r_loglike,
                       logprior = r_logprior,
                       burnin = 1e2,
                       samples = 3e3,
                       rungs = 2,
                       chains = 2,
                       silent = TRUE)
  
  # expect no output (messages or warnings) in all standard plotting functions
  expect_silent(plot_autocorrelation(mcmc_out))
  
  expect_silent(plot_mc_acceptance(mcmc_out))
  expect_silent(plot_mc_acceptance(mcmc_out, chain = 1))
  expect_silent(plot_mc_acceptance(mcmc_out, x_axis_type = 2))
  
  expect_silent(plot_trace(mcmc_out, display = FALSE))
  expect_silent(plot_trace(mcmc_out, show = "mu", display = FALSE))
  expect_silent(plot_trace(mcmc_out, hide = "mu", display = FALSE))
  expect_silent(plot_trace(mcmc_out, phase = "both", display = FALSE))
  
  expect_silent(plot_scatter(mcmc_out, parameter1 = "mu", parameter2 = "sigma"))
  expect_silent(plot_scatter(mcmc_out, parameter1 = "mu", parameter2 = "sigma", phase = "both"))
  
  expect_silent(plot_cor_mat(mcmc_out))
  expect_silent(plot_cor_mat(mcmc_out, show = c("mu", "sigma")))
  expect_silent(plot_cor_mat(mcmc_out, phase = "both"))
  
  expect_silent(plot_credible(mcmc_out))
  expect_silent(plot_credible(mcmc_out, show = "mu"))
  expect_silent(plot_credible(mcmc_out, phase = "both"))
  
  expect_silent(plot_rung_loglike(mcmc_out))
  expect_silent(plot_rung_loglike(mcmc_out, x_axis_type = 2))
  expect_silent(plot_rung_loglike(mcmc_out, y_axis_type = 2))
  expect_silent(plot_rung_loglike(mcmc_out, y_axis_type = 3))
  
  # repeat run with single rung
  mcmc_out <- run_mcmc(data = data_list,
                       df_params = df_params,
                       loglike = r_loglike,
                       logprior = r_logprior,
                       burnin = 1e2,
                       samples = 3e3,
                       rungs = 1,
                       chains = 2,
                       silent = TRUE)
  
  # further tests
  expect_error(plot_mc_acceptance(mcmc_out))
  
})
