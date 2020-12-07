
#------------------------------------------------
test_that("plots do not produce errors", {
  
  # draw example data
  data_list <- list(x = rnorm(10))
  
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
                       samples = 1e2,
                       rungs = 2,
                       silent = TRUE)
  
  # expect no output (messages or warnings) in all standard plotting functions
  expect_silent(plot_autocorrelation(mcmc_out))
  expect_silent(plot_cor(mcmc_out, parameter1 = "mu", parameter2 = "sigma"))
  expect_silent(plot_cor_mat(mcmc_out))
  expect_silent(plot_credible(mcmc_out))
  expect_silent(plot_mc_acceptance(mcmc_out))
  expect_silent(plot_par(mcmc_out, display = FALSE))
  expect_silent(plot_rung_loglike(mcmc_out))
  
})
