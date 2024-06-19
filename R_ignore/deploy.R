
# deploy.R
#
# Author: Bob Verity
# Date: 2024-06-19
#
# Purpose:
# Test package functions.

# ------------------------------------------------------------------

library(tidyverse)


# define R loglike function
r_loglike <- function(params, data, misc) {
  sum(dnorm(data$x, mean = params["mu"], sd = params["sigma"], log = TRUE))
}

# define R logprior function
r_logprior <- function(params, misc) {
  dunif(params["mu"], -10, 10, log = TRUE) + dnorm(params["sigma"], 0, 1, log = TRUE)
}

# source C++ likelihood and prior
#cpp11::cpp_source("ignore/source/normal_model.cpp")


# ------------------------------------------------------------------

set.seed(1)

# define data
x <- rnorm(1e2, mean = 1, sd = 4)

# define parameters dataframe
df_params <- data.frame() %>%
  add_parameter(name = "mu", min = -10, max = 10) %>%
  add_parameter(name = "sigma", min = 0, max = Inf)

my_mcmc <- dj$new(data = list(x = x),
                  df_params = df_params,
                  loglikelihood = loglike_cpp11,
                  logprior = logprior_cpp11)

my_mcmc$tune(iterations = 1e3, target_rung_acceptance = 0.25)
my_mcmc

my_mcmc$plot_tuning_rejection_rate()
my_mcmc$plot_local_communication_barrier()

my_mcmc$burn(1e3)
my_mcmc$sample(1e3)
my_mcmc$plot_mc_acceptance_rate(beta_axis = FALSE)

my_mcmc$plot_par("mu")
my_mcmc$plot_par("sigma")

z <- my_mcmc$output()

