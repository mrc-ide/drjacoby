
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
df_params <- define_params(name = "mu", min = -10, max = 10,
                           name = "sigma", min = 0, max = Inf)


mcmc <- run_mcmc(data = list(x = x),
                 df_params = df_params,
                 loglike = r_loglike,
                 logprior = r_logprior,
                 burnin = 1e3,
                 samples = 1e3)

plot_trace(mcmc, show = "mu", phase = "burnin")
