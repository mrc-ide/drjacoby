devtools::load_all()

# Source our CPP file
Rcpp::sourceCpp("test_scripts/test_cpp.cpp")

# define true parameter values
mu_true <- 3
sigma_true <- 2

# draw example data
data_list <- list(x = rnorm(10, mean = mu_true, sd = sigma_true))
# define parameters dataframe
df_params <- define_params(name = "mu", min = -10, max = 10,
                           name = "sigma", min = 0, max = Inf)

# define log-prior function
r_logprior <- function(params, misc) {
  
  # extract parameter values
  mu <- params["mu"]
  sigma <- params["sigma"]
  
  # calculate log-prior
  ret <- dunif(mu, min = -10, max = 10, log = TRUE) +
    dlnorm(sigma, meanlog = 0, sdlog = 1.0, log = TRUE)
  
  # return
  return(ret)
}
#Run
mcmc <- run_mcmc(data = data_list,
                 df_params = df_params,
                 loglike = "loglike",
                 logprior = r_logprior,
                 burnin = 1e3,
                 samples = 1e3,
                 pb_markdown = TRUE)

plot_par(mcmc, show = "mu", phase = "sampling")
