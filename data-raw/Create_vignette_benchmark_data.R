# Create vignette benchmarking data
devtools::load_all()
# Set up the example from the "Example MCMC Implementation" vignette
# set random seed
set.seed(1)
# define true parameter values
mu_true <- 3
sigma_true <- 2
# draw example data
x <- rnorm(10, mean = mu_true, sd = sigma_true)
# define parameters dataframe
df_params <- data.frame(name = c("mu", "sigma"),
                        min = c(-10, 0),
                        max = c(10, Inf),
                        init = c(5, 1))

# R log likelihood and prior functions
r_loglike <- function(params, x) {
  mu <- params[1]
  sigma <- params[2]
  ret <- sum(dnorm(x, mean = mu, sd = sigma, log = TRUE))
  return(ret)
}
r_logprior <- function(params) {
  mu <- params[1]
  sigma <- params[2]
  ret <- dunif(mu, min = -10, max = 10, log = TRUE) + dlnorm(sigma, meanlog = 0, sdlog = 1.0, log = TRUE)
  return(ret)
}
# Rcpp log likelihood and prior functions
cpp_loglike <- "SEXP loglike(std::vector<double> params, std::vector<double> x) {
  double mu = params[0];
  double sigma = params[1];
  double ret = 0.0;
  for (size_t i = 0; i < sizeof(x); ++i) {
    ret += -0.5*log(2*M_PI*sigma*sigma) - (x[i] - mu)*(x[i] - mu)/(2*sigma*sigma);
  }
  return Rcpp::wrap(ret);
}"
cpp_logprior <- "SEXP logprior(std::vector<double> params){
  double sigma = params[1];
  double ret = -log(20.0) - log(sigma) - 0.5*log(2*M_PI*1.0*1.0) - (log(sigma) - 0)*(log(sigma) - 0)/(2*1.0*1.0);
  return Rcpp::wrap(ret);
}"

## Run benchmarking
bm <- microbenchmark::microbenchmark(r_version = run_mcmc(data = x,
                                                          df_params = df_params,
                                                          loglike = r_loglike,
                                                          logprior = r_logprior,
                                                          burnin = 1e3,
                                                          samples = 1e5,
                                                          silent = TRUE),
                                     cpp_version = run_mcmc(data = x,
                                                            df_params = df_params,
                                                            loglike = cpp_loglike,
                                                            logprior = cpp_logprior,
                                                            burnin = 1e3,
                                                            samples = 1e5,
                                                            silent = TRUE),
                                     mixed_version = run_mcmc(data = x,
                                                              df_params = df_params,
                                                              loglike = cpp_loglike,
                                                              logprior = r_logprior,
                                                              burnin = 1e3,
                                                              samples = 1e5,
                                                              silent = TRUE),
                                     times = 10)
summary(bm)
bm$expr <- forcats::fct_recode(bm$expr,
                               "R" = "r_version",
                               "Mixed" = "mixed_version",
                               "Cpp" = "cpp_version") %>%
  forcats::fct_relevel(c("R", "Mixed", "Cpp"))
# Save the benchmarking data internally (Not exported)
usethis::use_data(bm, internal = TRUE, overwrite = TRUE)
