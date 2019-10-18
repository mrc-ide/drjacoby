context("test-mcmc-with-cpp-likelihood-and-prior")

test_that("Cpp likelihood and prior", {
  set.seed(1)
  # define true parameter values
  mu_true <- 3
  sigma_true <- 2
  # draw example data
  x <- rnorm(1000, mean = mu_true, sd = sigma_true)
  # define parameters dataframe
  df_params <- data.frame(name = c("mu", "sigma"),
                          min = c(-10, 0),
                          max = c(10, Inf),
                          init = c(5, 1))
  # Null log likelihood
  cpp_loglike_null <- "SEXP loglike(std::vector<double> params, std::vector<double> x) {
  double out = 0.0;
  return Rcpp::wrap(out);
    }"
  # Log likelihood
  cpp_loglike <- "SEXP loglike(std::vector<double> params, std::vector<double> x) {
  double ret = 0.0;
  for (int i = 0; i < x.size(); ++i) {
  ret += R::dnorm(x[i], params[0], params[1], 1);
  }
  return Rcpp::wrap(ret);
    }"
  # Log prior
  cpp_logprior_strong <- "SEXP logprior(std::vector<double> params){
  
  // calculate logprior
  double ret = R::dnorm(params[0], 6, 0.1, 1) +
      R::dnorm(params[1], 1, 0.1, 1);
  
  // return as SEXP
  return Rcpp::wrap(ret);
    }"

  # Null log prior
  cpp_logprior_null <- "SEXP loglike(std::vector<double> params) {
  double out = 0.0;
  return Rcpp::wrap(out);
    }"
  
  # Compilation checks
  expect_null(check_likelihood_compilation(cpp_loglike))
  expect_null(check_prior_compilation(cpp_logprior_strong))
  expect_error(check_likelihood_compilation(1:2))
  expect_error(check_prior_compilation(1:2))
  
  cpp_mcmc_null <- run_mcmc(data = x,
                            df_params = df_params,
                            loglike = cpp_loglike_null,
                            logprior = cpp_logprior_strong,
                            burnin = 1e3,
                            samples = 1e3,
                            silent = TRUE)
  
  posterior_estimate <- apply(cpp_mcmc_null$chain1$theta_sampling$rung1, 2, median)
  expect_lt(posterior_estimate[1] - 6, 0.1)
  expect_lt(posterior_estimate[2] - 1, 0.1)
  
  cpp_mcmc_data <- run_mcmc(data = x,
                            df_params = df_params,
                            loglike = cpp_loglike,
                            logprior = cpp_logprior_null,
                            burnin = 1e3,
                            samples = 1e4,
                            silent = TRUE)
  posterior_estimate2 <- apply(cpp_mcmc_data$chain1$theta_sampling$rung1, 2, median)
  expect_lt(posterior_estimate2[1] - 3, 0.1)
  expect_lt(posterior_estimate2[2] - 2, 0.1)
})
