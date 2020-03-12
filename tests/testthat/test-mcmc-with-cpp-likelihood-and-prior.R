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
                            chains = 1,
                            silent = TRUE)
  
  pe <- dplyr::filter(cpp_mcmc_null$output, stage == "sampling") %>%
    dplyr::select(mu, sigma)
  posterior_estimate <- apply(pe, 2, median)
  expect_lt(posterior_estimate[1] - 6, 0.1)
  expect_lt(posterior_estimate[2] - 1, 0.1)
  
  cpp_mcmc_data <- run_mcmc(data = x,
                            df_params = df_params,
                            loglike = cpp_loglike,
                            logprior = cpp_logprior_null,
                            burnin = 1e3,
                            samples = 1e4,
                            silent = TRUE)
  pe <- dplyr::filter(cpp_mcmc_data$output, stage == "sampling") %>%
    dplyr::select(mu, sigma)
  posterior_estimate2 <- apply(pe, 2, median)
  expect_lt(posterior_estimate2[1] - 3, 0.1)
  expect_lt(posterior_estimate2[2] - 2, 0.1)
  
  ## Multiple chains
  cpp_mcmc_chains <- run_mcmc(data = x,
                            df_params = df_params,
                            loglike = cpp_loglike,
                            logprior = cpp_logprior_null,
                            chains = 2,
                            burnin = 1e3,
                            samples = 1e4,
                            silent = TRUE)
  expect_length(cpp_mcmc_chains, 3)
  pe <- dplyr::filter(cpp_mcmc_chains$output, stage == "sampling", chain == "chain1") %>%
    dplyr::select(mu, sigma)
  posterior_estimate3a <- apply(pe, 2, median)
  expect_lt(posterior_estimate3a[1] - 3, 0.1)
  expect_lt(posterior_estimate3a[2] - 2, 0.1)
  pe <- dplyr::filter(cpp_mcmc_chains$output, stage == "sampling", chain == "chain2") %>%
    dplyr::select(mu, sigma)
  posterior_estimate3b <- apply(pe, 2, median)
  expect_lt(posterior_estimate3b[1] - 3, 0.1)
  expect_lt(posterior_estimate3b[2] - 2, 0.1)
  
  ## Metropolis coupling
  mcmc_out_MC <- run_mcmc(data = x,
                       df_params = df_params,
                       loglike = cpp_loglike,
                       logprior = cpp_logprior_null,
                       burnin = 1e3,
                       samples = 1e4,
                       rungs = 4,
                       silent = TRUE)
  pe <- dplyr::filter(mcmc_out_MC$output, stage == "sampling", chain == "chain1", rung == "rung1") %>%
    dplyr::select(mu, sigma)
  posterior_estimate4 <- apply(pe, 2, median)
  expect_lt(posterior_estimate4[1] - 3, 0.1)
  expect_lt(posterior_estimate4[2] - 2, 0.1)
})
