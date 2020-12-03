context("test-mcmc-with-cpp-likelihood-and-prior")

test_that("Cpp likelihood and prior", {
  set.seed(1)
  
  # define true parameter values
  mu_true <- 3
  sigma_true <- 2
  
  # draw example data
  x <- rnorm(1000, mean = mu_true, sd = sigma_true)
  data_list <- list(x = x)
  
  # define parameters dataframe
  df_params <- data.frame(name = c("mu", "sigma"),
                          min = c(-10, 0),
                          max = c(10, 10),
                          init = c(5, 1))
  
  # Null log likelihood
  cpp_loglike_null <- "SEXP loglike(Rcpp::NumericVector params, int param_i, Rcpp::List data, Rcpp::List misc) {
    return Rcpp::wrap(0.0);
  }"
  
  # Log likelihood
  cpp_loglike <- "SEXP loglike(Rcpp::NumericVector params, int param_i, Rcpp::List data, Rcpp::List misc) {
    
    // unpack data
    std::vector<double> x = Rcpp::as< std::vector<double> >(data[\"x\"]);
    
    // unpack params
    double mu = params[\"mu\"];
    double sigma = params[\"sigma\"];
    
    // calculate log-likelihood
    double ret = 0.0;
    for (unsigned int i = 0; i < x.size(); ++i) {
      ret += R::dnorm(x[i], mu, sigma, 1);
    }
    
    // catch underflow
    if (!std::isfinite(ret)) {
      const double OVERFLO_DOUBLE = DBL_MAX/100.0;
      ret = -OVERFLO_DOUBLE;
    }
    
    return Rcpp::wrap(ret);
  }"
  
  # Log prior
  cpp_logprior_strong <- "SEXP logprior(Rcpp::NumericVector params, int param_i, Rcpp::List misc) {
    
    // unpack params
    double mu = params[\"mu\"];
    double sigma = params[\"sigma\"];
    
    // calculate logprior
    double ret = R::dnorm(mu, 6, 0.1, 1) + R::dnorm(sigma, 1, 0.1, 1);
    
    // catch underflow
    if (!std::isfinite(ret)) {
      const double OVERFLO_DOUBLE = DBL_MAX/100.0;
      ret = -OVERFLO_DOUBLE;
    }
    
    // return as SEXP
    return Rcpp::wrap(ret);
  }"

  # Null log prior
  cpp_logprior_null <- "SEXP logprior(Rcpp::NumericVector params, int param_i, Rcpp::List misc) {
    return Rcpp::wrap(0.0);
  }"
  
  # Compilation checks
  expect_null(check_likelihood_compilation(cpp_loglike))
  expect_null(check_prior_compilation(cpp_logprior_strong))
  expect_error(check_likelihood_compilation(1:2))
  expect_error(check_prior_compilation(1:2))
  
  # run MCMC
  cpp_mcmc_null <- run_mcmc(data = data_list,
                            df_params = df_params,
                            loglike = cpp_loglike_null,
                            logprior = cpp_logprior_strong,
                            burnin = 1e3,
                            samples = 1e3,
                            chains = 1,
                            silent = TRUE)
  
  # subset output
  pe <- dplyr::filter(cpp_mcmc_null$output, stage == "sampling") %>%
    dplyr::select(mu, sigma)
  
  # check posterior estimates
  posterior_estimate <- apply(pe, 2, median)
  expect_lt(posterior_estimate["mu"] - 6, 0.1) # should approximately follow the values in logprior_strong
  expect_lt(posterior_estimate["sigma"] - 1, 0.1) # should approximately follow the values in logprior_strong
  
  # run MCMC with null prior
  cpp_mcmc_data <- run_mcmc(data = data_list,
                            df_params = df_params,
                            loglike = cpp_loglike,
                            logprior = cpp_logprior_null,
                            burnin = 1e3,
                            samples = 1e4,
                            silent = TRUE)
  
  # subset output
  pe <- dplyr::filter(cpp_mcmc_data$output, stage == "sampling") %>%
    dplyr::select(mu, sigma)
  
  # check posterior estimates
  posterior_estimate2 <- apply(pe, 2, median)
  expect_lt(posterior_estimate2["mu"] - mu_true, 0.1)
  expect_lt(posterior_estimate2["sigma"] - sigma_true, 0.1)
  
  ## Multiple chains
  cpp_mcmc_chains <- run_mcmc(data = data_list,
                              df_params = df_params,
                              loglike = cpp_loglike,
                              logprior = cpp_logprior_null,
                              chains = 2,
                              burnin = 1e3,
                              samples = 1e4,
                              silent = TRUE)
  expect_length(cpp_mcmc_chains, 3)
  
  # subset output to chain1 and check posterior estimates
  pe <- dplyr::filter(cpp_mcmc_chains$output, stage == "sampling", chain == "chain1") %>%
    dplyr::select(mu, sigma)
  posterior_estimate3a <- apply(pe, 2, median)
  expect_lt(posterior_estimate3a["mu"] - mu_true, 0.1)
  expect_lt(posterior_estimate3a["sigma"] - sigma_true, 0.1)
  
  # subset output to chain2 and check posterior estimates
  pe <- dplyr::filter(cpp_mcmc_chains$output, stage == "sampling", chain == "chain2") %>%
    dplyr::select(mu, sigma)
  posterior_estimate3b <- apply(pe, 2, median)
  expect_lt(posterior_estimate3b["mu"] - mu_true, 0.1)
  expect_lt(posterior_estimate3b["sigma"] - sigma_true, 0.1)
  
  ## Metropolis coupling
  mcmc_out_MC <- run_mcmc(data = data_list,
                          df_params = df_params,
                          loglike = cpp_loglike,
                          logprior = cpp_logprior_null,
                          burnin = 1e3,
                          samples = 1e4,
                          rungs = 4,
                          silent = TRUE)
  
  # subset output
  pe <- dplyr::filter(mcmc_out_MC$output, stage == "sampling", chain == "chain1") %>%
    dplyr::select(mu, sigma)
  
  # check posterior estimates
  posterior_estimate4 <- apply(pe, 2, median)
  expect_lt(posterior_estimate4["mu"] - mu_true, 0.1)
  expect_lt(posterior_estimate4["sigma"] - sigma_true, 0.1)
})
