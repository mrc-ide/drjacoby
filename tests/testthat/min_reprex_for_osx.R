context("test-mcmc-on-osx-with-clang-Xcode")

test_that("R MCMC likelihood runs", {
  set.seed(1)
  
  #...................... 
  # LL and LP
  #......................
  logpriorfunc <- "SEXP logprior(Rcpp::NumericVector params, int param_i, Rcpp::List misc) { double ret = -1.0; return Rcpp::wrap(ret);}"
  loglikfunc <- "SEXP loglike(Rcpp::NumericVector params, int param_i, Rcpp::List data, Rcpp::List misc) {
              double cfr = 0.5;
              double mod = 18.8;
              double sod = -0.45;

              std::vector<int> death = Rcpp::as< std::vector<int> >(data[\"death_interval\"]);
              double loglik = 0.0;
              for (int j = 0; j < death.size(); j++) {
                double l1 = 0.6;
                double l2 = 0.5;
                loglik += log(cfr) + log(l1 - l2);
              }

              if (!std::isfinite(loglik)) {
                const double OVERFLO_DOUBLE = DBL_MAX/100.0;
                loglik = -OVERFLO_DOUBLE;
              }
              return Rcpp::wrap(loglik);
            }"
  
  #...................... 
  # set up and run Dr. Jacoby
  #......................
  dth_intrvl <- round(runif(2e3, 1, 50), 0)
  data_list <- list("death_interval" = dth_intrvl)
  
  # dr jacoby
  df_params <- rbind.data.frame(list("mod", 2, 30, 17),
                                list("sod", 0, 3, 0.6),
                                list("cfr", 0, 1, 0.5))
  names(df_params) <- c("name", "min", "max", "init")
  misc_list <- list()
  
  mcmcout <- drjacoby::run_mcmc(data = data_list,
                                df_params = df_params,
                                misc = misc_list,
                                loglike = loglikfunc,
                                logprior = logpriorfunc,
                                burnin = 1e3,
                                samples = 1e3,
                                chains = 3,
                                rungs = 1,
                                silent = T)
  
  # run multiple seeds to find a mem fault
  wrap_seedling <- function(x) {
    cat(c(x, "\n"))
    mcmcout <- drjacoby::run_mcmc(data = data_list,
                                  df_params = df_params,
                                  misc = misc_list,
                                  loglike = loglikfunc,
                                  logprior = logpriorfunc,
                                  burnin = 1e3,
                                  samples = 1e3,
                                  chains = 3)
    return(mcmcout)
  }
  
  seeds <- round(runif(1e2, 10, 5e3), 0)
  testwrap <- lapply(seeds, wrap_seedling)
  expect_equal(length(testwrap), 100)
  
})
