#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
SEXP loglike(Rcpp::NumericVector params, Rcpp::List data, Rcpp::List misc) {
  
  // unpack parameters
  int pi = 0;
  std::vector<double> mu(6);
  for (int i = 0; i < 6; ++i) {
    mu[i] = params[pi++];
  }
  double sigma = params["sigma"];
  double phi = params["phi"];
  double tau = params["tau"];
  
  // get current update block
  int block = misc["block"];
  
  // distinct method for first 6 blocks, vs 7th
  double ret = 0.0;
  if (block == 7) {  // likelihood for the global parameters
    
    // calculate likelihood component
    for (int i = 0; i < 6; ++i) {
      ret += R::dnorm(mu[i], phi, tau, true);
    }
    
  } else {  // likelihood for each of the 6 data groups
    
    // get data for this group
    NumericVector x = data[block - 1];
    
    // calculate likelihood component
    for (int i = 0; i < x.size(); ++i) {
      ret += R::dnorm(x[i], mu[block - 1], sigma, true);
    }
    
  }
  
  return Rcpp::wrap(ret);
}

// [[Rcpp::export]]
SEXP logprior(Rcpp::NumericVector params, Rcpp::List misc) {
  
  // unpack parameters
  double sigma = params["sigma"];
  double phi = params["phi"];
  double tau = params["tau"];
  
  // apply priors
  double ret = R::dgamma(sigma, 0.01, 100.0, true) +
    R::dnorm(phi, 0, 1000, true) +
    R::dgamma(tau, 0.01, 100.0, true);
  
  return Rcpp::wrap(ret);
}

// [[Rcpp::export]]  
SEXP create_xptr(std::string function_name) {  
  typedef SEXP (*funcPtr_likelihood)(Rcpp::NumericVector params, Rcpp::List data, Rcpp::List misc);  
  typedef SEXP (*funcPtr_prior)(Rcpp::NumericVector params, Rcpp::List misc);  
  
  if (function_name == "loglike"){
    return(Rcpp::XPtr<funcPtr_likelihood>(new funcPtr_likelihood(&loglike)));
  } 
  if (function_name == "logprior"){
    return(Rcpp::XPtr<funcPtr_prior>(new funcPtr_prior(&logprior)));
  } 
  
  stop("cpp function %i not found", function_name);
}
