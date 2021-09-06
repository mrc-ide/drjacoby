#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
SEXP loglike(Rcpp::NumericVector params, Rcpp::List data, Rcpp::List misc) {
  
  int pi = 0;
  std::vector<double> mu(5);
  for (int i = 0; i < 5; ++i) {
    mu[i] = params[pi++];
  }
  
  int block = misc["block"];
  
  double ret = 0.0;
  if (block == 6) {  // likelihood for the global parameters
    for (int i = 0; i < 5; ++i) {
      ret += R::dnorm(mu[i], 0, 1.0, true);
    }
  } else {  // likelihood for each of the 6 data groups
    NumericVector x = data[block - 1];
    for (int i = 0; i < x.size(); ++i) {
      ret += R::dnorm(x[i], mu[block - 1], 1.0, true);
    }
  }
  
  return Rcpp::wrap(ret);
}

// [[Rcpp::export]]
SEXP logprior(Rcpp::NumericVector params, Rcpp::List misc) {
  return Rcpp::wrap(0.0);
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
