#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
SEXP loglike(Rcpp::NumericVector params, Rcpp::List data, Rcpp::List misc) {
  // Insert loglikelihood code here:
  
  double loglikelihood = 0.0;
  return Rcpp::wrap(loglikelihood);
}

// [[Rcpp::export]]
SEXP logprior(Rcpp::NumericVector params, Rcpp::List misc) {
  // Insert loglikelihood code here:
  
  double logprior = 0.0;
  return Rcpp::wrap(logprior);
}


// NOTE: Do not edit this function name
// [[Rcpp::export]]  
SEXP create_xptr(std::string function_name) {  
  typedef SEXP (*funcPtr_likelihood)(Rcpp::NumericVector params, Rcpp::List data, Rcpp::List misc);  
  typedef SEXP (*funcPtr_prior)(Rcpp::NumericVector params, Rcpp::List misc);  
  
  // NOTE: If your loglikelihood function is not called "loglike" please edit:
  if (function_name == "loglike"){
    return(Rcpp::XPtr<funcPtr_likelihood>(new funcPtr_likelihood(&loglike)));
  } 
  // NOTE: If your logprior function is not called "logprior" please edit:
  if (function_name == "logprior"){
    return(Rcpp::XPtr<funcPtr_prior>(new funcPtr_prior(&logprior)));
  } 
  
  stop("cpp function %i not found", function_name);
}

