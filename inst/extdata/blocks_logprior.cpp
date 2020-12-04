#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
// drjacoby_function_start
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
// drjacoby_function_end
