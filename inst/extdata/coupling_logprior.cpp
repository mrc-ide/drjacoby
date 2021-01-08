#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
// drjacoby_function_start
SEXP logprior(Rcpp::NumericVector params, Rcpp::List misc) {
  
  // unpack parameters
  double epsilon = params["epsilon"];
  
  // calculate logprior
  double ret = -log(20.0) - log(10.0) + R::dnorm(epsilon, 0.0, 1.0, true);
  
  // return as SEXP
  return Rcpp::wrap(ret);
}
// drjacoby_function_end
