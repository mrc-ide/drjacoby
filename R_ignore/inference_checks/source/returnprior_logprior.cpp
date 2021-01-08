#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
// drjacoby_function_start
SEXP logprior(Rcpp::NumericVector params, Rcpp::List misc) {
  
  double ret = R::dnorm(params["real_line"], 0.0, 1.0, true) +
    R::dgamma(-params["neg_line"], 5.0, 5.0, true) +
    R::dgamma(params["pos_line"], 5.0, 5.0, true) +
    R::dbeta(params["unit_interval"], 3.0, 3.0, true);
  
  return Rcpp::wrap(ret);
}
// drjacoby_function_end
