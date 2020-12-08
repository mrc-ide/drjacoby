#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
// drjacoby_function_start
SEXP loglike(Rcpp::NumericVector params, Rcpp::List data, Rcpp::List misc) {
  double mu = params["mu"];
  double gamma = params["gamma"];
  double ret = -gamma*(mu*mu - 1.0)*(mu*mu - 1.0);
  
  return Rcpp::wrap(ret);
}
// drjacoby_function_end
