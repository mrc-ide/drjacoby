#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
// drjacoby_function_start
SEXP loglike(Rcpp::NumericVector params, Rcpp::List data, Rcpp::List misc) {
  
  // unpack data
  std::vector<double> x = Rcpp::as< std::vector<double> >(data["x"]);
  
  // unpack parameters
  double alpha = params["alpha"];
  double beta = params["beta"];
  double epsilon = params["epsilon"];
  
  // sum log-likelihood over all data
  double mean = alpha*alpha*beta + epsilon;
  double ret = 0.0;
  for (unsigned int i = 0; i < x.size(); ++i) {
    ret += R::dnorm(x[i], mean, 1.0, true);
  }
  
  // return as SEXP
  return Rcpp::wrap(ret);
}
// drjacoby_function_end
