#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
// drjacoby_function_start
SEXP loglike(std::vector<double> params, std::vector<double> x) {
  
  // extract parameters
  double mu = params[0];
  double sigma = params[1];
  
  // sum log-likelihood over all data
  double ret = 0.0;
  for (unsigned int i = 0; i < x.size(); ++i) {
    ret += R::dnorm(x[i], mu, sigma, true);
  }
  
  // return as SEXP
  return Rcpp::wrap(ret);
}

