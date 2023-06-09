#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double loglike_normal_rcpp(Rcpp::NumericVector params, Rcpp::List data, Rcpp::List misc) {
  
  // extract parameters
  double mu = params["mu"];
  double sigma = params["sigma"];
  
  // unpack data
  std::vector<double> x = Rcpp::as< std::vector<double> >(data["x"]);
  
  // sum log-likelihood over all data
  double ret = 0.0;
  for (unsigned int i = 0; i < x.size(); ++i) {
    ret += R::dnorm(x[i], mu, sigma, true);
  }
  
  // return as SEXP
  return ret;
}

// [[Rcpp::export]]
double logprior_null_rcpp(Rcpp::NumericVector params, Rcpp::List misc) {
  
  // extract parameters
  double ret = 0.0;
  
  // return as SEXP
  return ret;
}
