#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
// drjacoby_function_start
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
// drjacoby_function_end
