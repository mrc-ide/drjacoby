#include "cpp11.hpp"
#include "cpp11/doubles.hpp"
#include "Rmath.h"
// #include <vector>
using namespace cpp11;
// namespace writable = cpp11::writable;

[[cpp11::register]]
double loglike(doubles params, list data, list misc){
  // extract parameters
  double mu = params["mu"];
  double sigma = params["sigma"];
  
  // unpack data
  doubles x = data["x"];
  
  // sum log-likelihood over all data
  double ret = 0.0;
  for (int i = 0; i < x.size(); ++i) {
    ret += Rf_dnorm4(x[i], mu, sigma, 1);
  }
  
  // return as SEXP
  return(ret);
}

[[cpp11::register]]
double logprior(doubles params, list misc){
  
  // extract parameters
  double sigma = params["sigma"];
  
  // calculate logprior
  double ret = -log(20.0) + Rf_dnorm4(sigma, 0.0, 1.0, true);
  
  // return as SEXP
  return(ret);
}
