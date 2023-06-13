#include "cpp11.hpp"
#include "cpp11/doubles.hpp"
#include "Rmath.h"
using namespace cpp11;

[[cpp11::register]]
double loglike_cpp11(const doubles params, const list data, const list misc){
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
  
  return ret;
}

[[cpp11::register]]
double logprior_cpp11(const doubles params, const list misc){
  // extract parameters
  double sigma = params["sigma"];
  
  // calculate logprior
  double ret = -log(20.0) + dlnorm(sigma, 0.0, 1.0, true);
  
  return ret;
}
