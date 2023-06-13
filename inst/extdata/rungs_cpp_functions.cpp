#include "cpp11.hpp"
#include "cpp11/doubles.hpp"
#include "Rmath.h"
using namespace cpp11;

[[cpp11::register]]
double rung_loglike_cpp11(const doubles params, const list data, const list misc){
  // unpack data
  doubles x = data["x"];
  
  // unpack parameters
  double alpha = params["alpha"];
  double beta = params["beta"];
  double epsilon = params["epsilon"];
  
  // sum log-likelihood over all data
  double mean = alpha*alpha*beta + epsilon;
  double ret = 0.0;
  for (int i = 0; i < x.size(); ++i) {
    ret += Rf_dnorm4(x[i], mean, 1.0, true);
  }
  
  return ret;
}

[[cpp11::register]]
double rung_logprior_cpp11(const doubles params, const list misc){
  // unpack parameters
  double epsilon = params["epsilon"];
  
  // calculate logprior
  double ret = -log(20.0) - log(10.0) + Rf_dnorm4(epsilon, 0.0, 1.0, true);
  
  
  return ret;
}

