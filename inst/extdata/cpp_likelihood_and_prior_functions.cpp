#include "cpp11.hpp"
#include "cpp11/doubles.hpp"
#include "Rmath.h"
// #include <vector>
using namespace cpp11;
// namespace writable = cpp11::writable;

[[cpp11::register]]
double loglike_normal_cpp11(const doubles params, const list data, const list misc){
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
  
  return(ret);
}

[[cpp11::register]]
double logprior_null_cpp11(const doubles params, const list misc){
  double ret = 0;
  return(ret);
}

[[cpp11::register]]
double loglike_null_cpp11(const doubles params, const list data, const list misc){
  double ret = 0;
  return(ret);
}
