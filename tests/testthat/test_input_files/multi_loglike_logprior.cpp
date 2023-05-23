#include "cpp11.hpp"
#include "cpp11/doubles.hpp"
#include "Rmath.h"
using namespace cpp11;


// Silly example so we can show that our loglike function can call others internally
double none(){
  double x = 0.0;
  return x;
}

[[cpp11::register]]
double logprior_normal_cpp11_multi(const doubles params, const list data, const list misc) {
  // extract parameters
  double mu = params["mu"];
  double sigma = params["sigma"];
  
  // unpack data
  doubles x = data["x"];
  
  // sum log-likelihood over all data
  double ret = 0.0;
  for (int i = 0; i < x.size(); ++i) {
    ret += Rf_dnorm4(x[i], mu, sigma, 1) + none();
  }
  
  return ret;
}
