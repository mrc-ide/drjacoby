#include "cpp11.hpp"
#include "cpp11/doubles.hpp"
#include "Rmath.h"
using namespace cpp11;

[[cpp11::register]]
double logprior_normal_cpp11(const doubles params, const list misc) {
  
  // unpack params
  double mu = params["mu"];
  double sigma = params["sigma"];
  
  // calculate logprior
  double ret = Rf_dnorm4(mu, 6, 0.1, 1) + Rf_dnorm4(sigma, 1, 0.1, 1);
  
  // return as SEXP
  return ret;
}
