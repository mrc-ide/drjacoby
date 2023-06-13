#include "cpp11.hpp"
#include "cpp11/doubles.hpp"
#include "Rmath.h"
using namespace cpp11;
namespace writable = cpp11::writable;

[[cpp11::register]]
double block_loglike_cpp11(const doubles params, const list data, const list misc){
  
  // unpack parameters
  int pi = 0;
  writable::doubles mu(6);
  for (int i = 0; i < 6; ++i) {
    mu[i] = params[pi++];
  }
  double sigma = params["sigma"];
  double phi = params["phi"];
  double tau = params["tau"];
  
  // get current update block
  int block = cpp11::as_cpp<int>(misc["block"]);
  
  // distinct method for first 6 blocks, vs 7th
  double ret = 0.0;
  if (block == 7) {  // likelihood for the global parameters
    
    // calculate likelihood component
    for (int i = 0; i < 6; ++i) {
      ret += Rf_dnorm4(mu[i], phi, tau, true);
    }
    
  } else {  // likelihood for each of the 6 data groups
    
    // get data for this group
    doubles x = data[block - 1];
    
    // calculate likelihood component
    for (int i = 0; i < x.size(); ++i) {
      ret += Rf_dnorm4(x[i], mu[block - 1], sigma, true);
    }
    
  }
  
  
  return ret;
}

[[cpp11::register]]
double block_logprior_cpp11(const doubles params, const list misc){
  // unpack parameters
  double sigma = params["sigma"];
  double phi = params["phi"];
  double tau = params["tau"];
  
  // apply priors
  double ret = dgamma(sigma, 0.01, 100.0, true) +
    Rf_dnorm4(phi, 0, 1000, true) +
    dgamma(tau, 0.01, 100.0, true);
  
  return ret;
}

