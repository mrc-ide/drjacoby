
#pragma once

#include "System.h"
#include "misc_v4.h"
#include "probability.h"

#include <Rcpp.h>

//------------------------------------------------
// class defining MCMC particle
class Particle {
  
public:
  // PUBLIC OBJECTS
  
  // pointer to system object
  System * s_ptr;
  
  // local copies of some parameters for convenience
  int d;
  
  // beta_raised stores values of beta (the thermodynamic power), raised to the
  // power GTI_pow
  double beta_raised;
  
  // theta is the parameter vector in natural space
  std::vector<double> theta;
  std::vector<double> theta_prop;
  
  // adjustment factor to account for transformation
  double adj;
  
  // phi is a vector of transformed parameters
  std::vector<double> phi;
  std::vector<double> phi_prop;
  
  // proposal standard deviations
  double propSD;
  
  // likelihoods and priors
  double loglike;
  double loglike_prop;
  double logprior;
  double logprior_prop;
  
  // store acceptance rates
  int accept;
  
  
  // PUBLIC FUNCTIONS
  
  // constructors
  Particle() {};
  
  // initialise everything EXCEPT FOR likelihood and prior values
  void init(System &s, double beta_raised);
  
  // initialise likelihood and prior values
  template<class TYPE1, class TYPE2>
  void init_like(TYPE1 get_loglike, TYPE2 get_logprior) {
    loglike = Rcpp::as<double>(get_loglike(theta, s_ptr->x));
    logprior = Rcpp::as<double>(get_logprior(theta));
  }
  
  // update using flexible likelihood function
  template<class TYPE1, class TYPE2>
  void update(TYPE1 get_loglike, TYPE2 get_logprior) {
    
    // generate new phi_prop
    propose_phi();
    
    // transform phi_prop to theta_prop
    phi_prop_to_theta_prop();
    
    // calculate adjustment factor, taking into account forwards and backwards
    // moves
    get_adjustment();
    
    // calculate likelihood and prior of proposed theta
    loglike_prop = Rcpp::as<double>(get_loglike(theta_prop, s_ptr->x));
    logprior_prop = Rcpp::as<double>(get_logprior(theta_prop));
    
    // calculate Metropolis-Hastings ratio
    double MH = (loglike_prop + logprior_prop) - (loglike + logprior) + adj;
    
    // accept or reject move
    if (log(runif_0_1()) < MH) {
      
      // update theta and phi
      theta = theta_prop;
      phi = phi_prop;
      
      // update likelihoods
      loglike = loglike_prop;
      logprior = logprior_prop;
      
    }
    
  }  // end update function
  
  // other public methods
  void propose_phi();
  void phi_prop_to_theta_prop();
  void theta_to_phi();
  void get_adjustment();
  
};
