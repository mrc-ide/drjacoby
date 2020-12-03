
#pragma once

#include "System.h"
#include "misc_v10.h"
#include "probability_v3.h"

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
  
  // thermodynamic power
  double beta_raised;
  
  // theta is the parameter vector in natural space
  Rcpp::NumericVector theta;
  Rcpp::NumericVector theta_prop;
  
  // phi is a vector of transformed parameters
  Rcpp::NumericVector phi;
  Rcpp::NumericVector phi_prop;
  
  // proposal parameters
  std::vector<double> bw;
  std::vector<int> bw_index;
  double bw_stepsize;
  
  // likelihoods and priors
  double loglike;
  double loglike_prop;
  double logprior;
  double logprior_prop;
  
  // store acceptance rates
  int accept_count;
  
  
  // PUBLIC FUNCTIONS
  
  // constructors
  Particle() {};
  
  // initialise everything except for likelihood and prior values
  void init(System &s, double beta_raised);
  
  // initialise likelihood and prior values
  template<class TYPE1, class TYPE2>
  void init_like(TYPE1 get_loglike, TYPE2 get_logprior) {
    loglike = Rcpp::as<double>(get_loglike(theta, -1, s_ptr->x, s_ptr->misc));
    logprior = Rcpp::as<double>(get_logprior(theta, -1, s_ptr->misc));
  }
  
  // update theta[i] via univariate Metropolis-Hastings
  template<class TYPE1, class TYPE2>
  void update(TYPE1 get_loglike, TYPE2 get_logprior) {
    
    // set theta_prop and phi_prop to current values of theta and phi
    theta_prop = Rcpp::clone(theta);
    phi_prop = Rcpp::clone(phi);
    
    // loop through parameters
    for (int i = 0; i < d; ++i) {
      if (s_ptr->skip_param[i]) {
        continue;
      }
      
      // generate new phi_prop[i]
      propose_phi(i);
      
      // transform phi_prop[i] to theta_prop[i]
      phi_prop_to_theta_prop(i);
      
      // calculate adjustment factor, taking into account forwards and backwards
      // moves
      double adj = get_adjustment(i);
      
      // calculate likelihood and prior of proposed theta
      loglike_prop = Rcpp::as<double>(get_loglike(theta_prop, i, s_ptr->x, s_ptr->misc));
      logprior_prop = Rcpp::as<double>(get_logprior(theta_prop, i, s_ptr->misc));
      
      // calculate Metropolis-Hastings ratio
      double MH = beta_raised*(loglike_prop - loglike) + (logprior_prop - logprior) + adj;
      
      // accept or reject move
      bool MH_accept = (log(runif_0_1()) < MH);
      
      // implement changes
      if (MH_accept) {
        
        // update theta and phi
        theta[i] = theta_prop[i];
        phi[i] = phi_prop[i];
        
        // update likelihoods
        loglike = loglike_prop;
        logprior = logprior_prop;
        
        // Robbins-Monro positive update  (on the log scale)
        bw[i] = exp(log(bw[i]) + bw_stepsize*(1 - 0.234)/sqrt(bw_index[i]));
        bw_index[i]++;
        
        // add to acceptance rate count
        accept_count++;
        
      } else {
        
        // reset theta_prop and phi_prop
        theta_prop[i] = theta[i];
        phi_prop[i] = phi[i];
        
        // Robbins-Monro negative update (on the log scale)
        bw[i] = exp(log(bw[i]) - bw_stepsize*0.234/sqrt(bw_index[i]));
        bw_index[i]++;
        
      } // end MH step
      
    }  // end loop over parameters
    
  }  // end update_univar function
  
  // other public methods
  void propose_phi(int i);
  void phi_prop_to_theta_prop(int i);
  void theta_to_phi();
  double get_adjustment(int i);
  
};
