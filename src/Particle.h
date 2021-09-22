
#pragma once

#include "System.h"
#include "misc.h"

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
  std::vector<double> loglike_block;
  std::vector<double> loglike_prop_block;
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
    for (int i = 0; i < d; ++i) {
      for (unsigned int j = 0; j < s_ptr->block[i].size(); ++j) {
        int this_block = s_ptr->block[i][j];
        s_ptr->misc["block"] = this_block;
        loglike_block[this_block - 1] = Rcpp::as<double>(get_loglike(theta, s_ptr->x, s_ptr->misc));
      }
    }
    loglike = sum(loglike_block);
    logprior = Rcpp::as<double>(get_logprior(theta, s_ptr->misc));
    // Catch for -Inf in likelihood or prior given init theta
    if(loglike == R_NegInf || logprior == R_NegInf){
      Rcpp::Rcerr << "\n Current theta " << theta << std::endl;
      Rcpp::stop("Starting values result in -Inf in likelihood or prior. Consider setting inital values in the parameters data.frame.");
    }
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
      
      // calculate loglikelihood in each block
      loglike_prop_block = loglike_block;
      for (unsigned int j = 0; j < s_ptr->block[i].size(); ++j) {
        int this_block = s_ptr->block[i][j];
        s_ptr->misc["block"] = this_block;
        loglike_prop_block[this_block - 1] = Rcpp::as<double>(get_loglike(theta_prop, s_ptr->x, s_ptr->misc));
      }
      
      // calculate overall likelihood and prior of proposed theta
      loglike_prop = sum(loglike_prop_block);
      
      logprior_prop = Rcpp::as<double>(get_logprior(theta_prop, s_ptr->misc));
      
      // Check for NA/NaN/Inf in likelihood or prior
      if(R_IsNaN(loglike_prop) || loglike_prop == R_PosInf || R_IsNA(loglike_prop)){
        Rcpp::Rcerr << "\n Current theta " << theta_prop << std::endl;
        Rcpp::stop("NA, NaN or Inf in likelihood");
      }
      if(R_IsNaN(logprior_prop) || logprior_prop == R_PosInf || R_IsNA(logprior_prop)){
        Rcpp::Rcerr << "\n Current theta " << theta_prop << std::endl;
        Rcpp::stop("NA, NaN or Inf in prior");
      }
      
      // calculate Metropolis-Hastings ratio
      double MH;
      if(beta_raised == 0.0){
        MH = (logprior_prop - logprior) + adj;
      } else {
        MH = beta_raised*(loglike_prop - loglike) + (logprior_prop - logprior) + adj;
      }

      // accept or reject move
      bool MH_accept = (log(R::runif(0,1)) < MH);
      
      // implement changes
      if (MH_accept) {
        // update theta and phi
        theta[i] = theta_prop[i];
        phi[i] = phi_prop[i];
        
        // update likelihoods
        loglike_block = loglike_prop_block;
        loglike = loglike_prop;
        logprior = logprior_prop;
        
        // Robbins-Monro positive update  (on the log scale)
        bw[i] = exp(log(bw[i]) + bw_stepsize*(1 - 0.234) / sqrt(bw_index[i]));
        bw_index[i]++;
        
        // add to acceptance rate count
        accept_count++;
        
      } else {
        // reset theta_prop and phi_prop
        theta_prop[i] = theta[i];
        phi_prop[i] = phi[i];
        
        // Robbins-Monro negative update (on the log scale)
        bw[i] = exp(log(bw[i]) - bw_stepsize*0.234 / sqrt(bw_index[i]));
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
