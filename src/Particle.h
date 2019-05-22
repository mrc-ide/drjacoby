
#pragma once

#include "System.h"
#include "misc_v7.h"
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
  
  // beta_raised stores values of beta (the thermodynamic power), raised to the
  // power GTI_pow
  double beta_raised;
  
  // theta is the parameter vector in natural space
  std::vector<double> theta;
  std::vector<double> theta_prop;
  
  // phi is a vector of transformed parameters
  std::vector<double> phi;
  std::vector<double> phi_prop;
  
  // proposal parameters
  std::vector<double> bw;
  double bw_multi;
  std::vector<int> bw_index;
  double bw_stepsize;
  std::vector<double> phi_sum;
  std::vector<std::vector<double>> phi_sumsq;
  std::vector<std::vector<double>> phi_cov;
  std::vector<std::vector<double>> phi_cholesky;
  
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
  
  // initialise everything EXCEPT FOR likelihood and prior values
  void init(System &s, double beta_raised);
  
  // initialise likelihood and prior values
  template<class TYPE1, class TYPE2>
  void init_like(TYPE1 get_loglike, TYPE2 get_logprior) {
    loglike = Rcpp::as<double>(get_loglike(theta, s_ptr->x));
    logprior = Rcpp::as<double>(get_logprior(theta));
  }
  
  // update theta[i] via univariate Metropolis-Hastings
  template<class TYPE1, class TYPE2>
  void update_univar(TYPE1 get_loglike, TYPE2 get_logprior, bool update_bw) {
    
    // set theta_prop and phi_prop
    theta_prop = theta;
    phi_prop = phi;
    
    // loop through parameters
    for (int i = 0; i < d; ++i) {
      
      // generate new phi_prop[i]
      propose_phi_univar(i);
      
      // transform phi_prop[i] to theta_prop[i]
      phi_prop_to_theta_prop(i);
      
      // calculate adjustment factor, taking into account forwards and backwards
      // moves
      double adj = get_adjustment(i);
      
      // calculate likelihood and prior of proposed theta
      loglike_prop = Rcpp::as<double>(get_loglike(theta_prop, s_ptr->x));
      logprior_prop = Rcpp::as<double>(get_logprior(theta_prop));
      
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
        if (update_bw) {
          bw[i] = exp(log(bw[i]) + bw_stepsize*(1 - 0.234)/sqrt(bw_index[i]));
          bw_index[i]++;
        }
        
        // add to acceptance rate count
        accept_count++;
        
      } else {
        
        // reset theta_prop and phi_prop
        theta_prop[i] = theta[i];
        phi_prop[i] = phi[i];
        
        // Robbins-Monro negative update (on the log scale)
        if (update_bw) {
          bw[i] = exp(log(bw[i]) - bw_stepsize*0.234/sqrt(bw_index[i]));
          bw_index[i]++;
        }
        
      } // end MH step
      
    }  // end loop over parameters
    
  }  // end update_univar function
  
  // update all theta via multivariate Metropolis-Hastings
  template<class TYPE1, class TYPE2>
  void update_multivar(TYPE1 get_loglike, TYPE2 get_logprior, bool update_bw) {
    
    // generate new phi_prop
    propose_phi_multivar();
    
    // transform phi_prop to theta_prop
    for (int i = 0; i < d; ++i) {
      phi_prop_to_theta_prop(i);
    }
    
    // calculate adjustment factor, taking into account forwards and backwards
    // moves
    double adj = 0;
    for (int i = 0; i < d; ++i) {
      adj += get_adjustment(i);
    }
    
    // calculate likelihood and prior of proposed theta
    loglike_prop = Rcpp::as<double>(get_loglike(theta_prop, s_ptr->x));
    logprior_prop = Rcpp::as<double>(get_logprior(theta_prop));
    
    // calculate Metropolis-Hastings ratio
    double MH = beta_raised*(loglike_prop - loglike) + (logprior_prop - logprior) + adj;
    
    // accept or reject move
    bool MH_accept = (log(runif_0_1()) < MH);
    
    // implement changes
    if (MH_accept) {
      
      // update theta and phi
      theta = theta_prop;
      phi = phi_prop;
      
      // update likelihoods
      loglike = loglike_prop;
      logprior = logprior_prop;
      
      // Robbins-Monro positive update  (on the log scale)
      if (update_bw) {
        bw_multi = exp(log(bw_multi) + bw_stepsize*(1 - 0.234)/sqrt(bw_index[0]));
        bw_index[0]++;
      }
      
      // add to acceptance rate count
      accept_count++;
      
    } else {
      
      // Robbins-Monro negative update (on the log scale)
      if (update_bw) {
        bw_multi = exp(log(bw_multi) - bw_stepsize*0.234/sqrt(bw_index[0]));
        bw_index[0]++;
      }
      
    }
    
  }  // end update_multivar function
  
  // other public methods
  void propose_phi_univar(int i);
  void propose_phi_multivar();
  void phi_prop_to_theta_prop(int i);
  void theta_to_phi();
  double get_adjustment(int i);
  void update_phi_sumsq();
  void get_phi_cov(int n);
  
};
