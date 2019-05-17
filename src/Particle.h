
#pragma once

#include <Rcpp.h>

// specify exact pattern that loglike and logprior function must take in C++
typedef double (*pattern_cpp_loglike)(std::vector<double>, std::vector<double>);
typedef double (*pattern_cpp_logprior)(std::vector<double>);

//------------------------------------------------
// class defining MCMC particle
class Particle {
  
public:
  // PUBLIC OBJECTS
  
  // pointer to data
  std::vector<double>* x_ptr;
  
  // theta is the parameter vector in natural space
  std::vector<double> theta;
  std::vector<double> theta_min;
  std::vector<double> theta_max;
  std::vector<double> theta_prop;
  int d;  // number of parameters
  
  // TODO - delete once chosen best method
  int input_type;
  pattern_cpp_loglike cpp_get_loglike;
  pattern_cpp_logprior cpp_get_logprior;
  
  // the type of transformation applied to each element of theta. See main.R for
  // a key
  std::vector<int> trans_type;
  
  // adjustment factor to account for reparameterisation
  double adj;
  
  // phi is a vector of transformed parameters
  std::vector<double> phi;
  std::vector<double> phi_prop;
  
  // beta_raised stores values of beta (the thermodynamic power), raised to the
  // power GTI_pow
  double beta_raised;
  
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
  
  // public methods
  void init(std::vector<double> &x,
            std::vector<double> &theta_init,
            std::vector<double> &theta_min,
            std::vector<double> &theta_max,
            std::vector<int> &trans_type,
            double beta_raised,
            int input_type,
            Rcpp::Function get_loglike,
            Rcpp::Function get_logprior,
            pattern_cpp_loglike cpp_get_loglike,
            pattern_cpp_logprior cpp_get_logprior);
  
  void update(Rcpp::Function get_loglike,
              Rcpp::Function get_logprior);
  
  void propose_phi();
  void phi_prop_to_theta_prop();
  void theta_to_phi();
  void get_adjustment();
  
  template<class TYPE>
  double get_loglike0(TYPE func, std::vector<double> &theta, std::vector<double> &x) {
    return func(theta, x);
  }
  
};
