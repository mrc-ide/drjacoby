
#pragma once

#include <Rcpp.h>

//------------------------------------------------
// class defining MCMC particle
class Particle {
  
public:
  // PUBLIC OBJECTS
  
  // parameter values
  std::vector<double> theta;
  std::vector<double> phi;
  
  // beta_raised stores values of beta (the thermodynamic power), raised to the 
  // power GTI_pow
  double beta_raised;
  
  // proposal standard deviations
  double propSD;
  
  // likelihoods
  double loglike_old;
  double loglike_new;
  
  // store acceptance rates
  int accept;
  
  
  // PUBLIC FUNCTIONS
  
  // constructors
  Particle() {};
  Particle(double beta_raised);
  
  // public methods
  void reset(double beta_raised);
  void update();
  void calculate_loglike();
  
};
