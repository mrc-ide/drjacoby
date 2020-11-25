
#pragma once

#include <Rcpp.h>

#include <vector>

//------------------------------------------------
// class holding all data, parameters and functions
class System {
  
public:
  // PUBLIC OBJECTS
  
  // data
  Rcpp::List x;
  
  // misc object
  Rcpp::List misc;
  
  // model parameters
  Rcpp::NumericVector theta_vector;
  std::vector<double> theta_min;
  std::vector<double> theta_max;
  std::vector<int> trans_type;
  std::vector<bool> skip_param;
  int d;
  
  // MCMC parameters
  int burnin;
  int samples;
  int rungs;
  bool coupling_on;
  std::vector<double> beta_raised;
  int chain;
  
  // misc parameters
  bool pb_markdown;
  bool silent;
  
  // PUBLIC FUNCTIONS
  
  // constructors
  System() {};
  
  // public methods
  void load(Rcpp::List args);
  
};
