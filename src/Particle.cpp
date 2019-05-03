
#include "Particle.h"
#include "misc_v4.h"
#include "probability.h"

using namespace std;

//------------------------------------------------
// initialise/reset particle
void Particle::init(vector<double> &x,
                    vector<double> &theta_init,
                    vector<double> &theta_min,
                    vector<double> &theta_max,
                    vector<int> &trans_type,
                    double beta_raised,
                    Rcpp::Function get_loglike,
                    Rcpp::Function get_logprior) {
  
  // pointer to data
  x_ptr = &x;
  
  // theta is the parameter vector in natural space
  theta = theta_init;
  this->theta_min = theta_min;
  this->theta_max = theta_max;
  d = theta.size();
  theta_prop = vector<double>(d);
  
  // the type of transformation applied to each element of theta. See main.R for
  // a key
  this->trans_type = trans_type;
  
  // adjustment factor to account for reparameterisation
  adj = 0;
  
  // phi is a vector of transformed parameters
  phi = vector<double>(d);
  theta_to_phi();
  phi_prop = vector<double>(d);
  
  // beta_raised stores values of beta (the thermodynamic power), raised to the
  // power GTI_pow
  this->beta_raised = beta_raised;
  
  // proposal standard deviations
  propSD = 1;
  
  // likelihoods and priors
  loglike = rcpp_to_double(get_loglike(theta, *x_ptr));
  loglike_prop = 0;
  logprior = rcpp_to_double(get_logprior(theta));;
  logprior_prop = 0;
  
  // acceptance rates
  accept = 0;
  
}

//------------------------------------------------
// update particle
void Particle::update(Rcpp::Function get_loglike, Rcpp::Function get_logprior) {
  
  // propose new phi
  propose_phi();
  
  // transform phi_prop to theta_prop
  phi_prop_to_theta_prop();
  
  // calculate adjustment factor forwards and backwards
  get_adjustment();
  
  // calculate loglikelihood of proposed theta
  loglike_prop = rcpp_to_double(get_loglike(theta_prop, *x_ptr));
  
  // calculate logprior of proposed theta
  logprior_prop = rcpp_to_double(get_logprior(theta_prop));
  
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
  
}

//------------------------------------------------
// propose new value of phi from multivariate normal distribution centred on
// current phi
void Particle::propose_phi() {
  
  // TODO - this function currently ignores correlations
  
  for (int i=0; i<d; ++i) {
    phi_prop[i] = rnorm1(phi[i], propSD);
  }
  
}

//------------------------------------------------
// transform phi_prop to theta_prop. See main.R for a key to transformation
// types
void Particle::phi_prop_to_theta_prop() {
  
  for (int i=0; i<d; ++i) {
    switch(trans_type[i]) {
    case 0:
      theta_prop[i] = phi_prop[i];
      break;
    case 1:
      theta_prop[i] = theta_max[i] - exp(phi_prop[i]);
      break;
    case 2:
      theta_prop[i] = exp(phi_prop[i]) + theta_min[i];
      break;
    case 3:
      theta_prop[i] = (theta_max[i]*exp(phi_prop[i]) + theta_min[i]) / (1 + exp(phi_prop[i]));
      break;
    default:
      Rcpp::stop("trans_type invalid");
    }
  }
  
}

//------------------------------------------------
// transform theta to phi. See main.R for a key to transformation types
void Particle::theta_to_phi() {
  
  for (int i=0; i<d; ++i) {
    switch(trans_type[i]) {
    case 0:
      phi[i] = theta[i];
      break;
    case 1:
      phi[i] = log(theta_max[i] - theta[i]);
      break;
    case 2:
      phi[i] = log(theta[i] - theta_min[i]);
      break;
    case 3:
      phi[i] = log(theta[i] - theta_min[i]) - log(theta_max[i] - theta[i]);
      break;
    default:
      Rcpp::stop("trans_type invalid");
    }
  }
  
}

//------------------------------------------------
// get adjustment factor to account for reparameterisation
void Particle::get_adjustment() {
  
  double ret = 0;
  for (int i=0; i<d; ++i) {
    switch(trans_type[i]) {
    case 0:
      // (no adjustment needed)
      break;
    case 1:
      ret += log(theta_prop[i] - theta_max[i]) - log(theta[i] - theta_max[i]);
      break;
    case 2:
      ret += log(theta_prop[i] - theta_min[i]) - log(theta[i] - theta_min[i]);
      break;
    case 3:
      ret += log(theta_max[i] - theta_prop[i]) + log(theta_prop[i] - theta_min[i]) - log(theta_max[i] - theta[i]) - log(theta[i] - theta_min[i]);
      break;
    default:
      Rcpp::stop("trans_type invalid");
    }
  }
  
}
