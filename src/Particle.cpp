
#include "Particle.h"

using namespace std;

//------------------------------------------------
// initialise/reset particle
void Particle::init(System &s, double beta_raised) {
  
  // pointer to system object
  this->s_ptr = &s;
  
  // local copies of some parameters for convenience
  d = s_ptr->d;
  
  // beta_raised stores values of beta (the thermodynamic power), raised to the
  // power GTI_pow
  this->beta_raised = beta_raised;
  
  // theta is the parameter vector in natural space
  theta = s_ptr->theta_init;
  theta_prop = vector<double>(d);
  
  // phi is a vector of transformed parameters
  phi = vector<double>(d);
  theta_to_phi();
  phi_prop = vector<double>(d);
  
  // proposal parameters
  bw = vector<double>(d, 1.0);
  bw_multi = 1.0;
  bw_index = vector<int>(d, 1);
  bw_stepsize = 1.0;
  phi_sum = vector<double>(d);
  phi_sumsq = vector<vector<double>>(d, vector<double>(d));
  phi_cov = vector<vector<double>>(d, vector<double>(d));
  phi_cholesky = vector<vector<double>>(d, vector<double>(d));
  for (int i = 0; i < d; ++i) {
    phi_cholesky[i][i] = 1.0;
  }
  
  // likelihoods and priors
  loglike = 0;
  loglike_prop = 0;
  logprior = 0;
  logprior_prop = 0;
  
  // acceptance rates
  accept_count = 0;
  
}

//------------------------------------------------
// propose new value of phi[i] from univariate normal distribution
void Particle::propose_phi_univar(int i) {
  phi_prop[i] = rnorm1(phi[i], bw[i]);
}

//------------------------------------------------
// propose new value of phi from multivariate normal distribution centred on
// current phi
void Particle::propose_phi_multivar() {
  rmnorm1(phi_prop, phi, phi_cholesky, bw_multi);
}

//------------------------------------------------
// transform phi_prop to theta_prop. See main.R for a key to transformation
// types
void Particle::phi_prop_to_theta_prop(int i) {
  
  switch(s_ptr->trans_type[i]) {
  case 0:
    theta_prop[i] = phi_prop[i];
    break;
  case 1:
    theta_prop[i] = s_ptr->theta_max[i] - exp(phi_prop[i]);
    break;
  case 2:
    theta_prop[i] = exp(phi_prop[i]) + s_ptr->theta_min[i];
    break;
  case 3:
    theta_prop[i] = (s_ptr->theta_max[i]*exp(phi_prop[i]) + s_ptr->theta_min[i]) / (1 + exp(phi_prop[i]));
    break;
  default:
    Rcpp::stop("trans_type invalid");
  }
  
}

//------------------------------------------------
// transform theta to phi. See main.R for a key to transformation types
void Particle::theta_to_phi() {
  
  for (int i = 0; i < d; ++i) {
    switch(s_ptr->trans_type[i]) {
    case 0:
      phi[i] = theta[i];
      break;
    case 1:
      phi[i] = log(s_ptr->theta_max[i] - theta[i]);
      break;
    case 2:
      phi[i] = log(theta[i] - s_ptr->theta_min[i]);
      break;
    case 3:
      phi[i] = log(theta[i] - s_ptr->theta_min[i]) - log(s_ptr->theta_max[i] - theta[i]);
      break;
    default:
      Rcpp::stop("trans_type invalid");
    }
  }
  
}

//------------------------------------------------
// get adjustment factor to account for reparameterisation
double Particle::get_adjustment(int i) {
  
  double ret = 0;
  switch(s_ptr->trans_type[i]) {
  case 0:
    // (no adjustment needed)
    break;
  case 1:
    ret = log(theta_prop[i] - s_ptr->theta_max[i]) - log(theta[i] - s_ptr->theta_max[i]);
    break;
  case 2:
    ret = log(theta_prop[i] - s_ptr->theta_min[i]) - log(theta[i] - s_ptr->theta_min[i]);
    break;
  case 3:
    ret = log(s_ptr->theta_max[i] - theta_prop[i]) + log(theta_prop[i] - s_ptr->theta_min[i]) - log(s_ptr->theta_max[i] - theta[i]) - log(theta[i] - s_ptr->theta_min[i]);
    break;
  default:
    Rcpp::stop("trans_type invalid");
  }
  return ret;
  
}

//------------------------------------------------
// add current phi values to running sum and sum of squares
void Particle::update_phi_sumsq() {
  
  // update phi_sum and phi_sumsq
  for (int i = 0; i < d; ++i) {
    phi_sum[i] += phi[i];
    for (int j = i; j < d; ++j) {
      phi_sumsq[i][j] += phi[i]*phi[j];
    }
  }
  
}

//------------------------------------------------
// calculate phi covariance matrix and Cholesky decomposition of covariance
// matrix from sum and sum-of-squares. n is the number of values that went into
// the sum
void Particle::get_phi_cov(int n) {
  
  // calculate covariance matrix
  for (int i = 0; i < d; ++i) {
    for (int j = i; j < d; ++j) {
      phi_cov[i][j] = phi_sumsq[i][j]/double(n) - (phi_sum[i]/double(n))*(phi_sum[j]/double(n));
    }
  }
  
  // calculate Cholesky decomposition
  cholesky(phi_cholesky, phi_cov);
  
}
