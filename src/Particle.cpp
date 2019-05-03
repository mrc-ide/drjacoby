
#include "Particle.h"
//#include "misc_v4.h"
//#include "probability.h"

using namespace std;

//------------------------------------------------
// constructor for Particle class
Particle::Particle(double beta_raised) {
  
  // beta_raised stores values of beta (the thermodynamic power), raised to the
  // power GTI_pow
  this->beta_raised = beta_raised;
  
  // initialise proposal standard deviations
  propSD = 1;
  
  // likelihood
  loglike_old = 0;
  loglike_new = 0;
  
  // store acceptance rates
  accept = 0;
  
}

//------------------------------------------------
// reset particle
void Particle::reset(double beta_raised) {
  
  // beta_raised
  this->beta_raised = beta_raised;
  
}

//------------------------------------------------
// update particle
void Particle::update() {
  
}

