
#include "System.h"
#include "Particle.h"

#include <Rcpp.h>

//------------------------------------------------
// main Rcpp function, deployed from R
Rcpp::List main_cpp(Rcpp::List args);

//------------------------------------------------
// run MCMC
template<class TYPE1, class TYPE2>
Rcpp::List run_mcmc(Rcpp::List args, TYPE1 get_loglike, TYPE2 get_logprior);

//------------------------------------------------
// Metropolis-coupling over temperature rungs
void coupling(std::vector<Particle> &particle_vec, std::vector<int> &mc_accept);
