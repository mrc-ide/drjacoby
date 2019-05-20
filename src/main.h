
#include "System.h"

#include <Rcpp.h>

//------------------------------------------------
// main Rcpp function, deployed from R
Rcpp::List main_cpp(Rcpp::List args);

//------------------------------------------------
// run MCMC
template<class TYPE1, class TYPE2>
Rcpp::List run_mcmc(Rcpp::List args, TYPE1 get_loglike, TYPE2 get_logprior);
