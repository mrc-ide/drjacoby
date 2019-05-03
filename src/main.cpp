
#include "main.h"
#include "misc_v4.h"
#include "probability.h"

#include <unistd.h>  // TODO - delete when not sleeping

using namespace std;

//------------------------------------------------
// Dummy function to test Rcpp working as expected
// [[Rcpp::export]]
Rcpp::List run_mcmc_cpp(Rcpp::List args) {
  
  // split argument lists
  Rcpp::List args_params = args["args_params"];
  Rcpp::List args_functions = args["args_functions"];
  Rcpp::List args_progress = args["args_progress"];
  
  // extract inputs from Rcpp format to base C++ format
  vector<double> x = rcpp_to_vector_double(args_params["x"]);
  vector<double> theta_init = rcpp_to_vector_double(args_params["init"]);
  vector<int> trans_type = rcpp_to_vector_int(args_params["trans_type"]);
  int burnin = rcpp_to_int(args_params["burnin"]);
  int samples = rcpp_to_int(args_params["samples"]);
  int rungs = rcpp_to_int(args_params["rungs"]);
  bool autoconverge_on = rcpp_to_bool(args_params["autoconverge_on"]);
  bool coupling_on = rcpp_to_bool(args_params["coupling_on"]);
  bool pb_markdown = rcpp_to_bool(args_params["pb_markdown"]);
  bool silent = rcpp_to_bool(args_params["silent"]);
  
  // extract R functions
  Rcpp::Function loglike = args_functions["loglike"];
  Rcpp::Function test_convergence = args_functions["test_convergence"];
  Rcpp::Function update_progress = args_functions["update_progress"];
  
  //double ll = rcpp_to_double(loglike(theta_init, x));
  
  // initialise vector of particles
  
  
  
  // ---------- burn-in MCMC ----------
  
  // print message to console
  if (!silent) {
    print("running MCMC");
    print("Burn-in phase");
  }
  
  // loop through burn-in iterations
  for (int rep=0; rep<burnin; ++rep) {
    
    usleep(1000);
    
    // update progress bars
    if (!silent) {
      int remainder = rep % int(ceil(double(burnin)/100));
      if ((remainder==0 && !pb_markdown) || ((rep+1) == burnin)) {
        update_progress(args_progress, "pb_burnin", rep+1, burnin);
      }
    }
    
  }  // end burn-in MCMC loop
  
  
  // ---------- sampling MCMC ----------
  
  // print message to console
  if (!silent) {
    print("Sampling phase");
  }
  
  // loop through sampling iterations
  for (int rep=0; rep<samples; ++rep) {
    
    usleep(1000);
    
    // update progress bars
    if (!silent) {
      int remainder = rep % int(ceil(double(samples)/100));
      if ((remainder==0 && !pb_markdown) || ((rep+1) == samples)) {
        update_progress(args_progress, "pb_samples", rep+1, samples);
      }
    }
    
  }  // end sampling MCMC loop
  
  
  // ---------- return ----------
  
  // return as Rcpp list
  Rcpp::List ret = Rcpp::List::create(Rcpp::Named("bar") = trans_type);
  return ret;
}
