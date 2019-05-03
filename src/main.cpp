
#include "main.h"
#include "misc_v4.h"
#include "probability.h"
#include "Particle.h"

#include <chrono>

using namespace std;

//------------------------------------------------
// Dummy function to test Rcpp working as expected
// [[Rcpp::export]]
Rcpp::List run_mcmc_cpp(Rcpp::List args) {
  
  // start timer
  chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();
  
  // split argument lists
  Rcpp::List args_params = args["args_params"];
  Rcpp::List args_functions = args["args_functions"];
  Rcpp::List args_progress = args["args_progress"];
  
  // extract inputs from Rcpp format to base C++ format
  vector<double> x = rcpp_to_vector_double(args_params["x"]);
  vector<double> theta_init = rcpp_to_vector_double(args_params["theta_init"]);
  vector<double> theta_min = rcpp_to_vector_double(args_params["theta_min"]);
  vector<double> theta_max = rcpp_to_vector_double(args_params["theta_max"]);
  vector<int> trans_type = rcpp_to_vector_int(args_params["trans_type"]);
  int d = theta_init.size();
  int burnin = rcpp_to_int(args_params["burnin"]);
  int samples = rcpp_to_int(args_params["samples"]);
  int rungs = rcpp_to_int(args_params["rungs"]);
  bool autoconverge_on = rcpp_to_bool(args_params["autoconverge_on"]);
  int converge_test = rcpp_to_int(args_params["converge_test"]);
  bool coupling_on = rcpp_to_bool(args_params["coupling_on"]);
  bool pb_markdown = rcpp_to_bool(args_params["pb_markdown"]);
  bool silent = rcpp_to_bool(args_params["silent"]);
  
  // extract R functions
  Rcpp::Function get_loglike = args_functions["loglike"];
  Rcpp::Function get_logprior = args_functions["logprior"];
  Rcpp::Function test_convergence = args_functions["test_convergence"];
  Rcpp::Function update_progress = args_functions["update_progress"];
  
  // initialise vector of particles
  vector<Particle> particle_vec(rungs);
  for (int r=0; r<rungs; ++r) {
    particle_vec[r].init(x, theta_init, theta_min, theta_max, trans_type, 1.0, get_loglike, get_logprior);
  }
  
  // objects for storing results
  vector<vector<double>> loglike_store(rungs, vector<double>(burnin+samples));
  vector<vector<double>> theta_store(burnin+samples, vector<double>(d));
  
  // specify stored values at first iteration
  for (int r=0; r<rungs; ++r) {
    loglike_store[r][0] = particle_vec[r].loglike;
  }
  theta_store[0] = particle_vec[0].theta;
  
  // store if convergence reached in each rung, and overall
  vector<bool> convergence_reached(rungs, false);
  bool all_convergence_reached = false;
  
  
  // ---------- burn-in MCMC ----------
  
  // print message to console
  if (!silent) {
    print("running MCMC");
    print("Burn-in phase");
  }
  
  // loop through burn-in iterations
  for (int rep=1; rep<burnin; ++rep) {
    
    // loop through rungs
    for (int r=0; r<rungs; ++r) {
      
      // update particles
      particle_vec[r].update(get_loglike, get_logprior);
      
      // store results
      loglike_store[r][rep] = particle_vec[r].loglike;
    }
    
    // store results
    theta_store[rep] = particle_vec[0].theta;
    
    // update progress bars
    if (!silent) {
      int remainder = rep % int(ceil(double(burnin)/100));
      if ((remainder==0 && !pb_markdown) || ((rep+1) == burnin)) {
        update_progress(args_progress, "pb_burnin", rep+1, burnin);
      }
    }
    
    // check for convergence
    if (autoconverge_on && (rep % converge_test) == 0) {
      
      // check for convergence of all rungs
      all_convergence_reached = true;
      for (int r=0; r<rungs; r++) {
        if (!convergence_reached[r]) {
          convergence_reached[r] = rcpp_to_bool(test_convergence(loglike_store[r], rep+1));
          if (!convergence_reached[r]) {
            all_convergence_reached = false;
          }
        }
      }
      if (all_convergence_reached) {
        if (!silent) {
          update_progress(args_progress, "pb_burnin", burnin, burnin);
          print("   converged within", rep, "iterations");
        }
        break;
      }
      
    }  // end check for convergence
    
  }  // end burn-in MCMC loop
  
  // warning if still not converged
  if (!all_convergence_reached && !silent) {
    print("   Warning: convergence still not reached within", burnin, "iterations");
  }
  
  
  // ---------- sampling MCMC ----------
  
  // print message to console
  if (!silent) {
    print("Sampling phase");
  }
  
  // loop through sampling iterations
  for (int rep=burnin; rep<(burnin+samples); ++rep) {
    
    // loop through rungs
    for (int r=0; r<rungs; ++r) {
      
      // update particles
      particle_vec[r].update(get_loglike, get_logprior);
      
      // store results
      loglike_store[r][rep] = particle_vec[r].loglike;
    }
    
    // store results
    theta_store[rep] = particle_vec[0].theta;
    
    // update progress bars
    if (!silent) {
      int remainder = rep % int(ceil(double(samples)/100));
      if ((remainder==0 && !pb_markdown) || ((rep+1) == samples)) {
        update_progress(args_progress, "pb_samples", rep+1, samples);
      }
    }
    
  }  // end sampling MCMC loop
  
  
  // ---------- return ----------
  
  // end timer
  chrono::high_resolution_clock::time_point t2 = chrono::high_resolution_clock::now();
  chrono::duration<double> time_span = chrono::duration_cast< chrono::duration<double> >(t2-t1);
  if (!silent) {
    print("   completed in", time_span.count(), "seconds\n");
  }
  
  // return as Rcpp list
  Rcpp::List ret = Rcpp::List::create(Rcpp::Named("loglike") = loglike_store,
                                      Rcpp::Named("theta") = theta_store);
  return ret;
}
