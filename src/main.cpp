
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
  vector<int> burnin = rcpp_to_vector_int(args_params["burnin"]);
  int samples = rcpp_to_int(args_params["samples"]);
  int rungs = rcpp_to_int(args_params["rungs"]);
  int burnin_phases = rcpp_to_int(args_params["burnin_phases"]);
  vector<bool> autoconverge_on = rcpp_to_vector_bool(args_params["autoconverge_on"]);
  vector<bool> converge_test = rcpp_to_vector_bool(args_params["converge_test"]);
  vector<bool> bw_update = rcpp_to_vector_bool(args_params["bw_update"]);
  vector<bool> cov_update = rcpp_to_vector_bool(args_params["cov_update"]);
  vector<bool> coupling_on = rcpp_to_vector_bool(args_params["coupling_on"]);
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
    double beta_raised = 1.0;  // TODO - define value per rung
    particle_vec[r].init(x, theta_init, theta_min, theta_max, trans_type, beta_raised, get_loglike, get_logprior);
  }
  
  // store loglikelihood and theta values
  vector<vector<vector<double>>> loglike_burnin(burnin_phases);
  vector<vector<vector<vector<double>>>> theta_burnin(burnin_phases);
  for (int i=0; i<burnin_phases; ++i) {
    loglike_burnin[i] = vector<vector<double>>(rungs, vector<double>(burnin[i]));
    theta_burnin[i] = vector<vector<vector<double>>>(rungs, vector<vector<double>>(burnin[i], vector<double>(d)));
  }
  vector<vector<double>> loglike_sampling(rungs, vector<double>(samples));
  vector<vector<vector<double>>> theta_sampling(rungs, vector<vector<double>>(samples, vector<double>(d)));
  
  // specify stored values at first iteration
  for (int r=0; r<rungs; ++r) {
    loglike_burnin[0][r][0] = particle_vec[r].loglike;
    theta_burnin[0][r][0] = particle_vec[r].theta;
  }
  
  
  // ---------- burn-in MCMC ----------
  
  // print message to console
  if (!silent) {
    print("running MCMC");
    print("Burn-in phase");
  }
  
  // loop through burn-in phases
  for (int phase=0; phase<burnin_phases; ++phase) {
    
    // store if convergence reached in each rung, and overall
    vector<bool> convergence_reached(rungs, false);
    bool all_convergence_reached = false;
    
    // loop through burn-in iterations
    for (int rep=0; rep<burnin[phase]; ++rep) {
      
      // skip over if first iteration of first phase, so that user-defined
      // initial values are the first stored values
      if (phase == 0 && rep == 0) {
        continue;
      }
      
      // loop through rungs
      for (int r=0; r<rungs; ++r) {
        
        // update particles
        particle_vec[r].update(get_loglike, get_logprior);
        
        // store results
        loglike_burnin[phase][r][rep] = particle_vec[r].loglike;
        theta_burnin[phase][r][rep] = particle_vec[r].theta;
      }
      
      // update progress bars
      if (!silent) {
        int remainder = rep % int(ceil(double(burnin[phase])/100));
        if ((remainder == 0 && !pb_markdown) || ((rep+1) == burnin[phase])) {
          update_progress(args_progress, "pb_burnin", rep+1, burnin[phase]);
        }
      }
      
      // check for convergence
      if (autoconverge_on[phase] && (rep % converge_test[phase]) == 0) {
        
        // check for convergence of all rungs
        all_convergence_reached = true;
        for (int r=0; r<rungs; ++r) {
          if (!convergence_reached[r]) {
            convergence_reached[r] = rcpp_to_bool(test_convergence(loglike_burnin[phase][r], rep+1));
            if (!convergence_reached[r]) {
              all_convergence_reached = false;
            }
          }
        }
        if (all_convergence_reached) {
          if (!silent) {
            update_progress(args_progress, "pb_burnin", burnin[phase], burnin[phase]);
            print("   converged within", rep, "iterations");
          }
          break;
        }
        
      }  // end check for convergence
      
    }  // end burn-in MCMC loop
    
    // warning if still not converged
    if (!all_convergence_reached && !silent) {
      print("   Warning: convergence still not reached within", burnin[phase], "iterations");
    }
    
  }  // end loop over burn-in phases
  
  /*
  // ---------- sampling MCMC ----------
  
  // print message to console
  if (!silent) {
    print("Sampling phase");
  }
  
  // loop through sampling iterations
  for (int rep=0; rep<samples; ++rep) {
    
    // loop through rungs
    for (int r=0; r<rungs; ++r) {
      
      // update particles
      particle_vec[r].update(get_loglike, get_logprior);
      
      // store results
      //loglike_store[r][rep] = particle_vec[r].loglike;
    }
    
    // store results
    //theta_store[rep] = particle_vec[0].theta;
    
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
  */
  Rcpp::List ret = Rcpp::List::create(Rcpp::Named("loglike") = -9);
    
  // return as Rcpp list
  //Rcpp::List ret = Rcpp::List::create(Rcpp::Named("loglike") = loglike_store,
  //                                    Rcpp::Named("theta") = theta_store);
  return ret;
}
