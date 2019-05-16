
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
  Rcpp::List args_progress_burnin = args_progress["pb_burnin"];
  
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
  vector<bool> bw_update = rcpp_to_vector_bool(args_params["bw_update"]);
  vector<bool> cov_update = rcpp_to_vector_bool(args_params["cov_update"]);
  vector<bool> coupling_on = rcpp_to_vector_bool(args_params["coupling_on"]);
  double GTI_pow = rcpp_to_double(args_params["GTI_pow"]);
  bool pb_markdown = rcpp_to_bool(args_params["pb_markdown"]);
  bool silent = rcpp_to_bool(args_params["silent"]);
  int chain = rcpp_to_int(args_params["chain"]);
  int input_type = rcpp_to_int(args_params["input_type"]);
  
  // extract R functions
  Rcpp::Function r_get_loglike = args_functions["r_loglike"];
  Rcpp::Function r_get_logprior = args_functions["r_logprior"];
  Rcpp::Function test_convergence = args_functions["test_convergence"];
  Rcpp::Function update_progress = args_functions["update_progress"];
  
  // extract C functions
  pattern_c_loglike* c_get_loglike = (pattern_c_loglike*) R_ExternalPtrAddrFn(args_functions["c_loglike"]);
  pattern_c_logprior* c_get_logprior = (pattern_c_logprior*) R_ExternalPtrAddrFn(args_functions["c_logprior"]);
  
  // extract C++ functions
  SEXP cpp_loglike = args_functions["cpp_loglike"];
  pattern_cpp_loglike cpp_get_loglike = *Rcpp::XPtr<pattern_cpp_loglike>(cpp_loglike);
  
  SEXP cpp_logprior = args_functions["cpp_logprior"];
  pattern_cpp_logprior cpp_get_logprior = *Rcpp::XPtr<pattern_cpp_logprior>(cpp_logprior);
  
  
  // initialise vector of particles
  vector<Particle> particle_vec(rungs);
  for (int r=0; r<rungs; ++r) {
    double beta_raised = (rungs == 1) ? 1 : pow((r + 1)/double(rungs), GTI_pow);
    particle_vec[r].init(x, theta_init, theta_min, theta_max, trans_type, beta_raised,
                         input_type, r_get_loglike, r_get_logprior, c_get_loglike, c_get_logprior,
                         cpp_get_loglike, cpp_get_logprior);
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
  
  // specify stored values at first iteration. Ensures that user-defined initial
  // values are the first stored values
  for (int r=0; r<rungs; ++r) {
    loglike_burnin[0][r][0] = particle_vec[r].loglike;
    theta_burnin[0][r][0] = particle_vec[r].theta;
  }
  
  
  // ---------- burn-in MCMC ----------
  
  // print message to console
  if (!silent) {
    print("MCMC chain", chain);
  }
  
  // loop through burn-in phases
  for (int phase=0; phase<burnin_phases; ++phase) {
    
    // print message to console
    if (!silent) {
      print("burn-in phase", phase+1);
    }
    
    // loop through burn-in iterations
    for (int rep=0; rep<burnin[phase]; ++rep) {
      
      // skip over if first iteration of first phase. Ensures that user-defined
      // initial values are the first stored values
      if (phase == 0 && rep == 0) {
        continue;
      }
      
      // loop through rungs
      for (int r=0; r<rungs; ++r) {
        
        // update particles
        particle_vec[r].update(r_get_loglike, r_get_logprior);
        
        // store results
        loglike_burnin[phase][r][rep] = particle_vec[r].loglike;
        theta_burnin[phase][r][rep] = particle_vec[r].theta;
      }
      
      // update progress bars
      if (!silent) {
        int remainder = rep % int(ceil(double(burnin[phase])/100));
        if ((remainder == 0 && !pb_markdown) || ((rep+1) == burnin[phase])) {
          update_progress(args_progress_burnin, phase+1, rep+1, burnin[phase], false);
          if ((rep+1) == burnin[phase]) {
            print("");
          }
        }
      }
      
    }  // end burn-in MCMC loop
    
  }  // end loop over burn-in phases
  
  
  // ---------- sampling MCMC ----------
  
  // print message to console
  if (!silent) {
    print("sampling phase");
  }
  
  // loop through sampling iterations
  for (int rep=0; rep<samples; ++rep) {
    
    // loop through rungs
    for (int r=0; r<rungs; ++r) {
      
      // update particles
      particle_vec[r].update(r_get_loglike, r_get_logprior);
      
      // store results
      loglike_sampling[r][rep] = particle_vec[r].loglike;
      theta_sampling[r][rep] = particle_vec[r].theta;
    }
    
    // update progress bars
    if (!silent) {
      int remainder = rep % int(ceil(double(samples)/100));
      if ((remainder==0 && !pb_markdown) || ((rep+1) == samples)) {
        update_progress(args_progress, "pb_samples", rep+1, samples, false);
        if ((rep+1) == samples) {
          print("");
        }
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
  
  //Rcpp::List ret = Rcpp::List::create(Rcpp::Named("loglike") = -9);
  
  // return as Rcpp list
  Rcpp::List ret = Rcpp::List::create(Rcpp::Named("loglike_burnin") = loglike_burnin,
                                      Rcpp::Named("theta_burnin") = theta_burnin,
                                      Rcpp::Named("loglike_sampling") = loglike_sampling,
                                      Rcpp::Named("theta_sampling") = theta_sampling);
  return ret;
}
