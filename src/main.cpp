
#include "main.h"
#include "misc_v4.h"
#include "probability.h"
#include "System.h"
#include "Particle.h"

#include <chrono>

using namespace std;

// specify exact pattern that loglike and logprior function must take in C++
typedef SEXP (*pattern_cpp_loglike)(std::vector<double>, std::vector<double>);
typedef SEXP (*pattern_cpp_logprior)(std::vector<double>);

//------------------------------------------------
// Dummy function to test Rcpp working as expected
// [[Rcpp::export]]
Rcpp::List main_cpp(Rcpp::List args) {
  
  // get flags for R vs. C++ likelihood and prior functions
  Rcpp::List args_params = args["args_params"];
  bool loglike_use_cpp = rcpp_to_bool(args_params["loglike_use_cpp"]);
  bool logprior_use_cpp = rcpp_to_bool(args_params["logprior_use_cpp"]);
  
  // extract function args
  Rcpp::List args_functions = args["args_functions"];
  
  // run MCMC with either C++ or R likelihood and prior
  Rcpp::List ret;
  if (loglike_use_cpp) {
    
    // extract likelihood function
    SEXP cpp_loglike = args_functions["loglike"];
    pattern_cpp_loglike get_loglike = *Rcpp::XPtr<pattern_cpp_loglike>(cpp_loglike);
    
    if (logprior_use_cpp) {
      
      // extract prior function
      SEXP cpp_logprior = args_functions["logprior"];
      pattern_cpp_logprior get_logprior = *Rcpp::XPtr<pattern_cpp_logprior>(cpp_logprior);
      
      // run MCMC with selected functions
      ret = run_mcmc(args, get_loglike, get_logprior);
      
    } else {
      
      // extract prior function
      Rcpp::Function get_logprior = args_functions["logprior"];
      
      // run MCMC with selected functions
      ret = run_mcmc(args, get_loglike, get_logprior);
    }
    
  } else {
    
    // extract likelihood function
    Rcpp::Function get_loglike = args_functions["loglike"];
    
    if (logprior_use_cpp) {
      
      // extract prior function
      SEXP cpp_logprior = args_functions["logprior"];
      pattern_cpp_logprior get_logprior = *Rcpp::XPtr<pattern_cpp_logprior>(cpp_logprior);
      
      // run MCMC with selected functions
      ret = run_mcmc(args, get_loglike, get_logprior);
      
    } else {
      
      // extract prior function
      Rcpp::Function get_logprior = args_functions["logprior"];
      
      // run MCMC with selected functions
      ret = run_mcmc(args, get_loglike, get_logprior);
    }
    
  }
  
  // return list
  return ret;
}

//------------------------------------------------
// run MCMC
template<class TYPE1, class TYPE2>
Rcpp::List run_mcmc(Rcpp::List args, TYPE1 get_loglike, TYPE2 get_logprior) {
  
  // start timer
  chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();
  
  // create sytem object and load args
  System s;
  s.load(args);
  
  // extract R utility functions that will be called from within MCMC
  Rcpp::List args_functions = args["args_functions"];
  Rcpp::Function test_convergence = args_functions["test_convergence"];
  Rcpp::Function update_progress = args_functions["update_progress"];
  
  // extract progress bar objects
  Rcpp::List args_progress = args["args_progress"];
  Rcpp::List args_progress_burnin = args_progress["pb_burnin"];
  
  // local copies of some parameters for convenience
  int d = s.d;
  int rungs = s.rungs;
  
  // initialise vector of particles
  vector<Particle> particle_vec(rungs);
  for (int r = 0; r < rungs; ++r) {
    
    // calculate thermodynamic power of this rung
    double beta_raised = (rungs == 1) ? 1 : pow((r + 1)/double(rungs), s.GTI_pow);
    
    // initialise particle with system objects
    particle_vec[r].init(s, beta_raised);
    
    // initialise particle initial likelihood and prior values
    particle_vec[r].init_like(get_loglike, get_logprior);
  }
  
  // objects for storing loglikelihood and theta values over iterations
  vector<vector<vector<double>>> loglike_burnin(s.burnin_phases);
  vector<vector<vector<vector<double>>>> theta_burnin(s.burnin_phases);
  for (int i = 0; i < s.burnin_phases; ++i) {
    loglike_burnin[i] = vector<vector<double>>(rungs, vector<double>(s.burnin[i]));
    theta_burnin[i] = vector<vector<vector<double>>>(rungs, vector<vector<double>>(s.burnin[i], vector<double>(d)));
  }
  vector<vector<double>> loglike_sampling(rungs, vector<double>(s.samples));
  vector<vector<vector<double>>> theta_sampling(rungs, vector<vector<double>>(s.samples, vector<double>(d)));
  
  // specify stored values at first iteration. Ensures that user-defined initial
  // values are the first stored values
  for (int r = 0; r < rungs; ++r) {
    loglike_burnin[0][r][0] = particle_vec[r].loglike;
    theta_burnin[0][r][0] = particle_vec[r].theta;
  }
  
  
  // ---------- burn-in MCMC ----------
  
  // print message to console
  if (!s.silent) {
    print("MCMC chain", s.chain);
  }
  
  // loop through burn-in phases
  for (int phase = 0; phase < s.burnin_phases; ++phase) {
    
    // print message to console
    if (!s.silent) {
      print("burn-in phase", phase+1);
    }
    
    // reset bandwidth of all rungs
    if (s.bw_reset[phase]) {
      for (int r = 0; r < rungs; ++r) {
        particle_vec[r].propSD = s.bw_init;
      }
    }
    
    // reset acceptance count of all rungs
    for (int r = 0; r < rungs; ++r) {
      particle_vec[r].accept_count = 0;
    }
    
    // loop through burn-in iterations
    for (int rep = 0; rep < s.burnin[phase]; ++rep) {
      
      // skip over if first iteration of first phase. Ensures that user-defined
      // initial values are the first stored values
      if (phase == 0 && rep == 0) {
        continue;
      }
      
      // loop through rungs
      for (int r = 0; r < rungs; ++r) {
        
        // update particles
        particle_vec[r].update(get_loglike, get_logprior,
                               rep+1, s.bw_update[phase]);
        
        // store results
        loglike_burnin[phase][r][rep] = particle_vec[r].loglike;
        theta_burnin[phase][r][rep] = particle_vec[r].theta;
      }
      
      // update progress bars
      if (!s.silent) {
        int remainder = rep % int(ceil(double(s.burnin[phase])/100));
        if ((remainder == 0 && !s.pb_markdown) || ((rep+1) == s.burnin[phase])) {
          update_progress(args_progress_burnin, phase+1, rep+1, s.burnin[phase], false);
          if ((rep+1) == s.burnin[phase]) {
            print("");
          }
        }
      }
      
    }  // end burn-in MCMC loop
    
    // print phase diagnostics
    if (!s.silent) {
      double accept_rate = round(particle_vec[0].accept_count/double(s.burnin[phase]) * 1000)/10.0;
      Rcpp::Rcout << "bandwidth: " << particle_vec[0].propSD << ", acceptance rate: " << accept_rate << "%\n";
    }
    
  }  // end loop over burn-in phases
  
  
  // ---------- sampling MCMC ----------
  
  // print message to console
  if (!s.silent) {
    print("sampling phase");
  }
  
  // reset acceptance count of all rungs
  for (int r = 0; r < rungs; ++r) {
    particle_vec[r].accept_count = 0;
  }
  
  // loop through sampling iterations
  for (int rep = 0; rep < s.samples; ++rep) {
    
    // loop through rungs
    for (int r = 0; r < rungs; ++r) {
      
      // update particles
      particle_vec[r].update(get_loglike, get_logprior,
                             rep+1, false);
      
      // store results
      loglike_sampling[r][rep] = particle_vec[r].loglike;
      theta_sampling[r][rep] = particle_vec[r].theta;
    }
    
    // update progress bars
    if (!s.silent) {
      int remainder = rep % int(ceil(double(s.samples)/100));
      if ((remainder==0 && !s.pb_markdown) || ((rep+1) == s.samples)) {
        update_progress(args_progress, "pb_samples", rep+1, s.samples, false);
        if ((rep+1) == s.samples) {
          print("");
        }
      }
    }
    
  }  // end sampling MCMC loop
  
  // print final diagnostics
  if (!s.silent) {
    double accept_rate = round(particle_vec[0].accept_count/double(s.samples) * 1000)/10.0;
    Rcpp::Rcout << "acceptance rate: " << accept_rate << "%\n";
  }
  
  
  // ---------- return ----------
  
  // end timer
  chrono::high_resolution_clock::time_point t2 = chrono::high_resolution_clock::now();
  chrono::duration<double> time_span = chrono::duration_cast< chrono::duration<double> >(t2-t1);
  if (!s.silent) {
    print("\ncompleted in", time_span.count(), "seconds\n");
  }
  
  // return as Rcpp list
  Rcpp::List ret = Rcpp::List::create(Rcpp::Named("loglike_burnin") = loglike_burnin,
                                      Rcpp::Named("theta_burnin") = theta_burnin,
                                      Rcpp::Named("loglike_sampling") = loglike_sampling,
                                      Rcpp::Named("theta_sampling") = theta_sampling);
  return ret;
}

