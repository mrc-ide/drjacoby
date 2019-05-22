
#include "main.h"
#include "misc_v7.h"
#include "probability_v3.h"
#include "System.h"

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
    
    // initialise particle
    particle_vec[r].init(s, beta_raised);
    
    // initialise particle initial likelihood and prior values
    particle_vec[r].init_like(get_loglike, get_logprior);
  }
  
  // initialise vector for storing the order of temperature rungs. When swapping
  // particles through Metropolis coupling it is faster to update this vector to
  // represent swapping temperatures, rather than swapping parameters around
  // within particles.
  vector<int> rung_order = seq_int(0,rungs-1);
  
  // objects for storing loglikelihood and theta values over iterations
  vector<vector<vector<double>>> loglike_burnin(s.burnin_phases);
  vector<vector<vector<vector<double>>>> theta_burnin(s.burnin_phases);
  vector<vector<vector<vector<double>>>> phi_burnin(s.burnin_phases);
  for (int i = 0; i < s.burnin_phases; ++i) {
    loglike_burnin[i] = vector<vector<double>>(rungs, vector<double>(s.burnin[i]));
    theta_burnin[i] = vector<vector<vector<double>>>(rungs, vector<vector<double>>(s.burnin[i], vector<double>(d)));
    phi_burnin[i] = theta_burnin[i];
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
    
    // reset bandwidth Robbins-Monro index of all rungs
    if (s.bw_reset[phase]) {
      for (int r = 0; r < rungs; ++r) {
        particle_vec[r].bw_index = vector<int>(d, 1);
      }
    }
    
    // reset phi covariance elements of all rungs
    if (s.cov_recalc[phase]) {
      for (int r = 0; r < rungs; ++r) {
        particle_vec[r].phi_sum = vector<double>(d);
        particle_vec[r].phi_sumsq = vector<vector<double>>(d, vector<double>(d));
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
        if (s.prop_method[phase] == 1) {
          particle_vec[r].update_univar(get_loglike, get_logprior, s.bw_update[phase]);
        } else {
          particle_vec[r].update_multivar(get_loglike, get_logprior, s.bw_update[phase]);
        }
        
        // update phi sum and sum-of-squares
        if (s.cov_recalc[phase]) {
          particle_vec[r].update_phi_sumsq();
        }
        
        // store results
        loglike_burnin[phase][r][rep] = particle_vec[r].loglike;
        theta_burnin[phase][r][rep] = particle_vec[r].theta;
        phi_burnin[phase][r][rep] = particle_vec[r].phi;
      }
      
      // perform Metropolis coupling
      if (s.coupling_on) {
        coupling(particle_vec);
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
    
    // update phi covariance
    if (s.cov_recalc[phase]) {
      for (int r = 0; r < rungs; ++r) {
        particle_vec[r].get_phi_cov(s.burnin[phase]);
      }
    }
    
    // print phase diagnostics
    if (!s.silent) {
      double accept_rate = particle_vec[rungs-1].accept_count/double(s.burnin[phase]);
      if (s.prop_method[phase] == 1) {
        accept_rate /= d;
      }
      Rcpp::Rcout << "acceptance rate: " << round(accept_rate*1000)/10.0 << "%\n";
    }
    
  }  // end loop over burn-in phases
  
  
  // ---------- sampling MCMC ----------
  
  // print message to console
  if (!s.silent) {
    print("sampling phase");
  }
  
  // get final sampling method
  int prop_method_sampling = s.prop_method[s.burnin_phases-1];
  
  // reset acceptance count of all rungs
  for (int r = 0; r < rungs; ++r) {
    particle_vec[r].accept_count = 0;
  }
  
  // loop through sampling iterations
  for (int rep = 0; rep < s.samples; ++rep) {
    
    // loop through rungs
    for (int r = 0; r < rungs; ++r) {
      
      // update particles
      if (prop_method_sampling == 1) {
        particle_vec[r].update_univar(get_loglike, get_logprior, false);
      } else {
        particle_vec[r].update_multivar(get_loglike, get_logprior, false);
      }
      
      // store results
      loglike_sampling[r][rep] = particle_vec[r].loglike;
      theta_sampling[r][rep] = particle_vec[r].theta;
    }
    
    // perform Metropolis coupling
    if (s.coupling_on) {
      coupling(particle_vec);
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
    double accept_rate = particle_vec[rungs-1].accept_count/double(s.samples);
    if (prop_method_sampling == 1) {
      accept_rate /= d;
    }
    Rcpp::Rcout << "acceptance rate: " << round(accept_rate*1000)/10.0 << "%\n";
  }
  
  
  // ---------- return ----------
  
  // end timer
  if (!s.silent) {
    print("");
    chrono_timer(t1);
  }
  
  // return as Rcpp list
  Rcpp::List ret = Rcpp::List::create(Rcpp::Named("loglike_burnin") = loglike_burnin,
                                      Rcpp::Named("theta_burnin") = theta_burnin,
                                      Rcpp::Named("phi_burnin") = phi_burnin,
                                      Rcpp::Named("loglike_sampling") = loglike_sampling,
                                      Rcpp::Named("theta_sampling") = theta_sampling);
  return ret;
}

//------------------------------------------------
// Metropolis-coupling over temperature rungs
void coupling(vector<Particle> &particle_vec) {
  
  // get number of rungs
  int rungs = int(particle_vec.size());
  
  // loop over rungs, starting with the hottest chain and moving to the cold
  // chain. Each time propose a swap with the next rung up
  for (int i = 0; i < (rungs-1); ++i) {
    
    // define rungs of interest
    int rung1 = i;
    int rung2 = i+1;
    
    // get log-likelihoods and beta values of two chains in the comparison
    double loglike1 = particle_vec[rung1].loglike;
    double loglike2 = particle_vec[rung2].loglike;
    
    double beta_raised1 = particle_vec[rung1].beta_raised;
    double beta_raised2 = particle_vec[rung2].beta_raised;
    
    // calculate acceptance ratio (still in log space)
    double acceptance = (loglike2*beta_raised1 + loglike1*beta_raised2) - (loglike1*beta_raised1 + loglike2*beta_raised2);
    
    // accept or reject move
    bool accept_move = (log(runif_0_1()) < acceptance);
    
    // implement swap
    if (accept_move) {
      
      // swap parameter values
      vector<double> theta_tmp = particle_vec[rung1].theta;
      particle_vec[rung1].theta = particle_vec[rung2].theta;
      particle_vec[rung2].theta = theta_tmp;
      
      vector<double> phi_tmp = particle_vec[rung1].phi;
      particle_vec[rung1].phi = particle_vec[rung2].phi;
      particle_vec[rung2].phi = phi_tmp;
      
      // swap loglikelihoods
      double loglike_tmp = particle_vec[rung1].loglike;
      particle_vec[rung1].loglike = particle_vec[rung2].loglike;
      particle_vec[rung2].loglike = loglike_tmp;
      
      double logprior_tmp = particle_vec[rung1].logprior;
      particle_vec[rung1].logprior = particle_vec[rung2].logprior;
      particle_vec[rung2].logprior = logprior_tmp;
      
      // update acceptance rates
      //coupling_accept[i]++;
    }
    
  }  // end loop over rungs
  
}
