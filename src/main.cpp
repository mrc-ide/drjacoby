
#include "main.h"
#include "misc.h"
#include "System.h"

#include <chrono>

using namespace std;

// specify exact pattern that loglike and logprior function must take in C++
typedef SEXP (*pattern_cpp_loglike)(Rcpp::NumericVector, Rcpp::List, Rcpp::List);
typedef SEXP (*pattern_cpp_logprior)(Rcpp::NumericVector, Rcpp::List);


//------------------------------------------------
// main MCMC function
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
  chrono::high_resolution_clock::time_point t0 = chrono::high_resolution_clock::now();
  
  // create sytem object and load args
  System s;
  s.load(args);
  
  // extract R utility functions that will be called from within MCMC
  Rcpp::List args_functions = args["args_functions"];
  Rcpp::Function test_convergence = args_functions["test_convergence"];
  Rcpp::Function update_progress = args_functions["update_progress"];
  
  // extract progress bar objects
  Rcpp::List args_progress = args["args_progress"];
  
  // local copies of some parameters for convenience
  int d = s.d;
  int rungs = s.rungs;
  
  // initialise vector of particles
  vector<Particle> particle_vec(rungs);
  for (int r = 0; r < rungs; ++r) {
    
    // initialise particle
    particle_vec[r].init(s, s.beta_raised[r]);
    
    // initialise particle initial likelihood and prior values
    particle_vec[r].init_like(get_loglike, get_logprior);
  }
  
  // objects for storing loglikelihood and theta values over iterations
  vector<vector<double>> loglike_burnin(rungs, vector<double>(s.burnin));
  vector<vector<double>> logprior_burnin(rungs, vector<double>(s.burnin));
  int theta_rungs = s.save_hot_draws ? rungs : 1;
  vector<vector<vector<double>>> theta_burnin(theta_rungs, vector<vector<double>>(s.burnin, vector<double>(d)));
  vector<vector<double>> loglike_sampling(rungs, vector<double>(s.samples));
  vector<vector<double>> logprior_sampling(rungs, vector<double>(s.samples));
  vector<vector<vector<double>>> theta_sampling(theta_rungs, vector<vector<double>>(s.samples, vector<double>(d)));
  
  // specify stored values at first iteration. Ensures that user-defined initial
  // values are the first stored values
  for (int r = 0; r < rungs; ++r) {
    loglike_burnin[r][0] = particle_vec[r].loglike;
    logprior_burnin[r][0] = particle_vec[r].logprior;
    if (s.save_hot_draws) {
      theta_burnin[r][0] = Rcpp::as< vector<double> >(particle_vec[r].theta);
    }
  }
  if (!s.save_hot_draws) {
    theta_burnin[0][0] = Rcpp::as< vector<double> >(particle_vec[rungs - 1].theta);
  }
  
  // store Metropolis coupling acceptance rates
  vector<int> mc_accept_burnin(rungs - 1);
  vector<int> mc_accept_sampling(rungs - 1);
  
  
  // ---------- burn-in MCMC ----------
  
  // print message to console
  if (!s.silent) {
    print("MCMC chain", s.chain);
    print("burn-in");
  }
  
  // loop through burn-in iterations
  for (int rep = 1; rep < s.burnin; ++rep) {
    
    // allow user to exit on escape
    Rcpp::checkUserInterrupt();
    
    // loop through rungs
    for (int r = 0; r < rungs; ++r) {
      
      // update particles
      particle_vec[r].update(get_loglike, get_logprior);
      
      // store results
      loglike_burnin[r][rep] = particle_vec[r].loglike;
      logprior_burnin[r][rep] = particle_vec[r].logprior;
      if (s.save_hot_draws) {
        theta_burnin[r][rep] = Rcpp::as< vector<double> >(particle_vec[r].theta);
      }
    }
    if (!s.save_hot_draws) {
      theta_burnin[0][rep] = Rcpp::as< vector<double> >(particle_vec[rungs - 1].theta);
    }
    
    // perform Metropolis coupling
    if (s.coupling_on) {
      coupling(particle_vec, mc_accept_burnin);
    }
    
    // update progress bars
    if (!s.silent) {
      int remainder = rep % int(ceil(double(s.burnin) / 100));
      if ((remainder == 0 && !s.pb_markdown) || ((rep+1) == s.burnin)) {
        update_progress(args_progress, "pb_burnin", rep+1, s.burnin, false);
        if ((rep+1) == s.burnin) {
          print("");
        }
      }
    }
    
  }  // end burn-in MCMC loop
  
  // print phase diagnostics
  if (!s.silent) {
    double accept_rate = particle_vec[rungs - 1].accept_count / double(s.burnin*d);
    Rcpp::Rcout << "acceptance rate: " << round(accept_rate*1000) / 10.0 << "%\n";
  }
  
  
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
    
    // allow user to exit on escape
    Rcpp::checkUserInterrupt();
    
    // loop through rungs
    for (int r = 0; r < rungs; ++r) {
      
      // update particles
      particle_vec[r].update(get_loglike, get_logprior);
      
      // store results
      loglike_sampling[r][rep] = particle_vec[r].loglike;
      logprior_sampling[r][rep] = particle_vec[r].logprior;
      if (s.save_hot_draws) {
        theta_sampling[r][rep] = Rcpp::as< vector<double> >(particle_vec[r].theta);
      }
    }
    if (!s.save_hot_draws) {
      theta_sampling[0][rep] = Rcpp::as< vector<double> >(particle_vec[rungs - 1].theta);
    }
    
    // perform Metropolis coupling
    if (s.coupling_on) {
      coupling(particle_vec, mc_accept_sampling);
    }
    
    // update progress bars
    if (!s.silent) {
      int remainder = rep % int(ceil(double(s.samples) / 100));
      if ((remainder == 0 && !s.pb_markdown) || ((rep + 1) == s.samples)) {
        update_progress(args_progress, "pb_samples", rep + 1, s.samples, false);
        if ((rep+1) == s.samples) {
          print("");
        }
      }
    }
    
  }  // end sampling MCMC loop
  
  // print final diagnostics
  if (!s.silent) {
    double accept_rate = particle_vec[rungs-1].accept_count / double(s.samples*d);
    Rcpp::Rcout << "acceptance rate: " << round(accept_rate*1000) / 10.0 << "%\n";
  }
  
  
  // ---------- return ----------
  
  
  
  // end timer
  double t_diff = chrono_timer(t0, "chain completed in ", !s.silent);
  
  // return as Rcpp list
  Rcpp::List ret = Rcpp::List::create(Rcpp::Named("loglike_burnin") = loglike_burnin,
                                      Rcpp::Named("logprior_burnin") = logprior_burnin,
                                      Rcpp::Named("theta_burnin") = theta_burnin,
                                      Rcpp::Named("loglike_sampling") = loglike_sampling,
                                      Rcpp::Named("logprior_sampling") = logprior_sampling,
                                      Rcpp::Named("theta_sampling") = theta_sampling,
                                      Rcpp::Named("mc_accept_burnin") = mc_accept_burnin,
                                      Rcpp::Named("mc_accept_sampling") = mc_accept_sampling,
                                      Rcpp::Named("t_diff") = t_diff);
  return ret;
}

//------------------------------------------------
// Metropolis-coupling over temperature rungs
void coupling(vector<Particle> &particle_vec, vector<int> &mc_accept) {
  
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
    double acceptance;
    if(beta_raised1 == 0.0){
      acceptance = (loglike1*beta_raised2) - (loglike2*beta_raised2);
    } else {
      acceptance = (loglike2*beta_raised1 + loglike1*beta_raised2) - (loglike1*beta_raised1 + loglike2*beta_raised2);
    }
    
    // accept or reject move
    bool accept_move = (log(R::runif(0,1)) < acceptance);
    
    // implement swap
    if (accept_move) {
      
      // swap parameter values
      Rcpp::NumericVector theta_tmp = Rcpp::clone(particle_vec[rung1].theta);
      particle_vec[rung1].theta = Rcpp::clone(particle_vec[rung2].theta);
      particle_vec[rung2].theta = Rcpp::clone(theta_tmp);
      
      Rcpp::NumericVector phi_tmp = Rcpp::clone(particle_vec[rung1].phi);
      particle_vec[rung1].phi = Rcpp::clone(particle_vec[rung2].phi);
      particle_vec[rung2].phi = Rcpp::clone(phi_tmp);
      
      // swap loglikelihoods
      vector<double> loglike_block_tmp = particle_vec[rung1].loglike_block;
      particle_vec[rung1].loglike_block = particle_vec[rung2].loglike_block;
      particle_vec[rung2].loglike_block = loglike_block_tmp;
      
      double loglike_tmp = particle_vec[rung1].loglike;
      particle_vec[rung1].loglike = particle_vec[rung2].loglike;
      particle_vec[rung2].loglike = loglike_tmp;
      
      double logprior_tmp = particle_vec[rung1].logprior;
      particle_vec[rung1].logprior = particle_vec[rung2].logprior;
      particle_vec[rung2].logprior = logprior_tmp;
      
      // update acceptance rates
      mc_accept[i]++;
    }
    
  }  // end loop over rungs
  
}

