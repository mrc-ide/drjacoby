#include "cpp11.hpp"
#include "cpp11/matrix.hpp"
#include "cpp11/doubles.hpp"
#include "Rmath.h"
#include "utils.h"
#include "transform.h"
#include <RProgress.h>
#include <iostream>
#include <vector>
#include <dust/r/random.hpp>
#include <dust/random/normal.hpp>

using namespace cpp11;
namespace writable = cpp11::writable;

[[cpp11::register]]
list mcmc(
    // Options
    const int chain,
    const bool burnin,
    const int iterations,
    const bool silent,
    // Parameters
    doubles_matrix<> theta_init,
    strings theta_names,
    integers transform_type,
    doubles theta_min,
    doubles theta_max,
    integers infer_parameter,
    // Data, likelihood and prior
    list data,
    function ll_f,
    function lp_f,
    writable::list misc,
    // Tuning and acceptance
    doubles_matrix<> proposal_sd_init,
    integers_matrix<> acceptance_init,
    double target_acceptance,
    // PT
    int swap,
    doubles beta_init,
    const integers swap_acceptance_init,
    // Blocks
    list blocks_list,
    const int n_unique_blocks,
    const int iteration_counter_init,
    // RNG
    cpp11::sexp rng_ptr
) {
  
  // start timer
  std::chrono::high_resolution_clock::time_point t0 =  std::chrono::high_resolution_clock::now();
  RProgress::RProgress progress("Progress [:bar] Time remaining: :eta");
  progress.set_total(iterations);
  if(!silent){
    message("Chain " + std::to_string(chain));
  }
  
  // Initialise RNG ////////////////////////////////////////////////////////////
  auto rng = dust::random::r::rng_pointer_get<dust::random::xoshiro256plus>(rng_ptr);
  auto& state = rng->state(0);
  //////////////////////////////////////////////////////////////////////////////
  
  // Initialisise variables ////////////////////////////////////////////////////
  const int n_par = theta_names.size();
  const int n_rungs = beta_init.size();
  int iteration_counter = iteration_counter_init + 1;
  double mh;
  bool mh_accept;
  double adjustment;
  int block;
  misc.push_back({"block"_nm = 0});
  
  // Initialise matrix for theta
  // double theta[n_rungs][n_par];
  std::vector<std::vector<double>> theta;
  theta.resize(n_rungs, std::vector<double>(n_par));
  for(int i = 0; i < n_rungs; ++i){
    for(int j = 0; j < n_par; ++j){
      theta[i][j] = theta_init(i,j);
    }
  }
  // Initialise vector for proposal theta
  // Proposal theta is always the theta given to the likelihood function
  // and therefore must be a named vector, which is why it is writeable::doubles
  writable::doubles theta_prop(n_par);
  for(int p = 0; p < n_par; ++p){
    theta_prop[p] = theta_init(0,p);
  }
  theta_prop.names() = theta_names;
  
  // Initialise value for transformed theta: phi
  std::vector<std::vector<double>> phi;
  phi.resize(n_rungs, std::vector<double>(n_par));
  for(int i = 0; i < n_rungs; ++i){
    for(int j = 0; j < n_par; ++j){
      phi[i][j] = theta_to_phi(theta[i][j], transform_type[j], theta_min[j], theta_max[j]);
    }
  }
  // Initialise vector for proposal phi
  std::vector<double> phi_prop(n_par);
  
  // Initialise vector to store blocked log likelihood
  std::vector<std::vector<double>> block_ll;
  block_ll.resize(n_rungs, std::vector<double>(n_unique_blocks));
  for(int b = 0; b < n_unique_blocks; ++b) {
    for(int i = 0; i < n_rungs; ++i){
      misc["block"] = as_sexp(b + 1);
      block_ll[i][b] = cpp11::as_cpp<double>(ll_f(theta_prop, data, misc));
    }
  }
  
  // Initialise vector to store proposal blocked log likelihood
  std::vector<double> block_ll_prop(n_unique_blocks);
  // Initialise vector to store rung log likelihood (summed over blocks)
  std::vector<double> ll(n_rungs);
  for(int r = 0; r < n_rungs; ++r){
    double sum_ll = 0;
    for(int b = 0; b < n_unique_blocks; ++b){
      sum_ll += block_ll[r][b];
    }
    ll[r] = sum_ll;
  }
  
  // Initialise log prior vector
  std::vector<double> lp(n_rungs);
  for(int r = 0; r < n_rungs; ++r){
    lp[r] = cpp11::as_cpp<double>(lp_f(theta_prop, misc));
  }
  // Initialise proposal log prior
  double lp_prop;
  //////////////////////////////////////////////////////////////////////////////
  
  // Outputs ///////////////////////////////////////////////////////////////////
  // Initialise output matrix
  writable::doubles_matrix<> out(iterations, n_par + 3);
  out(0, 0) = 1;
  for(int p = 0; p < n_par; ++p){
    out(0, p + 1) =  theta[0][p];
  }
  out(0, n_par + 1) = lp[0];
  out(0, n_par + 2) = ll[0];
  
  // Initialise matrix for theta all rung output (for restarting)
  // double theta_out[n_rungs][n_par];
  writable::doubles_matrix<> theta_out(n_rungs, n_par);
  writable::doubles_matrix<> proposal_sd_out(n_rungs, n_par);
  writable::integers_matrix<> acceptance_out(n_rungs, n_par);
  //////////////////////////////////////////////////////////////////////////////
  
  // Tuning ////////////////////////////////////////////////////////////////////
  // Initialise matrix for proposal sd
  // double[n_par][n_rungs]
  std::vector<std::vector<double>> proposal_sd;
  proposal_sd.resize(n_rungs, std::vector<double>(n_par));
  for(int i = 0; i < n_rungs; ++i){
    for(int j = 0; j < n_par; ++j){
      proposal_sd[i][j] = proposal_sd_init(i,j);
    }
  }

  // Initialise acceptance count matrix
  std::vector<std::vector<int>> acceptance;
  acceptance.resize(n_rungs, std::vector<int>(n_par));
  for(int i = 0; i < n_rungs; ++i){
    for(int j = 0; j < n_par; ++j){
      acceptance[i][j] = acceptance_init(i,j);
    }
  }
  
  // Initialise swap acceptance count vector
  std::vector<int> swap_acceptance(n_rungs - 1);
  for(int i = 0; i < (n_rungs - 1); ++i){
    swap_acceptance[i] = swap_acceptance_init[i];
  }
  //////////////////////////////////////////////////////////////////////////////
  
  // PT ////////////////////////////////////////////////////////////////////////
  double rung_beta;
  int index;
  
  // Index of rungs, 0 = cold rung, n_rungs - 1 = hot rung
  std::vector<int> rung_index(n_rungs);
  for(int r = 0; r < n_rungs; ++r){
    rung_index[r] = r;
  }
  // Betas for each rung
  std::vector<double> beta(n_rungs);
  for(int r = 0; r < n_rungs; ++r){
    beta[r] = beta_init[r];
  }
  //////////////////////////////////////////////////////////////////////////////
  
  // Run ///////////////////////////////////////////////////////////////////////
  for(int i = 1; i < iterations; ++i){
    if(!silent){
      progress.tick();
    }
    for(int r = 0; r < n_rungs; ++r){
      rung_beta = beta[r];
      index = rung_index[r];
      
      // Copy rung theta and phi
      for(int p = 0; p < n_par; ++p){
        theta_prop[p] = theta[index][p];
        phi_prop[p] = phi[index][p];
      }
      
      lp_prop = lp[r];
      
      for(int p = 0; p < n_par; ++p){
        if(infer_parameter[p] == 1){
          // Propose new value of phi
          phi_prop[p] = dust::random::normal(state, phi[index][p], proposal_sd[r][p]);
          // Transform new phi_prop -> theta_prop
          theta_prop[p] = phi_to_theta(phi_prop[p], transform_type[p], theta_min[p], theta_max[p]);
          
          // Copy blocked log-likelihood
          for(int b = 0; b < n_unique_blocks; ++b){
            block_ll_prop[b] = block_ll[index][b];
          }
          // For any block this parameter is associated with, update blocked log-likelihood
          writable::integers blocks(blocks_list[p]);
          for(int b = 0; b < blocks.size(); ++b) {
            block = blocks[b];
            misc["block"] = as_sexp(block);
            block_ll_prop[block - 1] = cpp11::as_cpp<double>(ll_f(theta_prop, data, misc));
          }
          // Update proposal prior
          lp_prop = cpp11::as_cpp<double>(lp_f(theta_prop, misc));
          
          // Check for NA/NaN/+Inf in likelihood or prior
          if(value_check(lp_prop) ||  value_check(sum(block_ll_prop))){
             // Extract the theta, proposal sd and acceptance rates for output
             for(int i = 0; i < n_rungs; ++i){
               for(int j = 0; j < n_par; ++j){
                 index = rung_index[i];
                 theta_out(i,j) = theta[index][j];
                 proposal_sd_out(i,j) = proposal_sd[i][j];
                 acceptance_out(i,j) = acceptance[i][j];
               }
             }
             
            return writable::list({
                "error"_nm = "NA or Inf returned by log prior function",
                "phi_prop"_nm = phi_prop,
                "theta_prop"_nm = theta_prop,
                "rung_index"_nm = rung_index,
                "proposal_sd"_nm = proposal_sd_out,
                "theta"_nm = theta_out,
                "output"_nm = out
            });
          }
          // get parameter transformation adjustment
          adjustment = get_adjustment(theta[index][p], theta_prop[p], transform_type[p], theta_min[p], theta_max[p]);
          // calculate Metropolis-Hastings ratio
          mh = rung_beta * (sum(block_ll_prop) - ll[r]) + (lp_prop - lp[r]) + adjustment;
          // accept or reject move
          mh_accept = log(dust::random::random_real<double>(state)) < mh;
          if(mh_accept){
            // Update theta
            theta[index][p] = theta_prop[p];
            // Update phi
            phi[index][p] = phi_prop[p];
            // Update blocked log likelihood
            for(int b = 0; b < n_unique_blocks; ++b){
              block_ll[index][b] = block_ll_prop[b];
            }
            // Update log likelihood
            ll[r] = sum(block_ll_prop);
            lp[r] = lp_prop;
            // Robbins monroe step
            if(burnin){
              proposal_sd[r][p] = exp(log(proposal_sd[r][p]) + (1 - target_acceptance) / sqrt(iteration_counter));
            }
            acceptance[r][p] = acceptance[r][p] + 1;
          } else {
            // Revert theta prop
            theta_prop[p] =  theta[index][p];
            // Revert phi prop
            phi_prop[p] = phi[index][p];
            // Robbins monroe step
            if(burnin){
              proposal_sd[r][p] = exp(log(proposal_sd[r][p]) - target_acceptance / sqrt(iteration_counter));
            }
          }
        }
        // Only store values for the cold chain
        if(r == 0){
          out(i, 0) = i + 1;
          // Record parameters
          for(int p = 0; p < n_par; ++p){
            out(i, p + 1) =  theta[index][p];
          }
          out(i, n_par + 1) = lp[0];
          out(i, n_par + 2) = ll[0];
        }
      }
    }
    
    // loop over rungs, starting with the hottest chain and moving to the cold
    // chain. Each time propose a swap with the next rung up
    if(swap == 1){
      for(int r = (n_rungs - 1); r >0; --r){
        double rung_beta1 = beta[r];
        double rung_beta2 = beta[r - 1];
        double loglike1 = ll[r];
        double loglike2 = ll[r - 1];
        
        double accept = (loglike2*rung_beta1 + loglike1*rung_beta2) - (loglike1*rung_beta1 + loglike2*rung_beta2);
        // accept or reject move
        bool accept_move = log(dust::random::random_real<double>(state)) < accept;
        
        if(accept_move){
          int ri1 = rung_index[r];
          int ri2 = rung_index[r - 1];
          swap_acceptance[r - 1] += 1;
          rung_index[r] = ri2;
          rung_index[r - 1] = ri1;
          ll[r] = loglike2;
          ll[r - 1] = loglike1;
        }
      }
    }
    if(swap == 2){
      // TODO: Bob, Implement efficient swapping routine
    }
    iteration_counter ++;
  }
  
  // Extract the theta, proposal sd and acceptance rates for output
  for(int i = 0; i < n_rungs; ++i){
    for(int j = 0; j < n_par; ++j){
      index = rung_index[i];
      theta_out(i,j) = theta[index][j];
      proposal_sd_out(i,j) = proposal_sd[i][j];
      acceptance_out(i,j) = acceptance[i][j];
    }
  }
  
  // end timer
  double dur = chrono_timer(t0, "\nDuration ", !silent);
  
  // Return outputs in a list
  return writable::list({
      "output"_nm = out,
      "theta"_nm = theta_out,
      "proposal_sd"_nm = proposal_sd_out,
      "iteration_counter"_nm = iteration_counter,
      "acceptance"_nm = acceptance_out,
      "rung_index"_nm = rung_index,
      "swap_acceptance"_nm = swap_acceptance,
      "dur"_nm = dur,
      "rng_ptr"_nm = rng_ptr
  });
}
