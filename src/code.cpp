#include "cpp11.hpp"
#include "cpp11/matrix.hpp"
#include "cpp11/doubles.hpp"
#include "Rmath.h"
#include "utils.h"
#include "transform.h"
#include <vector>
using namespace cpp11;
namespace writable = cpp11::writable;


[[cpp11::register]]
list mcmc(
    doubles theta_init,
    strings theta_names,
    integers transform_type,
    doubles theta_min,
    doubles theta_max,
    list blocks_list,
    int n_unique_blocks,
    list data,
    int burnin,
    int samples,
    function ll_f,
    function lp_f,
    double target_acceptance,
    writable::list misc,
    int n_rungs,
    doubles beta_init,
    bool swap) {


  int iterations = burnin + samples;
  int n_par = theta_init.size();

  // Initialisise variables, ////////////////////////////////////////////////////
  double mh;
  bool mh_accept;
  double adjustment;
  int block;
  misc.push_back({"block"_nm = 0});

  // Initialise matrix for theta
  double theta[n_rungs][n_par];
  for(int i = 0; i < n_rungs; ++i){
    for(int j = 0; j < n_par; ++j){
      theta[i][j] = theta_init[j];
    }
  }
  // Initialise vector for proposal theta
  // Proposal theta is always the theta given to the likelihood function
  // and therefore must be a named vectord, which is why it is writeable::doubles
  writable::doubles theta_prop(n_par);
  for(int p = 0; p < n_par; ++p){
    theta_prop[p] = theta_init[p];
  }
  theta_prop.names() = theta_names;

  // Initialise value for transformed theta: phi
  double phi[n_rungs][n_par];
  for(int i = 0; i < n_rungs; ++i){
    for(int j = 0; j < n_par; ++j){
      phi[i][j] = theta_to_phi(theta[i][j], transform_type[j], theta_min[j], theta_max[j]);
    }
  }
  // Initialise vector for proposal phi
  std::vector<double> phi_prop(n_par);

  // Initialise vector to store blocked log likelihood
  double block_ll[n_rungs][n_unique_blocks];
    for(int b = 0; b < n_unique_blocks; ++b) {
      for(int i = 0; i < n_rungs; ++i){
      misc["block"] = as_sexp(b + 1);
      block_ll[i][b] = ll_f(theta_prop, data, misc);
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
    lp[r] = lp_f(theta_prop);
  }
  // Initialise proposal log prior
  double lp_prop;
  //////////////////////////////////////////////////////////////////////////////

  // Outputs ///////////////////////////////////////////////////////////////////
  // Initialise output matrix
  writable::doubles_matrix<> out(iterations, n_par + 3);
  out(0, 0) = 1;
  for(int p = 1; p < (n_par + 1); ++p){
    out(0, p) =  theta[0][p];
  }
  out(0, n_par + 1) = lp[0];
  out(0, n_par + 2) = ll[0];
  //////////////////////////////////////////////////////////////////////////////

  // Tuning ////////////////////////////////////////////////////////////////////
  // Initialise matrix for proposal sd
  double proposal_sd[n_par][n_rungs];
  for(int i = 0; i < n_par; ++i){
    for(int j = 0; j < n_rungs; ++j){
      proposal_sd[i][j] = 0.1;
    }
  }

  // Initialise acceptance count matrix
  int acceptance[n_par][n_rungs];
  for(int i = 0; i < n_par; ++i){
    for(int j = 0; j < n_rungs; ++j){
      acceptance[i][j] = 0;
    }
  }

  // Initialise swap acceptance count vector
  std::vector<int> swap_acceptance(n_rungs - 1);
  for(int i = 0; i < (n_rungs - 1); ++i){
    swap_acceptance[i] = 0;
  }
  //////////////////////////////////////////////////////////////////////////////

  // PT ////////////////////////////////////////////////////////////////////////
  int rung;
  double rung_beta;
  double rung_proposal_sd;
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
        // Propose new value of phi
        phi_prop[p] = Rf_rnorm(phi[index][p], proposal_sd[p][r]);
        // Transform new phi_prop -> theta_prop
        theta_prop[p] = phi_to_theta(phi_prop[p], transform_type[p], theta_min[p], theta_max[p]);

        // Copy blocked log-likelihood
        for(int b = 0; b < n_unique_blocks; ++b){
          block_ll_prop[b] = block_ll[index][b];
        }
        // For any block this parameter is associated with, update blocked log-likelihood
        writable::integers blocks(blocks_list[p]);
        for(int b = 0; b < blocks.size(); ++b) {
          int block = blocks[b];
          misc["block"] = as_sexp(block);
          block_ll_prop[block - 1] = ll_f(theta_prop, data, misc);
        }
        // Update proposal prior
        lp_prop = lp_f(theta_prop);

        // get parameter transformation adjustment
        adjustment = get_adjustment(theta[index][p], theta_prop[p], transform_type[p], theta_min[p], theta_max[p]);
        // calculate Metropolis-Hastings ratio
        mh = rung_beta * (sum(block_ll_prop) - ll[r]) + (lp_prop - lp[r]) + adjustment;
        // accept or reject move
        mh_accept = log(Rf_runif(0, 1)) < mh;
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
          if(i <= burnin){
            proposal_sd[p][r] = exp(log(proposal_sd[p][r]) + (1 - target_acceptance) / sqrt(i));
          }
          acceptance[p][r] = acceptance[p][r] + 1;
        } else {
          // Revert theta prop
          theta_prop[p] =  theta[index][p];
          // Revert phi prop
          phi_prop[p] = phi[index][p];
          // Robbins monroe step
          if(i <= burnin){
            proposal_sd[p][r] = exp(log(proposal_sd[p][r]) - target_acceptance / sqrt(i));
          }
        }
      }
      // Only store values for the cold chain
      if(r == 0){
        out(i, 0) = i;
        // Record parameters
        for(int p = 0; p < n_par; ++p){
          out(i, p + 1) =  theta[index][p];
        }
        out(i, n_par + 1) = lp[0];
        out(i, n_par + 2) = ll[0];
      }
    }

    // loop over rungs, starting with the hottest chain and moving to the cold
    // chain. Each time propose a swap with the next rung up
    if(swap){
      for(int r = (n_rungs - 1); r >0; --r){
        double rung_beta1 = beta[r];
        double rung_beta2 = beta[r - 1];
        double loglike1 = ll[r];
        double loglike2 = ll[r - 1];

        double acceptance = (loglike2*rung_beta1 + loglike1*rung_beta2) - (loglike1*rung_beta1 + loglike2*rung_beta2);
        // accept or reject move
        bool accept_move = log(Rf_runif(0, 1)) < acceptance;

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
  }

  // Extract the cold chain proposal sd and acceptance rates for output
  std::vector<double> proposal_sd_out(n_par);
  std::vector<double> acceptance_out(n_par);
  for(int i = 0; i < n_par; ++i){
    proposal_sd_out[i] = proposal_sd[i][0];
    acceptance_out[i] = acceptance[i][0];
  }

  // Return outputs in a list
  return writable::list({
      "output"_nm = out,
      "proposal_sd"_nm = proposal_sd_out,
      "acceptance"_nm = acceptance_out,
      "rung_index"_nm = rung_index,
      "swap_acceptance"_nm = swap_acceptance
  });
}
