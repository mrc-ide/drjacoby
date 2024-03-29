
#include "System.h"
#include "misc.h"

using namespace std;

void System::load(Rcpp::List args) {
  
  // split argument lists
  Rcpp::List args_params = args["args_params"];
  Rcpp::List args_functions = args["args_functions"];
  Rcpp::List args_progress = args["args_progress"];
  Rcpp::List args_progress_burnin = args_progress["pb_burnin"];
  
  // data
  x = args_params["x"];
  
  // misc
  misc = args_params["misc"];
  
  // model parameters
  theta_vector = args_params["theta_vector"];
  theta_min = rcpp_to_vector_double(args_params["theta_min"]);
  theta_max = rcpp_to_vector_double(args_params["theta_max"]);
  block = rcpp_to_matrix_int(args_params["block"]);
  n_block = rcpp_to_int(args_params["n_block"]);
  trans_type = rcpp_to_vector_int(args_params["trans_type"]);
  skip_param = rcpp_to_vector_bool(args_params["skip_param"]);
  d = int(theta_min.size());
  target_acceptance = rcpp_to_double(args_params["target_acceptance"]);
  
  // MCMC parameters
  burnin = rcpp_to_int(args_params["burnin"]);
  samples = rcpp_to_int(args_params["samples"]);
  rungs = rcpp_to_int(args_params["rungs"]);
  coupling_on = rcpp_to_bool(args_params["coupling_on"]);
  beta_raised = rcpp_to_vector_double(args_params["beta_raised"]);
  chain = rcpp_to_int(args_params["chain"]);
  
  // misc parameters
  save_hot_draws = rcpp_to_bool(args_params["save_hot_draws"]);
  pb_markdown = rcpp_to_bool(args_params["pb_markdown"]);
  silent = rcpp_to_bool(args_params["silent"]);
  
}
