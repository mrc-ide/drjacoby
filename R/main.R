
#------------------------------------------------
#' @title Check that drjacoby package has loaded successfully
#'
#' @description Simple function to check that drjacoby package has loaded 
#'   successfully. Prints "drjacoby loaded successfully!" if so.
#'
#' @export

check_drjacoby_loaded <- function() {
  message("drjacoby loaded successfully!")
}

#------------------------------------------------
#' @title Run drjacoby MCMC
#'
#' @description Run flexible MCMC through drjacoby. 
#'
#' @param data a vector of data values. When using C++ versions of the
#'   likelihood and/or prior these values are treated internally as doubles, so
#'   while integer and boolean values can be used, keep in mind that these will
#'   be recast as doubles in the likelihood (i.e. \code{TRUE = 1.0}).
#' @param df_params a dataframe of parameters. Must contain the following
#'   elements:
#'   \itemize{
#'     \item \code{name} - the parameter name.
#'     \item \code{min} - the minimum value of the parameter. \code{-Inf} is
#'     allowed.
#'     \item \code{max} - the maximum value of the parameter. \code{Inf} is
#'     allowed.
#'     \item \code{init} - the initial value of the parameter.
#'   }
#' @param loglike TODO.
#' @param logprior TODO.
#' @param burnin the number of burn-in iterations.
#' @param samples the number of sampling iterations.
#' @param rungs the number of temperature rungs used in Metropolis coupling (see
#'   \code{coupling_on}).
#' @param chains the number of independent replicates of the MCMC chain to run.
#'   If a \code{cluster} object is defined then these chains are run in
#'   parallel, otherwise they are run in serial.
#' @param coupling_on whether to implement Metropolis-coupling over temperature 
#'   rungs.
#' @param GTI_pow the power used in the generalised thermodynamic integration 
#'   method.
#' @param cluster option to pass in a cluster environment, allowing chains to be
#'   run in parallel (see package "parallel").
#' @param pb_markdown whether to run progress bars in markdown mode, meaning
#'   they are only updated when they reach 100% to avoid large amounts of output
#'   being printed to markdown files.
#' @param silent whether to suppress all console output.
#'
#' @importFrom utils txtProgressBar
#' @export

run_mcmc <- function(data,
                     df_params,
                     loglike,
                     logprior,
                     burnin = 1e3,
                     samples = 1e4,
                     rungs = 1,
                     chains = 1,
                     coupling_on = TRUE,
                     GTI_pow = 3,
                     cluster = NULL,
                     pb_markdown = FALSE,
                     silent = FALSE) {
  
  # Cleanup pointers on exit
  on.exit(gc())
  
  # ---------- check inputs ----------
  # check data
  assert_vector(data)
  assert_numeric(data)
  
  # check df_params
  assert_dataframe(df_params)
  assert_in(c("name", "min", "max", "init"), names(df_params),
            message = "df_params must contain the columns 'name', 'min', 'max', 'init'")
  assert_numeric(df_params$min)
  assert_numeric(df_params$max)
  assert_le(df_params$min, df_params$max)
  assert_numeric(df_params$init)
  assert_gr(df_params$init, df_params$min)
  assert_le(df_params$init, df_params$max)
  
  # check loglikelihood and logprior functions
  assert_custom_class(loglike, c("function", "character"))
  assert_custom_class(logprior, c("function", "character"))
  # TODO - further checks that these functions are defined correctly?
  
  # check MCMC parameters
  assert_single_pos_int(burnin, zero_allowed = FALSE)
  assert_single_pos_int(samples, zero_allowed = FALSE)
  assert_single_pos_int(rungs, zero_allowed = FALSE)
  assert_single_pos_int(chains, zero_allowed = FALSE)
  assert_single_logical(coupling_on)
  assert_single_pos(GTI_pow, zero_allowed = FALSE)
  
  # check misc parameters
  if (!is.null(cluster)) {
    assert_custom_class(cluster, "cluster")
  }
  assert_single_logical(pb_markdown)
  assert_single_logical(silent)
  
  
  # ---------- pre-processing ----------
  
  # calculate transformation type for each parameter
  # 0 = [-Inf,Inf] -> phi = theta
  # 1 = [-Inf,b]   -> phi = log(b - theta)
  # 2 = [a,Inf]    -> phi = log(theta - a)
  # 3 = [a,b]      -> phi = log((theta - a)/(b - theta))
  df_params$trans_type <- 2*is.finite(df_params$min) + is.finite(df_params$max)
  
  # flag whether likelihood/prior are C++ functions
  loglike_use_cpp <- inherits(loglike, "character")
  logprior_use_cpp <- inherits(logprior, "character")
  
  
  # ---------- define argument lists ----------
  
  # parameters to pass to C++
  args_params <- list(x = data,
                      loglike_use_cpp = loglike_use_cpp,
                      logprior_use_cpp = logprior_use_cpp,
                      theta_init = df_params$init,
                      theta_min = df_params$min,
                      theta_max = df_params$max,
                      trans_type = df_params$trans_type,
                      burnin = burnin,
                      samples = samples,
                      rungs = rungs,
                      coupling_on = coupling_on,
                      GTI_pow = GTI_pow,
                      pb_markdown = pb_markdown,
                      silent = silent)
  
  # functions to pass to C++
  args_functions <- list(loglike = loglike,
                         logprior = logprior,
                         test_convergence = test_convergence,
                         update_progress = update_progress)
  
  # complete list of arguments
  args <- list(args_params = args_params,
               args_functions = args_functions)
  
  # replicate arguments over chains
  parallel_args <- replicate(chains, args, simplify = FALSE)
  for (i in 1:chains) {
    parallel_args[[i]]$args_params$chain <- i
  }
  
  
  # ---------- run MCMC ----------
  
  # split into parallel and serial implementations
  if (!is.null(cluster)) {
    
    # run in parallel
    parallel::clusterEvalQ(cluster, library(drjacoby))
    output_raw <- parallel::clusterApplyLB(cl = cluster, parallel_args, deploy_chain)
    
  } else {
    
    # run in serial
    output_raw <- lapply(parallel_args, deploy_chain)
  }
  
  # ---------- process output ----------
  
  # define names
  chain_names <- sprintf("chain%s", 1:chains)
  rung_names <- sprintf("rung%s", 1:rungs)
  param_names <- df_params$name
  
  # define processed output list
  output_processed <- replicate(chains, list())
  names(output_processed) <- chain_names
  
  # loop through chains
  for (c in 1:chains) {
    
    # get loglikelihoods
    loglike_burnin <- as.data.frame(t(rcpp_to_matrix(output_raw[[c]]$loglike_burnin)))
    loglike_sampling <- as.data.frame(t(rcpp_to_matrix(output_raw[[c]]$loglike_sampling)))
    names(loglike_burnin) <- names(loglike_sampling) <- rung_names
    
    # function for extracting theta into list of matrices over rungs
    get_theta_rungs <- function(theta_list) {
      ret <- mapply(function(x) {
        ret <- as.data.frame(rcpp_to_matrix(x))
        names(ret) <- param_names
        ret
      }, theta_list, SIMPLIFY = FALSE)
      names(ret) <- rung_names
      ret
    }
    
    # get theta values
    theta_burnin <- get_theta_rungs(output_raw[[c]]$theta_burnin)
    theta_sampling <- get_theta_rungs(output_raw[[c]]$theta_sampling)
    names(theta_burnin) <- names(theta_sampling) <- rung_names
    
    # get Metropolis coupling acceptance rates
    beta_raised_vec <- output_raw[[c]]$beta_raised_vec
    mc_accept_burnin <- output_raw[[c]]$mc_accept_burnin/burnin
    mc_accept_sampling <- output_raw[[c]]$mc_accept_sampling/samples

    # store in processed output list
    output_processed[[c]]$diagnostics <- list(beta_raised = beta_raised_vec,
                                              mc_accept_burnin = mc_accept_burnin,
                                              mc_accept_sampling = mc_accept_sampling,
                                              ess = apply(theta_sampling$rung1, 2, ess))
    output_processed[[c]]$loglike_burnin <- loglike_burnin
    output_processed[[c]]$loglike_sampling <- loglike_sampling
    output_processed[[c]]$theta_burnin <- theta_burnin
    output_processed[[c]]$theta_sampling <- theta_sampling
  }
  
  # Estimate rhat
  if(chains > 1){
    rhat_est <- set_rhat(output_processed, chains)
    # Add rhat etimate to each chain diagnostic element
    for (c in 1:chains) {
      output_processed[[c]]$diagnostics$rhat <- rhat_est
    }
  }
  
  # save output as custom class
  class(output_processed) <- "drjacoby_output"
  
  # return
  return(output_processed)
}

#------------------------------------------------
# deploy main_mcmc for this chain
#' @noRd
deploy_chain <- function(args) {
  
  # convert C++ functions to pointers
  if (args$args_params$loglike_use_cpp) {
    args$args_functions$loglike <- RcppXPtrUtils::cppXPtr(args$args_functions$loglike)
    RcppXPtrUtils::checkXPtr(args$args_functions$loglike, "SEXP", c("std::vector<double>",
                                                                    "std::vector<double>"))
  }
  if (args$args_params$logprior_use_cpp) {
    args$args_functions$logprior <- RcppXPtrUtils::cppXPtr(args$args_functions$logprior)
    RcppXPtrUtils::checkXPtr(args$args_functions$logprior, "SEXP", "std::vector<double>")
  }
  
  # get parameters
  burnin <- args$args_params$burnin
  samples <- args$args_params$samples
  
  # make progress bars
  pb_burnin <- txtProgressBar(min = 0, max = burnin, initial = NA, style = 3)
  pb_samples <- txtProgressBar(min = 0, max = samples, initial = NA, style = 3)
  args$args_progress <- list(pb_burnin = pb_burnin,
                             pb_samples = pb_samples)
  
  
  # run C++ function
  ret <- main_cpp(args)
  
  # remove arguments
  rm(args)

  return(ret)
}
