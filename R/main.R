
#------------------------------------------------
# link to Rcpp
#' @useDynLib drjacoby, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

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
#' @title Run MCMC
#'
#' @description Run MCMC.
#'
#' @param data TODO.
#' @param df_params a dataframe of parameters. Must contain the following
#'   elements: TODO
#' @param input_type switch between different methods of injecting likelihood
#'   and prior.
#' @param r_loglike TODO.
#' @param r_logprior TODO.
#' @param c_loglike TODO.
#' @param c_logprior TODO.
#' @param cpp_loglike TODO.
#' @param cpp_logprior TODO.
#' @param burnin the number of burn-in iterations (see also \code{burnin_phases}).
#' @param samples the number of sampling iterations.
#' @param rungs the number of temperature rungs used in Metropolis coupling (see
#'   \code{coupling_on}).
#' @param chains the number of independent replicates of the MCMC chain to run.
#'   If a \code{cluster} object is defined then these chains are run in
#'   parallel, otherwise they are run in serial.
#' @param burnin_phases the number of times a burn-in should be carried out
#'   prior to entering the sampling phase. Several parameters can be defined
#'   separately for each phase of the burn-in, giving greater control over which
#'   tuning methods are applied, and in which order. For example, by default the
#'   proposal bandwidth is updated in the first burn-in phase, followed by two
#'   burn-in phases in which both the proposal bandwidth and covariance are
#'   updated. The complete list of arguments that can be specified separately
#'   for each burn-in phase is: \code{burnin}, \code{bw_update},
#'   \code{cov_update}, \code{coupling_on}. Each of these arguments can be
#'   specified as a vector of length \code{burnin_phases}, or as a single value
#'   in which case the same value applies over every phase.
#' @param bw_update whether the proposal bandwidth should be updated dynamically
#'   via Robbins-Monro (see also \code{burnin_phases}).
#' @param cov_update whether the proposal covariance matrix should be updated
#'   dynamically from the posterior draws (see also \code{burnin_phases}).
#' @param coupling_on whether to implement Metropolis-coupling over temperature 
#'   rungs (see also \code{burnin_phases}).
#' @param GTI_pow the power used in the generalised thermodynamic integration 
#'   method. Must be greater than 1.1.
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
                     input_type = 1,
                     r_loglike,
                     r_logprior,
                     c_loglike,
                     c_logprior,
                     cpp_loglike,
                     cpp_logprior,
                     burnin = 1e3,
                     samples = 1e4,
                     rungs = 11,
                     chains = 1,
                     burnin_phases = 1,
                     bw_update = TRUE,
                     cov_update = FALSE,
                     coupling_on = TRUE,
                     GTI_pow = 3,
                     cluster = NULL,
                     pb_markdown = FALSE,
                     silent = FALSE) {
  
  
  # ---------- check inputs ----------
  
  # check df_params
  assert_dataframe(df_params)
  assert_in(c("name", "min", "max", "init"), names(df_params),
            message = "df_frame must contain the columns 'name', 'min', 'max', 'init'")
  assert_numeric(df_params$min)
  assert_numeric(df_params$max)
  assert_le(df_params$min, df_params$max)
  assert_numeric(df_params$init)
  assert_gr(df_params$init, df_params$min)
  assert_le(df_params$init, df_params$max)
  
  # check likelihoods and data
  # TODO - loglike
  # TODO - logprior
  # TODO - data?
  
  # check MCMC parameters
  assert_single_pos_int(burnin_phases, zero_allowed = FALSE)
  assert_pos_int(burnin, zero_allowed = TRUE)
  assert_single_pos_int(samples, zero_allowed = FALSE)
  assert_single_pos_int(chains, zero_allowed = FALSE)
  assert_single_pos_int(rungs, zero_allowed = FALSE)
  assert_logical(bw_update)
  assert_logical(cov_update)
  assert_logical(coupling_on)
  assert_single_pos(GTI_pow)
  assert_gr(GTI_pow, 1.1)
  
  # check misc parameters
  if (!is.null(cluster)) {
    assert_custom_class(cluster, "cluster")
  }
  assert_single_logical(pb_markdown)
  assert_single_logical(silent)
  
  # any arguments that can be defined either as single values or as vectors are
  # forced to vectors by repeating the same single value
  burnin <- force_vector(burnin, burnin_phases)
  bw_update <- force_vector(bw_update, burnin_phases)
  cov_update <- force_vector(cov_update, burnin_phases)
  coupling_on <- force_vector(coupling_on, burnin_phases)
  assert_length(burnin, burnin_phases)
  assert_same_length_multiple(burnin, bw_update, cov_update, coupling_on)
  
  # ---------- pre-processing ----------
  
  # calculate transformation type for each parameter
  # 0 = [-Inf,Inf] -> phi = theta
  # 1 = [-Inf,b]   -> phi = log(b - theta)
  # 2 = [a,Inf]    -> phi = log(theta - a)
  # 3 = [a,b]      -> phi = log((theta - a)/(b - theta))
  df_params$trans_type <- 2*is.finite(df_params$min) + is.finite(df_params$max)
  
  
  # ---------- define argument lists ----------
  
  # parameters to pass to C++
  args_params <- list(x = data,
                      theta_init = df_params$init,
                      theta_min = df_params$min,
                      theta_max = df_params$max,
                      trans_type = df_params$trans_type,
                      burnin = burnin,
                      samples = samples,
                      rungs = rungs,
                      burnin_phases = burnin_phases,
                      bw_update = bw_update,
                      cov_update = cov_update,
                      coupling_on = coupling_on,
                      GTI_pow = GTI_pow,
                      pb_markdown = pb_markdown,
                      silent = silent,
                      input_type = input_type)
  
  # functions to pass to C++
  args_functions <- list(r_loglike = r_loglike,
                         r_logprior = r_logprior,
                         c_loglike = c_loglike,
                         c_logprior = c_logprior,
                         cpp_loglike = cpp_loglike,
                         cpp_logprior = cpp_logprior,
                         test_convergence = test_convergence,
                         update_progress = update_progress)
  
  # progress bars
  pb_burnin <- list()
  for (i in 1:burnin_phases) {
    pb_burnin[[i]] <- txtProgressBar(min = 0, max = burnin[i], initial = NA, style = 3)
  }
  pb_samples <- txtProgressBar(min = 0, max = samples, initial = NA, style = 3)
  args_progress <- list(pb_burnin = pb_burnin,
                        pb_samples = pb_samples)
  
  # complete list of arguments
  args <- list(args_params = args_params,
               args_functions = args_functions,
               args_progress = args_progress)
  
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
    output_raw <- parallel::clusterApplyLB(cl = cluster, parallel_args, run_mcmc_cpp)
    
  } else {
    
    # run in serial
    output_raw <- lapply(parallel_args, run_mcmc_cpp)
  }
  
  # ---------- process output ----------
  
  # define names
  chain_names <- sprintf("chain%s", 1:chains)
  phase_names <- sprintf("phase%s", 1:burnin_phases)
  rung_names <- sprintf("rung%s", 1:rungs)
  param_names <- df_params$name
  
  # define processed output list
  output_processed <- replicate(chains, list())
  names(output_processed) <- chain_names
  
  # loop through chains
  for (c in 1:chains) {
    
    # get burnin loglikelihoods
    loglike_burnin <- mapply(function(x) {
      ret <- as.data.frame(t(rcpp_to_matrix(x)))
      names(ret) <- rung_names
      ret
    }, output_raw[[c]]$loglike_burnin, SIMPLIFY = FALSE)
    names(loglike_burnin) <- phase_names
    
    # get sampling loglikelihoods
    loglike_sampling <- as.data.frame(t(rcpp_to_matrix(output_raw[[c]]$loglike_sampling)))
    names(loglike_sampling) <- rung_names
    
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
    
    # get burnin theta
    theta_burnin <- mapply(get_theta_rungs, output_raw[[c]]$theta_burnin, SIMPLIFY = FALSE)
    names(theta_burnin) <- phase_names
    
    # get sampling theta
    theta_sampling <- get_theta_rungs(output_raw[[c]]$theta_sampling)
    
    # store in processed output list
    output_processed[[c]]$loglike_burnin <- loglike_burnin
    output_processed[[c]]$loglike_sampling <- loglike_sampling
    output_processed[[c]]$theta_burnin <- theta_burnin
    output_processed[[c]]$theta_sampling <- theta_sampling
  }
  
  # TODO - save output as custom class
  
  # return
  return(output_processed)
}
