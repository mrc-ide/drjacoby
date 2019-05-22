
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
#' @details TODO - explain df_burnin defaults.
#'
#' @param data TODO.
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
#' @param df_burnin a dataframe specifying the parameters of each phase of
#'   burn-in. Must contain the following elements:
#'   \itemize{
#'     \item \code{iterations} - the number of MCMC iterations this phase.
#'     \item \code{prop_method} - integer flag specifying which method to use
#'     when proposing new values this phase. 1 = independent proposal per
#'     parameter, 2 = multivariate proposal over all parameters.
#'     \item \code{bw_update} - boolean, whether to update the proposal
#'     bandwidth dynamically via Robbins-Monro this phase. If \code{TRUE} then
#'     for the independent proposal method this means updating a different
#'     bandwidth per parameter, for the multivariate proposal method this means
#'     updating a single bandwidth for all parameters, scaled by a covariance
#'     matrix.
#'     \item \code{bw_reset} - boolean, whether the Robbins-Monro index should
#'     be reset at the beginning of this phase. Resetting the index makes the
#'     proposal bandwidths more free to move to new values.
#'     \item \code{cov_recalc} boolean, whether the multivariate proposal
#'     covariance matrix should be re-calculated from the posterior draws at the
#'     end of this burn-in phase.
#'   }
#'   If \code{NULL} then a default dataframe is used (see details).
#' @param samples the number of sampling iterations.
#' @param rungs the number of temperature rungs used in Metropolis coupling (see
#'   \code{coupling_on}).
#' @param chains the number of independent replicates of the MCMC chain to run.
#'   If a \code{cluster} object is defined then these chains are run in
#'   parallel, otherwise they are run in serial.
#' @param coupling_on whether to implement Metropolis-coupling over temperature 
#'   rungs.
#' @param GTI_pow the power used in the generalised thermodynamic integration 
#'   method. Must be greater than 1.0.
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
                     df_burnin = NULL,
                     samples = 1e4,
                     rungs = 1,
                     chains = 1,
                     coupling_on = TRUE,
                     GTI_pow = 3,
                     cluster = NULL,
                     pb_markdown = FALSE,
                     silent = FALSE) {
  
  
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
  assert_custom_class(loglike, c("function", "XPtr"))
  # TODO - further checks that these functions are defined correctly?
  
  # define default df_burnin
  if (is.null(df_burnin)) {
    df_burnin <- data.frame(iterations = 1e3,
                            prop_method = c(1, 1, 2),
                            bw_update = TRUE,
                            bw_reset = TRUE,
                            cov_recalc = c(FALSE, TRUE, TRUE))
  }
  
  # check df_burnin
  assert_dataframe(df_burnin)
  assert_in(c("iterations", "prop_method", "bw_update", "bw_reset", "cov_recalc"), names(df_burnin),
            message = "df_burnin must contain the columns 'iterations', 'prop_method', 'bw_update', 'bw_reset', 'cov_recalc'")
  assert_pos_int(df_burnin$iterations, zero_allowed = FALSE)
  assert_in(df_burnin$prop_method, c(1,2))
  assert_logical(df_burnin$bw_update)
  assert_logical(df_burnin$bw_reset)
  assert_logical(df_burnin$cov_recalc)
  
  # check other MCMC parameters
  assert_single_pos_int(samples, zero_allowed = FALSE)
  assert_single_pos_int(rungs, zero_allowed = FALSE)
  assert_single_pos_int(chains, zero_allowed = FALSE)
  assert_single_logical(coupling_on)
  assert_single_pos(GTI_pow)
  assert_gr(GTI_pow, 1.0)
  
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
  loglike_use_cpp <- inherits(loglike, "XPtr")
  logprior_use_cpp <- inherits(logprior, "XPtr")
  
  
  # ---------- define argument lists ----------
  
  # parameters to pass to C++
  args_params <- list(x = data,
                      loglike_use_cpp = loglike_use_cpp,
                      logprior_use_cpp = logprior_use_cpp,
                      theta_init = df_params$init,
                      theta_min = df_params$min,
                      theta_max = df_params$max,
                      trans_type = df_params$trans_type,
                      burnin = df_burnin$iterations,
                      prop_method = df_burnin$prop_method,
                      bw_update = df_burnin$bw_update,
                      bw_reset = df_burnin$bw_reset,
                      cov_recalc = df_burnin$cov_recalc,
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
  
  # progress bars
  pb_burnin <- list()
  for (i in 1:nrow(df_burnin)) {
    pb_burnin[[i]] <- txtProgressBar(min = 0, max = df_burnin$iterations[i], initial = NA, style = 3)
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
    output_raw <- parallel::clusterApplyLB(cl = cluster, parallel_args, main_cpp)
    
  } else {
    
    # run in serial
    output_raw <- lapply(parallel_args, main_cpp)
  }
  
  
  # ---------- process output ----------
  
  # define names
  chain_names <- sprintf("chain%s", 1:chains)
  phase_names <- sprintf("phase%s", 1:nrow(df_burnin))
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
    
    # get burnin phi
    phi_burnin <- mapply(get_theta_rungs, output_raw[[c]]$phi_burnin, SIMPLIFY = FALSE)
    names(phi_burnin) <- phase_names
    
    # get sampling theta
    theta_sampling <- get_theta_rungs(output_raw[[c]]$theta_sampling)
    
    # store in processed output list
    output_processed[[c]]$loglike_burnin <- loglike_burnin
    output_processed[[c]]$loglike_sampling <- loglike_sampling
    output_processed[[c]]$theta_burnin <- theta_burnin
    output_processed[[c]]$phi_burnin <- phi_burnin
    output_processed[[c]]$theta_sampling <- theta_sampling
  }
  
  # TODO - save output as custom class
  
  # return
  return(output_processed)
}
