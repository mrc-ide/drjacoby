
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
#' @param df_params a dataframe of parameters.
#' @param loglike TODO.
#' @param logprior TODO.
#' @param data TODO.
#' @param burnin the number of burn-in iterations.
#' @param samples the number of sampling iterations.
#' @param chains the number of independent replicates of the MCMC chain to run.
#'   If a \code{cluster} object is defined then these chains are run in
#'   parallel, otherwise they are run in serial.
#' @param rungs the number of temperature rungs.
#' @param GTI_pow the power used in the generalised thermodynamic integration 
#'   method. Must be greater than 1.1.
#' @param autoconverge_on whether convergence should be assessed automatically
#'   every \code{converge_test} iterations, leading to termination of the
#'   burn-in phase. If \code{FALSE} then the full \code{burnin} iterations are
#'   used.
#' @param converge_test test for convergence every \code{convergence_test} 
#'   iterations if \code{auto_converge} is being used.
#' @param coupling_on whether to implement Metropolis-coupling over temperature 
#'   rungs.
#' @param cluster option to pass in a cluster environment, allowing chains to be
#'   run in parallel (see package "parallel").
#' @param pb_markdown whether to run progress bars in markdown mode, meaning
#'   they are only updated when they reach 100% to avoid large amounts of output
#'   being printed to markdown files.
#' @param silent whether to suppress all console output.
#'
#' @importFrom utils txtProgressBar
#' @export

run_mcmc <- function(df_params,
                     loglike,
                     logprior,
                     data,
                     burnin = 1e3,
                     samples = 1e4,
                     chains = 1,
                     rungs = 11,
                     GTI_pow = 3,
                     autoconverge_on = TRUE,
                     converge_test = 1e2,
                     coupling_on = FALSE,
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
  
  # check other inputs
  # TODO - loglike
  # TODO - logprior
  # TODO - data?
  assert_single_pos_int(burnin, zero_allowed = FALSE)
  assert_single_pos_int(samples, zero_allowed = FALSE)
  assert_single_pos_int(chains, zero_allowed = FALSE)
  assert_single_pos_int(rungs, zero_allowed = FALSE)
  assert_single_pos(GTI_pow)
  assert_gr(GTI_pow, 1.1)
  assert_single_logical(autoconverge_on)
  assert_single_logical(coupling_on)
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
  
  
  # ---------- define argument lists ----------
  
  # parameters to pass to C++
  args_params <- list(x = x,
                      theta_init = df_params$init,
                      theta_min = df_params$min,
                      theta_max = df_params$max,
                      trans_type = df_params$trans_type,
                      burnin = burnin,
                      samples = samples,
                      rungs = rungs,
                      autoconverge_on = autoconverge_on,
                      coupling_on = coupling_on,
                      pb_markdown = pb_markdown,
                      silent = silent)
  
  # functions to pass to C++
  args_functions <- list(loglike = loglike,
                         logprior = logprior,
                         test_convergence = test_convergence,
                         update_progress = update_progress)
  
  # progress bars
  pb_burnin <- txtProgressBar(min = 0, max = burnin, initial = NA, style = 3)
  pb_samples <- txtProgressBar(min = 0, max = samples, initial = NA, style = 3)
  args_progress <- list(pb_burnin = pb_burnin,
                        pb_samples = pb_samples)
  
  # complete list of arguments
  args <- list(args_params = args_params,
               args_functions = args_functions,
               args_progress = args_progress)
  
  # replicate arguments over chains
  parallel_args <- replicate(chains, args, simplify = FALSE)
  
  
  # ---------- run MCMC ----------
  
  # split into parallel and serial implementations
  if (!is.null(cluster)) {
    
    # run in parallel
    clusterEvalQ(cluster, library(drjacoby))
    output_raw <- clusterApplyLB(cl = cluster, parallel_args, run_mcmc_cpp)
    
  } else {
    
    # run in serial
    output_raw <- lapply(parallel_args, run_mcmc_cpp)
  }
  
  # ---------- process output ----------
  
  # define processed output list
  output_processed <- replicate(chains, list())
  names(output_processed) <- sprintf("chain%s", 1:chains)
  
  # loop through chains
  for (c in 1:chains) {
    
    # get theta values
    output_processed[[c]]$theta <- rcpp_to_matrix(output_raw[[c]]$theta)
  }
  
  # TODO - save output as custom class
  
  # return
  return(output_processed)
}
