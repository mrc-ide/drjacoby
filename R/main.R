
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
#' @description Run MCMC using defined data object, likelihood function, prior
#'   function and parameters.
#'
#' @param data a named list of numeric data values. When using C++ versions of
#'   the likelihood and/or prior these values are treated internally as doubles,
#'   so while integer and boolean values can be used, keep in mind that these
#'   will be recast as doubles in the likelihood (i.e. \code{TRUE = 1.0}).
#' @param df_params a dataframe of parameters. Must contain the following
#'   elements:
#'   \itemize{
#'     \item \code{name} - the parameter name.
#'     \item \code{min} - the minimum value of the parameter. \code{-Inf} is
#'     allowed.
#'     \item \code{max} - the maximum value of the parameter. \code{Inf} is
#'     allowed.
#'     \item \code{init} - the initial value of the parameter. If running
#'     multiple chains a list of initial values can be used to specify distinct
#'     values for each chain.
#'   }
#' @param misc optional list object passed to likelihood and prior.
#' @param loglike,logprior the log-likelihood and log-prior functions used in
#'   the MCMC. Can either be passed in as R functions, or as character strings
#'   which are compiled in C++ functions.
#' @param burnin the number of burn-in iterations.
#' @param samples the number of sampling iterations.
#' @param rungs the number of temperature rungs used in Metropolis coupling (see
#'   \code{coupling_on}).
#' @param chains the number of independent replicates of the MCMC chain to run.
#'   If a \code{cluster} object is defined then these chains are run in
#'   parallel, otherwise they are run in serial.
#' @param coupling_on whether to implement Metropolis-coupling over temperature 
#'   rungs.
#' @param GTI_pow values in the temperature ladder are raised to this power.
#'   Provides a convenient way of concentrating rungs towards one end of the
#'   temperature scale.
#' @param beta_manual option to manually define temperature ladder. These values
#'   are raised to the power \code{GTI_pow}, hence you should use \code{GTI_code
#'   = 1} if you want to fix powers exactly. If \code{NULL} then an equal
#'   spacing of length \code{rungs} is used between 0 and 1.
#' @param cluster option to pass in a cluster environment, allowing chains to be
#'   run in parallel (see package "parallel").
#' @param pb_markdown whether to run progress bars in markdown mode, meaning
#'   they are only updated when they reach 100% to avoid large amounts of output
#'   being printed to markdown files.
#' @param silent whether to suppress all console output.
#'
#' @importFrom utils txtProgressBar
#' @importFrom stats setNames
#' @export

run_mcmc <- function(data,
                     df_params,
                     misc = list(),
                     loglike,
                     logprior,
                     burnin = 1e3,
                     samples = 1e4,
                     rungs = 1,
                     chains = 5,
                     coupling_on = TRUE,
                     GTI_pow = 1.0,
                     beta_manual = NULL,
                     cluster = NULL,
                     pb_markdown = FALSE,
                     silent = FALSE) {
  
  # declare variables to avoid "no visible binding" issues
  stage <- rung <- value <- chain <- link <- NULL
  
  # Cleanup pointers on exit
  on.exit(gc())
  
  # ---------- check inputs ----------
  # check data
  assert_list(data)
  if (is.null(names(data)) | any(names(data) == "")) {
    stop("data must be a *named* list")
  }
  assert_numeric(unlist(data))
  
  # check df_params
  assert_dataframe(df_params)
  assert_in(c("name", "min", "max", "init"), names(df_params),
            message = "df_params must contain the columns 'name', 'min', 'max', 'init'")
  assert_numeric(df_params$min)
  assert_numeric(df_params$max)
  assert_leq(df_params$min, df_params$max)
  if (!is.list(df_params$init)) {
    df_params$init <- as.list(df_params$init)
  }
  mapply(function(i) {  # checks on each initial value
    assert_vector_numeric(df_params$init[[i]], message = "all df_params$init must be numeric")
    if (length(df_params$init[[i]]) != 1) {
      assert_length(df_params$init[[i]], chains, message = paste0("must define one df_params$init value per parameter, ",
                                                                  "or alternatively a list of values one for each chain"))
    }
    msg_range <- "all df_params$init must be within specified range"
    assert_greq(df_params$init[[i]], df_params$min[i], message = msg_range)
    assert_leq(df_params$init[[i]], df_params$max[i], message = msg_range)
  }, seq_along(df_params$init))
  
  # check misc
  assert_list(misc)
  
  # check loglikelihood and logprior functions
  assert_custom_class(loglike, c("function", "character"))
  assert_custom_class(logprior, c("function", "character"))
  
  # check MCMC parameters
  assert_single_pos_int(burnin, zero_allowed = FALSE)
  assert_single_pos_int(samples, zero_allowed = FALSE)
  assert_single_pos_int(rungs, zero_allowed = FALSE)
  assert_single_pos_int(chains, zero_allowed = FALSE)
  assert_single_logical(coupling_on)
  assert_single_pos(GTI_pow)
  
  # calculate/check final temperature vector
  if (is.null(beta_manual)) {
    beta_manual <- rev(seq(1, 0, l = rungs))
  }
  rungs <- length(beta_manual)
  assert_vector_bounded(beta_manual)
  assert_increasing(beta_manual)
  assert_eq(beta_manual[rungs], 1.0)
  
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
  
  # flag to skip over fixed parameters
  skip_param <- (df_params$min == df_params$max)
  
  # get initial values into matrix. Rows for parameters, columns for chains
  init_mat <- do.call(rbind, mapply(function(x) {
    if (length(x) == 1) {
      rep(x, chains)
    } else {
      x
    }
  }, df_params$init, SIMPLIFY = FALSE))
  
  # flag whether likelihood/prior are C++ functions
  loglike_use_cpp <- inherits(loglike, "character")
  logprior_use_cpp <- inherits(logprior, "character")
  
  # raise temperature vector to power (prepping for later version which will
  # implement generalised thermodynamic integration)
  beta_raised <- beta_manual^GTI_pow
  
  
  # ---------- define argument lists ----------
  
  # parameters to pass to C++
  args_params <- list(x = data,
                      misc = misc,
                      loglike_use_cpp = loglike_use_cpp,
                      logprior_use_cpp = logprior_use_cpp,
                      theta_min = df_params$min,
                      theta_max = df_params$max,
                      trans_type = df_params$trans_type,
                      skip_param = skip_param,
                      burnin = burnin,
                      samples = samples,
                      rungs = rungs,
                      coupling_on = coupling_on,
                      beta_raised = beta_raised,
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
  
  # create distinct argument sets over chains
  parallel_args <- replicate(chains, args, simplify = FALSE)
  for (i in 1:chains) {
    parallel_args[[i]]$args_params$chain <- i
    
    # create named vector object for passing internally within C++ functions.
    # Initial values defined separately for each chain
    parallel_args[[i]]$args_params$theta_vector <- setNames(init_mat[,i], df_params$name)
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
  
  # get raw output into dataframe
  df_output <- do.call(rbind, mapply(function(j) {
    do.call(rbind, mapply(function(i) {
      
      # concatenate burn-in and sampling logprior and loglikelihood
      logprior <- c(output_raw[[j]]$logprior_burnin[[i]], output_raw[[j]]$logprior_sampling[[i]])
      loglike <- c(output_raw[[j]]$loglike_burnin[[i]], output_raw[[j]]$loglike_sampling[[i]])
      
      # create dataframe of loglike and logprior
      ret <- data.frame(chain = chain_names[j],
                        rung = rung_names[i],
                        iteration = seq_along(loglike),
                        stage = rep(c("burnin", "sampling"), times = c(burnin, samples)),
                        logprior = logprior,
                        loglikelihood = loglike)
      
      # concatenate theta into dataframe and append
      theta <- as.data.frame(do.call(rbind, c(output_raw[[j]]$theta_burnin[[i]], output_raw[[j]]$theta_sampling[[i]])))
      names(theta) <- param_names
      ret <- cbind(ret, theta)
      
      return(ret)
    }, seq_len(rungs), SIMPLIFY = FALSE))
  }, seq_len(chains), SIMPLIFY = FALSE))
  
  # append to output list
  output_processed <- list(output = df_output)
  output_processed$diagnostics <- list()
  
  ## Diagnostics
  # Rhat (Gelman-Rubin diagnostic)
  if (chains > 1) {
    rhat_est <- c()
    for (p in seq_along(param_names)) {
      pm <- subset(output_processed$output, stage == "sampling", select = c("chain", param_names[p]))
      rhat_est[p] <- gelman_rubin(pm, chains, samples)
    }
    rhat_est[skip_param] <- NA
    output_processed$diagnostics$rhat <- rhat_est
  }
  
  # ESS
  # NOTE - some issues with line ess_est <- apply(output_sub, 2, coda::effectiveSize), causing tests to fail. Adding tryCatch line fixes the problem, even though problem line is unchanged. Leaving commented out for now so can proceed with development.
  #output_sub <- subset(output_processed$output, stage == "sampling" & rung == "rung1",
  #                     select = as.character(param_names))
  #tc <- tryCatch(apply(output_sub, 2, coda::effectiveSize))
  #ess_est <- apply(output_sub, 2, coda::effectiveSize)
  #ess_est[skip_param] <- NA
  #output_processed$diagnostics$ess <- ess_est
  
  # Thermodynamic power
  output_processed$diagnostics$rung_details <- data.frame(rung = 1:rungs,
                                                          thermodynamic_power = beta_raised)
  
  # Metropolis-coupling
  # store acceptance rates between pairs of rungs (links)
  mc_accept <- NA
  if (rungs > 1) {
    
    # MC accept
    mc_accept <- expand.grid(link = seq_len(rungs - 1), chain = chain_names)
    mc_accept_burnin <- unlist(lapply(output_raw, function(x){x$mc_accept_burnin})) / burnin
    mc_accept_sampling <- unlist(lapply(output_raw, function(x){x$mc_accept_sampling})) / samples
    mc_accept <- rbind(cbind(mc_accept, stage = "burnin", value = mc_accept_burnin),
                       cbind(mc_accept, stage = "sampling", value = mc_accept_sampling))
    
  }
  output_processed$diagnostics$mc_accept <- mc_accept
  
  ## Parameters
  output_processed$parameters <- list(data = data,
                                      df_params = df_params,
                                      loglike = loglike,
                                      logprior = logprior,
                                      burnin = burnin,
                                      samples = samples,
                                      rungs = rungs,
                                      chains = chains,
                                      coupling_on = coupling_on,
                                      GTI_pow = GTI_pow,
                                      beta_manual = beta_manual)
  
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
    RcppXPtrUtils::checkXPtr(args$args_functions$loglike, "SEXP", c("Rcpp::NumericVector",
                                                                    "int",
                                                                    "Rcpp::List",
                                                                    "Rcpp::List"))
  }
  if (args$args_params$logprior_use_cpp) {
    args$args_functions$logprior <- RcppXPtrUtils::cppXPtr(args$args_functions$logprior)
    RcppXPtrUtils::checkXPtr(args$args_functions$logprior, "SEXP", c("Rcpp::NumericVector",
                                                                     "int",
                                                                     "Rcpp::List"))
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

