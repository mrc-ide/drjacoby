
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
#' @title Define parameters dataframe
#'
#' @description Provides a convenient way of defining parameters in the format
#'   required by \code{run_mcmc()}. Each parameter must have the following three
#'   elements, defined in order:
#'   \itemize{
#'     \item \code{name} - the parameter name.
#'     \item \code{min} - the minimum value of the parameter. \code{-Inf} is
#'     allowed.
#'     \item \code{max} - the maximum value of the parameter. \code{Inf} is
#'     allowed.
#'   }
#'   There following arguments are also optional:
#'   \itemize{
#'     \item \code{init} - the initial value of the parameter. If running
#'     multiple chains a vector of initial values can be used to specify distinct
#'     values for each chain.
#'     \item \code{block} - which likelihood block(s) this parameter belongs to.
#'     See vignettes for instructions on using likelihood blocks.
#'   }
#'
#' @param ... a series of named input arguments.
#'   
#'
#' @export
#' @examples
#' define_params(name = "mu", min = -10, max = 10, init = 0,
#'               name = "sigma", min = 0, max = 5, init = c(1, 2))
#'               
#' define_params(name = "mu1", min = -10, max = 10, init = 0, block = 1,
#'               name = "mu2", min = -10, max = 10, init = 0, block = 2,
#'               name = "sigma", min = 0, max = 5, init = 1, block = c(1, 2))

define_params <- function(...) {
  x <- list(...)
  
  # check input format of arguments
  assert_gr(length(x), 0, message = "input cannot be empty")
  assert_in(names(x), c("name", "min", "max", "init", "block"))
  use_init <- ("init" %in% names(x))
  use_block <- ("block" %in% names(x))
  n_cols <- 3 + use_init + use_block
  if ((length(x) %% n_cols) != 0) {
    stop("must have the same number of inputs per parameter")
  }
  n_param <- length(x) / n_cols
  arg_names <-c("name", "min", "max")
  if (use_init) {
    arg_names <- c(arg_names, "init")
  }
  if (use_block) {
    arg_names <- c(arg_names, "block")
  }
  assert_eq(names(x), rep(arg_names, n_param))
  
  # create params dataframe
  v <- n_cols*(0:(n_param - 1))
  ret <- data.frame(name = unlist(x[1 + v]),
                    min = unlist(x[2 + v]),
                    max = unlist(x[3 + v]))
  if (use_init) {
    ret$init <- x[which(arg_names == "init") + v]
  }
  if (use_block) {
    ret$block <- x[which(arg_names == "block") + v]
  }
  
  # run checks and standardise format
  ret <- check_params(ret)
  
  return(ret)
}

#------------------------------------------------
# Check that params dataframe is formatted correctly, and return in standardised
# format (init and block coerced to list)
#' @noRd
check_params <- function(x) {
  
  # check dataframe has correct elements
  assert_dataframe(x)
  assert_in(c("name", "min", "max"), names(x),
            message = "df_params must contain the columns 'name', 'min', 'max'")
  if (any(duplicated(x$name))) {
    stop("parameter names must be unique")
  }
  use_init <- ("init" %in% names(x))
  use_block <- ("block" %in% names(x))
  
  # coerce init and block to list
  if (use_init) {
    if (!is.list(x$init)) {
      x$init <- as.list(x$init)
    }
  }
  if (use_block) {
    if (!is.list(x$block)) {
      x$block <- as.list(x$block)
    }
  }
  
  # check each row in turn
  for (i in seq_len(nrow(x))) {
    
    # check format
    assert_single_string(x$name[i], message = "parameter names must be character strings")
    assert_single_numeric(x$min[i], message = "min values must be single values")
    assert_single_numeric(x$max[i], message = "min values must be single values")
    if (use_init) {
      assert_vector_numeric(x$init[[i]], message = "init values must be numeric")
    }
    if (use_block) {
      assert_vector_numeric(x$block[[i]], message = "block values must be numeric")
    }
    
    # check order
    assert_leq(x$min[i], x$max[i], message = "min values must be less than or equal to max values")
    if (use_init) {
      assert_greq(x$init[[i]], x$min[i], message = "init values must be greater than or equal to min values")
      assert_leq(x$init[[i]], x$max[i], message = "init values must be less than or equal to max values")
    }
  }
  
  return(x)
}

#------------------------------------------------
#' @title Run drjacoby MCMC
#'
#' @description Run MCMC either with or without parallel tempering turned on.
#'   Minimum inputs include a data object, a data.frame of parameters, a
#'   log-likelihood function and a log-prior function. Produces an object of
#'   class \code{drjacoby_output}, which contains all MCMC output along with
#'   some diagnostics and a record of inputs.
#'   
#' @details Note that both \code{data} and \code{misc} are passed into
#'   log-likelihood/log-prior functions *by reference*. This means if you modify
#'   these objects inside the functions then any changes will persist.
#'
#' @param data a named list of numeric data values. When using C++ likelihood
#'   and/or prior these values are treated internally as doubles, so while
#'   integer and Boolean values can be used, keep in mind that these will be
#'   recast as doubles in the likelihood (i.e. \code{TRUE = 1.0}).
#' @param df_params a data.frame of parameters (see \code{?define_params}).
#' @param misc optional list object passed to likelihood and prior. This can be
#'   useful for passing values that are not strictly data, for example passing a
#'   lookup table to the prior function.
#' @param loglike,logprior the log-likelihood and log-prior functions used in
#'   the MCMC. Can either be passed in as R functions (not in quotes), or as
#'   character strings naming compiled C++ functions (in quotes).
#' @param burnin the number of burn-in iterations. Automatic tuning of proposal
#'   standard deviations is only active during the burn-in period.
#' @param samples the number of sampling iterations.
#' @param rungs the number of temperature rungs used in the parallel tempering
#'   method. By default, \eqn{\beta} values are equally spaced between 0 and 1,
#'   i.e. \eqn{\beta[i]=}\code{(i-1)/(rungs-1)} for \code{i} in \code{1:rungs}.
#'   The likelihood for the \out{i<sup>th</sup>} heated chain is raised to the
#'   power \eqn{\beta[i]^\alpha}, meaning we can use the \eqn{\alpha} parameter
#'   to concentrate rungs towards the start or the end of the interval (see the
#'   \code{alpha} argument).
#' @param chains the number of independent replicates of the MCMC to run. If a
#'   \code{cluster} object is defined then these chains are run in parallel,
#'   otherwise they are run in serial.
#' @param beta_manual vector of manually defined \eqn{\beta} values used in the
#'   parallel tempering approach. If defined, this overrides the spacing defined
#'   by \code{rungs}. Note that even manually defined \eqn{\beta} values are
#'   raised to the power \eqn{\alpha} internally, hence you should set
#'   \code{alpha = 1} if you want to fix \eqn{\beta} values exactly.
#' @param alpha the likelihood for the \out{i<sup>th</sup>} heated chain is
#'   raised to the power \eqn{\beta[i]^\alpha}, meaning we can use the
#'   \eqn{\alpha} parameter to concentrate rungs towards the start or the end of
#'   the temperature scale.
#' @param target_acceptance Target acceptance rate. Should be between 0 and 1.
#'   Default of 0.44, set as optimum for unvariate proposal distributions.
#' @param cluster option to pass in a cluster environment, allowing chains to be
#'   run in parallel (see package "parallel").
#' @param coupling_on whether to implement Metropolis-coupling over temperature
#'   rungs. The option of deactivating coupling has been retained for general
#'   interest and debugging purposes only. If this parameter is \code{FALSE}
#'   then parallel tempering will have no impact on MCMC mixing.
#' @param pb_markdown whether to run progress bars in markdown mode, meaning
#'   they are only updated when they reach 100\% to avoid large amounts of output
#'   being printed to markdown files.
#' @param save_data if \code{TRUE} (the default) the raw input data is stored
#'   for reference in the project output. This allows complete reproducibility
#'   from a project, but may be undesirable when datasets are very large.
#' @param save_hot_draws if \code{TRUE} the parameter draws relating to the hot
#'   chains are also stored inside the \code{pt} element of the project output.
#'   If \code{FALSE} (the default) only log-likelihoods and log-priors are
#'   stored from heated chains.
#' @param silent whether to suppress all console output.
#'
#' @importFrom utils txtProgressBar
#' @importFrom stats setNames var runif
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
                     beta_manual = NULL,
                     alpha = 1.0,
                     target_acceptance = 0.44,
                     cluster = NULL,
                     coupling_on = TRUE,
                     pb_markdown = FALSE,
                     save_data = TRUE,
                     save_hot_draws = FALSE,
                     silent = FALSE) {
  
  # declare variables to avoid "no visible binding" issues
  phase <- rung <- value <- chain <- link <- NULL
  
  # Cleanup pointers on exit
  on.exit(gc())
  
  # ---------- check inputs ----------
  
  # check data
  assert_list_named(data)
  assert_numeric(unlist(data))
  
  # check misc
  assert_list(misc)
  
  # check loglikelihood and logprior functions
  assert_class(loglike, c("function", "character"))
  assert_class(logprior, c("function", "character"))
  
  # check MCMC parameters
  assert_single_pos_int(burnin, zero_allowed = FALSE)
  assert_single_pos_int(samples, zero_allowed = FALSE)
  assert_single_pos_int(rungs, zero_allowed = FALSE)
  assert_single_pos_int(chains, zero_allowed = FALSE)
  assert_single_logical(coupling_on)
  assert_single_pos(alpha)
  assert_bounded(target_acceptance, 0, 1)
  
  # check df_params
  df_params <- check_params(df_params)
  use_init <- ("init" %in% names(df_params))
  use_block <- ("block" %in% names(df_params))
  if (use_init) {
    for (i in 1:nrow(df_params)) {
      if (length(df_params$init[[i]]) != 1) {
        assert_length(df_params$init[[i]], chains, message = paste0("must define one df_params$init value per parameter, ",
                                                                    "or alternatively a list of values one for each chain"))
      }
    }
  }
  
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
    assert_class(cluster, "cluster")
  }
  assert_single_logical(pb_markdown)
  assert_single_logical(save_data)
  assert_single_logical(save_hot_draws)
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
  
  # define default init values
  if (!use_init) {
    init_list <- list()
    for (i in 1:nrow(df_params)) {
      p <- runif(chains)
      if (df_params$trans_type[i] == 0) {
        init_list[[i]] <- log(p) - log(1 - p)
      } else if (df_params$trans_type[i] == 1) {
        init_list[[i]] <- log(p) + df_params$max[i]
      } else if (df_params$trans_type[i] == 2) {
        init_list[[i]] <- df_params$min[i] - log(p)
      } else if (df_params$trans_type[i] == 3) {
        init_list[[i]] <- df_params$min[i] + (df_params$max[i] - df_params$min[i])*p
      }
    }
    df_params$init <- init_list
  }
  
  # define default blocks
  if (!use_block) {
    df_params$block <- as.list(rep(1, nrow(df_params)))
  }
  
  # get initial values into matrix. Rows for parameters, columns for chains
  init_mat <- do.call(rbind, mapply(function(x) {
    if (length(x) == 1) {
      rep(x, chains)
    } else {
      x
    }
  }, df_params$init, SIMPLIFY = FALSE))
  
  # flag whether likelihood and/or prior are C++ functions
  loglike_use_cpp <- inherits(loglike, "character")
  logprior_use_cpp <- inherits(logprior, "character")
  
  # raise temperature vector to power
  beta_raised <- beta_manual^alpha
  
  # make sure "block" is not an element of misc already being used, and if not
  # create dummy element for storing current block
  if (length(misc) > 0) {
    assert_not_in("block", names(misc), message = "'block' is a reserved name within misc object")
  }
  misc$block <- -1
  
  
  # ---------- define argument lists ----------
  
  # parameters to pass to C++
  args_params <- list(x = data,
                      misc = misc,
                      loglike_use_cpp = loglike_use_cpp,
                      logprior_use_cpp = logprior_use_cpp,
                      theta_min = df_params$min,
                      theta_max = df_params$max,
                      block = df_params$block,
                      n_block = max(unlist(df_params$block)),
                      trans_type = df_params$trans_type,
                      skip_param = skip_param,
                      burnin = burnin,
                      samples = samples,
                      rungs = rungs,
                      coupling_on = coupling_on,
                      beta_raised = beta_raised,
                      save_hot_draws = save_hot_draws,
                      pb_markdown = pb_markdown,
                      silent = silent,
                      target_acceptance = target_acceptance)
  
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
  
  # print total runtime
  chain_runtimes <- mapply(function(x) x$t_diff, output_raw)
  if (!silent) {
    message(sprintf("total MCMC run-time: %s seconds", signif(sum(chain_runtimes), 3)))
  }
  
  
  # ---------- process output ----------
  
  # define names
  chain_names <- 1:chains
  rung_names <- 1:rungs
  param_names <- df_params$name
  
  # get parameter draws into dataframe. This will be over all rungs if
  # save_hot_draws is TRUE, otherwise it will only be over the cold chain
  df_theta <- do.call(rbind, mapply(function(j) {
    do.call(rbind, mapply(function(i) {
      
      theta_burnin <- do.call(rbind, output_raw[[j]]$theta_burnin[[i]]) %>%
        as.data.frame() %>%
        magrittr::set_colnames(param_names) %>%
        dplyr::mutate(chain = chain_names[j],
                      rung = rung_names[i],
                      phase = "burnin", .before = 1)
      
      theta_sampling <- do.call(rbind, output_raw[[j]]$theta_sampling[[i]]) %>%
        as.data.frame() %>%
        magrittr::set_colnames(param_names) %>%
        dplyr::mutate(chain = chain_names[j],
                      rung = rung_names[i],
                      phase = "sampling", .before = 1)
      
      ret <- theta_burnin %>%
        dplyr::bind_rows(theta_sampling) %>%
        dplyr::mutate(iteration = seq_along(phase), .after = "phase")
      
      return(ret)
    }, seq_along(output_raw[[j]]$theta_burnin), SIMPLIFY = FALSE))
  }, seq_along(output_raw), SIMPLIFY = FALSE))
  
  # fix rungs field if save_hot_draws is FALSE
  if (!save_hot_draws) {
    df_theta$rung <- rungs
  }
  
  # get likelihoods and priors over all rungs
  df_pt <- do.call(rbind, mapply(function(j) {
    do.call(rbind, mapply(function(i) {
      
      pt_burnin <- data.frame(chain = chain_names[j],
                              rung = rung_names[i],
                              phase = "burnin",
                              logprior = output_raw[[j]]$logprior_burnin[[i]],
                              loglikelihood = output_raw[[j]]$loglike_burnin[[i]])
      
      pt_sampling <- data.frame(chain = chain_names[j],
                                rung = rung_names[i],
                                phase = "sampling",
                                logprior = output_raw[[j]]$logprior_sampling[[i]],
                                loglikelihood = output_raw[[j]]$loglike_sampling[[i]])
      
      ret <- pt_burnin %>%
        dplyr::bind_rows(pt_sampling) %>%
        dplyr::mutate(iteration = seq_along(phase), .after = "phase")
      
      return(ret)
    }, seq_along(output_raw[[j]]$logprior_burnin), SIMPLIFY = FALSE))
  }, seq_along(output_raw), SIMPLIFY = FALSE))
  
  # merge loglike and logprior for cold chain into main output
  df_theta <- df_theta %>%
    dplyr::left_join(df_pt, by = c("chain", "rung", "phase", "iteration"))
  
  # if save_hot_draws = TRUE then merge theta values back into pt output
  if (save_hot_draws) {
    df_pt <- df_pt %>%
      dplyr::left_join(dplyr::select(df_theta, -.data$loglikelihood, -.data$logprior), by = c("chain", "rung", "phase", "iteration"))
  }
  
  # drop rungs field from main output
  df_output <- df_theta %>%
    dplyr::filter(.data$rung == max(rungs)) %>%
    dplyr::select(-.data$rung)
  
  # check for bad values in output
  if (!all(is.finite(unlist(df_output[, param_names])))) {
    stop("output contains non-finite values")
  }
  
  # append to output list
  output_processed <- list(output = df_output,
                           pt = df_pt)
  
  ## Diagnostics
  output_processed$diagnostics <- list()
  
  # run-times
  run_time <- data.frame(chain = chain_names,
                         seconds = chain_runtimes)
  output_processed$diagnostics$run_time <- run_time
  
  # Rhat (Gelman-Rubin diagnostic)
  if (chains > 1) {
    rhat_est <- c()
    for (p in seq_along(param_names)) {
      rhat_est[p] <- df_output %>%
        dplyr::filter(phase == "sampling") %>%
        dplyr::select(chain, param_names[p]) %>%
        gelman_rubin(chains = chains, samples = samples)
    }
    rhat_est[skip_param] <- NA
    names(rhat_est) <- param_names
    output_processed$diagnostics$rhat <- rhat_est
  }
  
  # ESS
  ess_est <- df_output %>%
    dplyr::filter(phase == "sampling") %>%
    dplyr::select(param_names) %>%
    apply(2, coda::effectiveSize)
  ess_est[skip_param] <- NA
  output_processed$diagnostics$ess <- ess_est
  
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
    mc_accept <- rbind(cbind(mc_accept, phase = "burnin", value = mc_accept_burnin),
                       cbind(mc_accept, phase = "sampling", value = mc_accept_sampling))
    
  }
  output_processed$diagnostics$mc_accept <- mc_accept
  
  # DIC
  DIC <- df_pt %>%
    dplyr::filter(.data$phase == "sampling" & .data$rung == rungs) %>%
    dplyr::select(.data$loglikelihood) %>%
    dplyr::mutate(deviance = -2*.data$loglikelihood) %>%
    dplyr::summarise(DIC = mean(.data$deviance) + 0.5*var(.data$deviance)) %>%
    dplyr::pull(.data$DIC)
  output_processed$diagnostics$DIC_Gelman <- DIC
  
  ## Parameters
  data_store <- NULL
  if (save_data) {
    data_store <- data
  }
  output_processed$parameters <- list(data = data_store,
                                      df_params = df_params,
                                      loglike = loglike,
                                      logprior = logprior,
                                      burnin = burnin,
                                      samples = samples,
                                      rungs = rungs,
                                      chains = chains,
                                      coupling_on = coupling_on,
                                      alpha = alpha,
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
  
  # Specify pointers to cpp functions
  if (args$args_params$loglike_use_cpp) {
    args$args_functions$loglike <- create_xptr(args$args_functions$loglike)
  }
  if (args$args_params$logprior_use_cpp) {
    args$args_functions$logprior <- create_xptr(args$args_functions$logprior)
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

#------------------------------------------------
# update progress bar
# pb_list = list of progress bar objects
# name = name of this progress bar
# i = new value of bar
# max_i = max value of bar (close when reach this value)
# close = whether to close when reach end
#' @importFrom utils setTxtProgressBar
#' @noRd
update_progress <- function(pb_list, name, i, max_i, close = TRUE) {
  setTxtProgressBar(pb_list[[name]], i)
  if (i == max_i & close) {
    close(pb_list[[name]])
  }
}

# Deal with user input cpp not being defined
globalVariables(c("create_xptr"))
