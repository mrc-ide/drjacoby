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
#' @param data a named list or data frame or data values.
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
run_mcmc <- function(
    data,
    df_params,
    loglike,
    logprior,
    burnin,
    samples,
    misc = list(),
    rungs = 1L,
    chains = 1L,
    beta_manual = NULL,
    alpha = 1,
    target_acceptance = 0.44,
    cluster = NULL,
    coupling_on = TRUE,
    silent = TRUE
){

  ### Inputs ###################################################################
  # Input checks
  stopifnot(is.list(data))
  stopifnot(is.list(df_params))
  stopifnot(is.integer(burnin))
  stopifnot(is.integer(samples))
  stopifnot(is.integer(rungs))
  stopifnot(is.integer(chains))
  
  use_init <- ("init" %in% names(df_params))
  if(use_init){
    check_init(df_params, chains)
  } else {
    df_params$init <- get_init(df_params, chains)
  }
  if(!"block" %in% names(df_params)){
    df_params$block <- 1
  }
  
  if(is.null(beta_manual)){
    beta_manual <- seq(1, 0, length.out = rungs)
  }
  beta_raised <- beta_manual^alpha
  
  infer_parameter <- as.integer(!df_params$max == df_params$min)
  
  # List inputs - to distribute if running in parallel
  input <- lapply(1:chains, function(x){
    input <- list(
      chain = x,
      theta_init = sapply(df_params$init, '[', x),
      theta_names = unlist(df_params$name),
      theta_min =  unlist(df_params$min),
      theta_max =  unlist(df_params$max),
      theta_transform_type = get_transform_type(unlist(df_params$min), unlist(df_params$max)),
      blocks_list = lapply(df_params$block, as.integer),
      n_unique_blocks = length(unique(unlist(lapply(df_params$block, as.integer)))),
      data = data,
      burnin = burnin,
      samples = samples,
      loglike = loglike,
      logprior = logprior,
      target_acceptance = target_acceptance,
      misc = misc,
      rungs = rungs,
      beta_init = beta_raised,
      swap = coupling_on,
      chains = chains,
      infer_parameter = infer_parameter,
      silent = silent
    )
    return(input)
  })
  # Currently need to do this to get progress bar to display on windows 
  # See https://github.com/r-lib/progress/issues/56
  if(!silent){
    Sys.setenv(RSTUDIO = "1")
  }
  ##############################################################################
  
  ### Run MCMC #################################################################
  if(is.null(cluster)){
    mcmc_runs <- lapply(
      X = input,
      FUN = run_internal
    )
  } else {
    parallel::clusterEvalQ(
      cl = cluster,
      expr = library(dj11)
    )
    mcmc_runs <- parallel::clusterApplyLB(
      cl = cluster,
      x = input,
      fun = run_internal
    )
  }
  ##############################################################################
  
  ### MCMC output ##############################################################
  if("error" %in% sapply(mcmc_runs, names)){
    warning("One or more chains produced errors, returning error output")
    return(mcmc_runs)
  }
  out <- list()
  out$output <- dplyr::bind_rows(sapply(mcmc_runs, '[', 'output'))
  ##############################################################################
  
  ### Diagnostics ##############################################################
  # DIC
  out$diagnostics$DIC_Gelman <- dic(out$output)
  
  # MC
  out$diagnostics$mc_accept <- NULL
  out$diagnostics$rung_index <- NULL
  if(rungs > 1){
    # MC acceptance
    out$diagnostics$mc_accept <- mc_acceptance(
      mcmc_runs = mcmc_runs,
      chains = chains,
      rungs = rungs,
      burnin = burnin,
      samples = samples
    )
    # Rung index
    out$diagnostics$rung_index <- dplyr::bind_rows(
      sapply(mcmc_runs, '[', 'rung_index')
    )
    # Thermodynamic power
    out$diagnostics$rung_details <- data.frame(rung = 1:rungs,
                                               thermodynamic_power = beta_raised)
  }
  # ESS
  out$diagnostics$ess <- ess(
    output = out$output,
    parameter_names = input[[1]]$theta_names
  )
  
  # Rhat (Gelman-Rubin diagnostic)
  out$diagnostics$rhat <- "rhat not relevant for a single chain"
  if(chains > 1){
    out$diagnostics$rhat <- rhat(
      output = out$output,
      parameter_names = input[[1]]$theta_names,
      n_chains = chains,
      samples= samples
    )
  }
  ##############################################################################
  
  ### Parameters ###############################################################
  out$parameters <- input[c("data",
                            "df_params",
                            "loglike",
                            "logprior",
                            "burnin",
                            "samples",
                            "n_rungs",
                            "chains",
                            "swap")]
  ##############################################################################
  
  out <- out[c("output", "diagnostics", "parameters")]
  class(out) <- "drjacoby_output"
  return(out)
}

run_internal <- function(input){
  # Convert cpp11 loglikelihood
  if(is.character(input$loglike)){
    input$loglike <- get(input$loglike)
  }
  # Convert cpp11 logprior
  if(is.character(input$logprior)){
    input$logprior <- get(input$logprior)
  }
  # Run mcmc
  mcmc_out <- mcmc(input$chain, input$theta_init, input$theta_names, input$theta_transform_type,  input$theta_min,  input$theta_max,
                   input$blocks_list, input$n_unique_blocks, input$data, input$burnin, input$samples, input$loglike, input$logprior,
                   input$target_acceptance, input$misc, input$rungs, input$beta_init, input$swap, input$infer_parameter, input$silent)
  
  if("error" %in% names(mcmc_out)){
    return(mcmc_out)
  }
  
  # Add Chain, burnin, sampling columns
  mcmc_out$output <- cbind(
    data.frame(chain = input$chain, phase = rep(c("burnin", "sampling"), c(input$burnin, input$samples))),
    mcmc_out$output
  )
  names(mcmc_out$output) <- c("chain", "phase", "iteration", input$theta_names, "logprior", "loglikelihood")
  
  return(mcmc_out)
}

get_transform_type <- function(theta_min, theta_max){
  as.integer(2 * is.finite(theta_min) + is.finite(theta_max))
}
