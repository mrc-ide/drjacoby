run_dj11 <- function(
    data,
    df_params,
    loglike,
    logprior,
    burnin,
    samples,
    misc = list(),
    rungs = 1L,
    chains = 1,
    beta_manual = NULL,
    target_acceptance = 0.44,
    cluster = NULL,
    coupling_on = TRUE,
    progress = TRUE
){
  # TODO: progress bar
  # TODO: beta manual
  # TODO: skip params
  # TODO: Inf likelihood
  
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
    df_params$init <- get_init(df_params)
  }
  
  # List inputs - to distribute if running in parallel
  input <- lapply(1:chains, function(x){
    input <- list()
    input$chain = x
    input$theta_init <- sapply(df_params$init, '[', x)
    input$theta_names <- unlist(df_params$name)
    input$theta_min <-  unlist(df_params$min)
    input$theta_max <-  unlist(df_params$max)
    input$theta_transform_type <- get_transform_type(input$theta_min, input$theta_max)
    input$blocks_list <- lapply(df_params$block, as.integer)
    input$n_unique_blocks <- length(unique(unlist(input$blocks_list)))
    input$data <- data
    input$burnin <- burnin
    input$samples <- samples
    input$loglike <- loglike
    input$logprior <- logprior
    input$target_acceptance <- target_acceptance
    input$misc <- misc
    input$n_rungs <- rungs
    input$beta_init <- seq(1, 0, length.out = input$n_rungs)
    input$swap <- coupling_on
    input$chains <- chains
    return(input)
  })
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
    out$diagnostics$mc_accept <- dplyr::bind_rows(
      sapply(mcmc_runs, '[', 'swap_acceptance')
    )
    # Rung index
    out$diagnostics$rung_index <- dplyr::bind_rows(
      sapply(mcmc_runs, '[', 'rung_index')
    )
  }
  
  # ESS
  out$diagnostics$ess <- ess(
    output = out$output,
    parameter_names = input$theta_names
  )
  
  # Rhat (Gelman-Rubin diagnostic)
  out$diagnostics$rhat <- NA
  if(chains > 1){
    out$diagnostics$rhat <- rhat(
      output = out$output,
      parameter_names = input$theta_names,
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
  mcmc_out <- mcmc(input$theta_init, input$theta_names, input$theta_transform_type,  input$theta_min,  input$theta_max,
                   input$blocks_list, input$n_unique_blocks, input$data, input$burnin, input$samples, input$loglike, input$logprior,
                   input$target_acceptance, input$misc, input$n_rungs, input$beta_init, input$swap)
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
