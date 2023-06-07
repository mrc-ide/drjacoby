#' @title R6 class representing a mcmc object
#'
#' @description 
#' Create a new dj object to run an mcmc
#' @export
dj <- R6::R6Class(
  "dj",
  
  private = list(
    ### private variables ###
    tune_called = FALSE,
    burn_called = FALSE,
    sample_called = FALSE,
    random_number_seed = NULL,
    
    data = NULL,
    df_params = NULL,
    misc = NULL,
    loglikelihood = NULL,
    logprior = NULL,
    
    theta_names = NULL,
    theta_min =  NULL,
    theta_max =  NULL,
    theta_transform_type = NULL,
    n_par = NULL,
    infer_parameter = NULL,
    
    chains = NULL,
    chain_objects = NULL,
    
    blocks = NULL,
    n_unique_blocks = NULL,
    
    iteration_counter = list(
      tune = 0L,
      burn = 0L,
      sample = 0L
    ),
    
    target_acceptance = NULL,
    
    rungs = list(
      tune = 1L,
      burn = 1L,
      sample = 1L
    ),
    beta = list(
      tune = 1,
      burn = 1,
      sample = 1
    ),
    swap = 0L,
    target_rung_acceptance = NULL,
    lambda = NULL,
    phases = c("tune", "burn", "sample"),
    
    output_df = NULL
  ),
  
  public = list(
    ### Public variables ###
    
    #' @field error_debug Diagnostic output in the case of mcmc failure
    error_debug = NULL,
    
    ### Initialisation ###
    
    #' @description
    #' Initialise the dj object.
    #' @param data a named list or data frame or data values.
    #' @param df_params a data.frame of parameters (see \code{?define_params}).
    #' @param loglikelihood the log-likelihood function used in
    #'   the MCMC. Can either be passed in as R functions (not in quotes), or as
    #'   character strings naming compiled cpp11 functions (in quotes).
    #' @param logprior the log-prior function used in
    #'   the MCMC. Can either be passed in as R functions (not in quotes), or as
    #'   character strings naming compiled cpp11 functions (in quotes).
    #' @param chains the number of independent replicates of the MCMC to run.
    #' @param misc optional list object passed to likelihood and prior. This can be
    #'   useful for passing values that are not strictly data, for example passing a
    #'   lookup table to the prior function.
    #' @return A new `dj` object.
    initialize = function(data, df_params, loglikelihood, logprior, chains = 1L, misc = list()){
      
      Sys.setenv(RSTUDIO = "1") # To get progress bar to show on Windows (progress package outstanding issue)
      
      # Input checks
      stopifnot(is.list(data))
      stopifnot(is.list(df_params))
      check_params(df_params)
      stopifnot(is.list(misc))
      stopifnot(is.integer(chains))
      stopifnot(chains >= 1)
      stopifnot(!"block" %in% names(misc))
      
      # Shared elements
      private$chains = chains
      private$data = data
      private$df_params = df_params
      private$misc = misc
      private$loglikelihood = loglikelihood
      private$logprior = logprior
      private$theta_names = unlist(df_params$name)
      private$n_par = length(private$theta_names)
      private$theta_min = unlist(df_params$min)
      private$theta_max = unlist(df_params$max)
      private$theta_transform_type = get_transform_type(private$theta_min,  private$theta_max)
      private$infer_parameter = as.integer(!(private$theta_min == private$theta_max))
      private$blocks = set_blocks(df_params)
      private$n_unique_blocks = length(unique(unlist(private$blocks)))
      private$output_df = vector("list", private$chains)
      
      # Chain-specific elements
      rungs = 1
      theta = initial(
        df_params = private$df_params,
        chains = private$chains,
        rungs = rungs
      )
      proposal_sd = create_proposal_sd(
        chains = private$chains,
        rungs = rungs,
        n_par = private$n_par
      )
      acceptance_counter = create_chain_phase_list(
        chains = private$chains,
        base = matrix(
          0L, 
          nrow = rungs, 
          ncol = private$n_par
        )
      )
      swap_acceptance_counter = create_chain_phase_list(
        chains = private$chains,
        base = rep(0L, rungs - 1)
      )
      duration = create_chain_phase_list(
        chains = private$chains,
        base = 0
      )
      rng_ptr = dust::dust_rng_distributed_pointer(
        seed = NULL,
        n_nodes = private$chains
      )
      
      chain_objects <- list()
      for(i in 1:private$chains){
        chain_objects[[i]] <- list(
          chain = i,
          theta = theta[[i]],
          proposal_sd = proposal_sd[[i]],
          acceptance_counter = acceptance_counter[[i]],
          swap_acceptance_counter = swap_acceptance_counter[[i]],
          duration = duration[[i]],
          rng_ptr = rng_ptr[[i]]
        )
      }
      
      private$chain_objects = chain_objects
    },
    
    ### Print ###
    #' @description
    #' Print mcmc object summary
    print = function(){
      # print summary
      cat("drjacoby object:\n")
      cat(" \U0001F453  Parameters: ", private$n_par, "\n", sep = "")
      cat(" \U0001F453  Chains: ", private$chains, "\n", sep = "")
      cat(" \U0001F453  Tuning rungs: ", private$rungs[["tune"]], "\n", sep = "")
      cat(" \U0001F453  Burn-in and sampling rungs: ", private$rungs[["sample"]], "\n", sep = "")
      cat(" \U0001F453  Tuning iterations: ", private$iteration_counter[["tune"]], "\n", sep = "")
      cat(" \U0001F453  Burn-in iterations: ", private$iteration_counter[["burn"]], "\n", sep = "")
      cat(" \U0001F453  Sampling iterations: ", private$iteration_counter[["sample"]], "\n", sep = "")
      cat(" \U0001F453  Total compute time: ", round(sum(self$timing()$seconds[4,]), 4), " seconds", "\n", sep = "")
      # return invisibly
      invisible(self)
    },
    
    #### Public functions ###
    
    ### Tuning ###
    #' @description
    #' Run tuning. This will dynamically select the optimum rung number and beta
    #' schedule. Proposal standard deviations will also be tuned in this phase but
    #' target acceptance rates and chain convergence are not guaranteed  - run
    #' \code{$burn()} to burn in further.
    #' @param iterations Number of tuning iterations to run
    #' @param swap integer 0 = no swapping, 1 = standard swapping, 2 = super cool new swapping
    #' @param beta Initial beta schedule.
    #' @param max_rungs The maximum number of rungs
    #' @param target_acceptance Target acceptance rate
    #' @param silent print progress (boolean)
    tune = function(target_rung_acceptance = 0.5, n_rungs = 50, iterations = 1000, initial_beta = NULL, swap = 2L,  target_acceptance = 0.44,  silent = FALSE){
      if(private$chains > 1){
        stop("To use parallel tempering please set the number of chains = 1")
      }
      if(private$burn_called | private$sample_called){
        stop("Cannot call tune after burn or sample have been called")
      }
      phase <- "tune"
      private$tune_called <- TRUE
      burnin <- TRUE
      
      # Tuning beta initialised with simple power law, unless overridden by user with initial_beta
      private$beta[[phase]] <- seq(1, 0, length.out = n_rungs)^2
      if(!is.null(initial_beta)){
        private$beta[[phase]] <- initial_beta
      }
      private$rungs[[phase]] <- length(private$beta[[phase]])
      
      private$swap <- swap
      private$target_acceptance <- target_acceptance
      private$target_rung_acceptance <- target_rung_acceptance
      
      # Adjust objects dependent on number of rungs
      private$chain_objects[[1]]$theta = initial(private$df_params, chains = private$chains, rungs = private$rungs[[phase]])[[1]]
      private$chain_objects[[1]]$proposal_sd = create_proposal_sd(private$chains, private$rungs[[phase]], private$n_par)[[1]]
      private$chain_objects[[1]]$acceptance_counter[[phase]] = matrix(0L, nrow = private$rungs[[phase]], ncol = private$n_par)
      private$chain_objects[[1]]$swap_acceptance_counter[[phase]] <- rep(0L, private$rungs[[phase]] - 1)
      
      # For tuning, chains always = 1, so never in parallel
      apply_func <- noclusterApply
      
      # Run chains
      chain_output <- apply_func(
        cl = NULL, 
        x = private$chain_objects,
        fun = run_mcmc,
        phase = phase,
        burnin = burnin,
        iterations = iterations,
        silent = silent,
        theta_names = private$theta_names,
        theta_transform_type = private$theta_transform_type,
        theta_min = private$theta_min,
        theta_max = private$theta_max,
        infer_parameter = private$infer_parameter,
        data = private$data,
        loglikelihood = private$loglikelihood,
        logprior = private$logprior,
        misc = private$misc,
        target_acceptance = private$target_acceptance,
        swap = private$swap,
        beta = private$beta[[phase]],
        blocks = private$blocks,
        n_unique_blocks = private$n_unique_blocks,
        iteration_counter = private$iteration_counter[[phase]]
      )
      
      # Error checking
      lapply(chain_output, function(x){
        if("error" %in% names(x)){
          self$error_debug = x
          stop("Error in mcmc, check $error_debug")
        }
      })
      
      # Update elemtents
      private$iteration_counter[[phase]] <- private$iteration_counter[[phase]] + iterations
      private$chain_objects[[1]]$duration[[phase]] = private$chain_objects[[1]]$duration[[phase]] + chain_output[[1]]$dur
      private$chain_objects[[1]]$acceptance_counter[[phase]] = chain_output[[1]]$acceptance
      private$chain_objects[[1]]$swap_acceptance_counter[[phase]] = chain_output[[1]]$swap_acceptance
      
      ### Define final beta ###
      # Check all rung pairs have achieved some swaps
      tune_rejection_rate <- 1 - self$mc_acceptance_rate("tune")
      if(any(tune_rejection_rate > 0.99)){
        stop("Tuning needs more rungs to achieve an accurate estimate of communication barrier")
      }
      private$lambda <-  sum(tune_rejection_rate)
      rungs <- ceiling(private$lambda / (1 - private$target_rung_acceptance))
      private$rungs[["burn"]] <- rungs
      private$rungs[["sample"]] <- rungs
      new_beta <- propose_new_beta(
        n = rungs,
        beta_mid = beta_mid(private$beta[[phase]]),
        rejection_rate = tune_rejection_rate,
        lambda = private$lambda
      )
      private$beta[["burn"]] <- new_beta
      private$beta[["sample"]] <- new_beta
      
      # Adjust objects dependent on number of rungs
      private$chain_objects[[1]]$theta = initial(private$df_params, chains = private$chains, rungs = rungs)[[1]]
      private$chain_objects[[1]]$proposal_sd = create_proposal_sd(private$chains, rungs, private$n_par)[[1]]
      private$chain_objects[[1]]$swap_acceptance_counter[["burn"]] = rep(0L, rungs - 1)
      private$chain_objects[[1]]$swap_acceptance_counter[["sample"]] = rep(0L, rungs - 1)
      private$chain_objects[[1]]$acceptance_counter[["burn"]] = matrix(0L, nrow = rungs, ncol = private$n_par)
      private$chain_objects[[1]]$acceptance_counter[["sample"]] = matrix(0L, nrow = rungs, ncol = private$n_par)
      
      # Best-guess to inform starting proposal_sd by matching new beta to closest tuning beta
      ic <- index_closest(new_beta, private$beta[[phase]])
      private$chain_objects[[1]]$proposal_sd = chain_output[[1]]$proposal_sd[ic,]
      
    },
    
    ### Burn in ###
    #' @description
    #' Run burn in. Runs the burn in phase of the MCMC to allow where 
    #' proposal standard deviations to be tuned towards a target acceptance rate
    #' and chains to converge.
    #' @param iterations Number of burn-in iterations to run
    #' @param target_acceptance Target acceptance rate
    #' @param silent print progress (boolean)
    #' @param cl parallel cluster object 
    burn = function(iterations, target_acceptance = 0.44, silent = FALSE, cl = NULL){
      stopifnot(is.integer(iterations))
      if(private$sample_called){
        stop("Cannot call burn after sample has been called")
      }
      
      # Set burn parameters
      private$burn_called <- TRUE
      private$target_acceptance <- target_acceptance
      burnin <- TRUE
      phase <- "burn"
      
      # Define the apply function
      is_parallel <- !is.null(cl)
      apply_func <- noclusterApply
      if(is_parallel){
        apply_func <- parallel::clusterApply
      }
      # Run chains
      chain_output <- apply_func(        
        cl = cl, 
        x = private$chain_objects,
        fun = run_mcmc,
        phase = phase,
        burnin = burnin,
        iterations = iterations,
        silent = silent,
        theta_names = private$theta_names,
        theta_transform_type = private$theta_transform_type,
        theta_min = private$theta_min,
        theta_max = private$theta_max,
        infer_parameter = private$infer_parameter,
        data = private$data,
        loglikelihood = private$loglikelihood,
        logprior = private$logprior,
        misc = private$misc,
        target_acceptance = private$target_acceptance,
        swap = private$swap,
        beta = private$beta[[phase]],
        blocks = private$blocks,
        n_unique_blocks = private$n_unique_blocks,
        iteration_counter = private$iteration_counter[[phase]]
      )
      
      # Error checking
      lapply(chain_output, function(x){
        if("error" %in% names(x)){
          self$error_debug = x
          stop("Error in mcmc, check $error_debug")
        }
      })
      
      private$iteration_counter[[phase]] <- private$iteration_counter[[phase]] + iterations
      
      # Update R6 object with mcmc outputs
      for(i in 1:private$chains){
        # Update user output
        private$output_df[[i]] <- append_output(
          current = private$output_df[[i]],
          new = chain_output[[i]]$output,
          phase = phase,
          theta_names = private$theta_names,
          chain = i
        )
        # Update internal states
        private$chain_objects[[i]]$duration[[phase]] = private$chain_objects[[i]]$duration[[phase]] + chain_output[[i]]$dur
        private$chain_objects[[i]]$proposal_sd = chain_output[[i]]$proposal_sd
        private$chain_objects[[i]]$acceptance_counter[[phase]] = chain_output[[i]]$acceptance
        private$chain_objects[[i]]$theta = chain_output[[i]]$theta
        if(private$rungs[[phase]] > 1){
          private$chain_objects[[i]]$swap_acceptance_counter[[phase]] = chain_output[[i]]$swap_acceptance
        }
        if(is_parallel){
          private$chain_objects[[i]]$rng_ptr <- chain_output[[i]]$rng_ptr
        }
      }
    },
    
    ### Sampling ###
    #' @description
    #' Run sampling in. Runs the sampling phase of the MCMC.
    #' @param iterations Number of sampling iterations to run
    #' @param silent print progress (boolean)
    #' @param cl parallel cluster object 
    sample = function(iterations, silent = FALSE, cl = NULL){
      stopifnot(is.integer(iterations))
      
      # Set sample parameters
      private$sample_called <- TRUE
      burnin = FALSE
      phase <- "sample"
      
      # Define the apply function
      is_parallel <- !is.null(cl)
      apply_func <- noclusterApply
      if(is_parallel){
        apply_func <- parallel::clusterApply
      }
      # Run chains
      chain_output <- apply_func(
        cl = cl, 
        x = private$chain_objects,
        fun = run_mcmc,
        phase = phase,
        burnin = burnin,
        iterations = iterations,
        silent = silent,
        theta_names = private$theta_names,
        theta_transform_type = private$theta_transform_type,
        theta_min = private$theta_min,
        theta_max = private$theta_max,
        infer_parameter = private$infer_parameter,
        data = private$data,
        loglikelihood = private$loglikelihood,
        logprior = private$logprior,
        misc = private$misc,
        target_acceptance = private$target_acceptance,
        swap = private$swap,
        beta = private$beta[[phase]],
        blocks = private$blocks,
        n_unique_blocks = private$n_unique_blocks,
        iteration_counter = private$iteration_counter[[phase]]
      )
      
      # Error checking
      lapply(chain_output, function(x){
        if("error" %in% names(x)){
          self$error_debug = x
          stop("Error in mcmc, check $error_debug")
        }
      })
      
      private$iteration_counter[[phase]] <- private$iteration_counter[[phase]] + iterations
      
      # Update R6 object with mcmc outputs
      for(i in 1:private$chains){
        # Update user output
        private$output_df[[i]] <- append_output(
          current = private$output_df[[i]],
          new = chain_output[[i]]$output,
          phase = phase,
          theta_names = private$theta_names,
          chain = i
        )
        # Update internal states
        private$chain_objects[[i]]$duration[[phase]] = private$chain_objects[[i]]$duration[[phase]] + chain_output[[i]]$dur
        private$chain_objects[[i]]$theta = chain_output[[i]]$theta
        private$chain_objects[[i]]$acceptance_counter[[phase]] = chain_output[[i]]$acceptance
        if(private$rungs[[phase]] > 1){
          private$chain_objects[[i]]$swap_acceptance_counter[[phase]] = chain_output[[i]]$swap_acceptance
        }
        if (is_parallel) {
          private$chain_objects[[i]]$rng_ptr <- chain_output[[i]]$rng_ptr
        }
      }
    },
    
    ### Output ###
    #' @description
    #' Get mcmc output data.frame
    #' @param chain option chain(s) selection
    #' @param phase optional phase selection, can be a vector chosen from "tune", "burn","sample"
    output = function(chain = NULL, phase = NULL){
      data <- list_r_bind(private$output_df)
      if(!is.null(phase)){
        data <- data[data$phase %in% phase, ]
      }
      if(!is.null(chain)){
        data <- data[data$chain %in% chain, ]
      }
      return(data)
    },
    
    ### Diagnostics ###
    #' @description
    #' Get acceptance rates
    acceptance_rate = function(phase = "sample", chain = NULL, rung = 1){
      acceptance_counter <- chain_element(private$chain_objects, "acceptance_counter")
      
      phases <- private$phases
      if(!private$tune_called){
        phases <- phases[!phases == "tune"]
      }
      if(!is.null(phase)){
        phases <- phases[phases %in% phase]
      }
      
      chains <- 1:private$chains
      if(!is.null(chain)){
        chains <- chains[chains %in% chain]
      }
      
      estimate_acceptance_rate(
        acceptance_counter = acceptance_counter,
        iteration_counter = private$iteration_counter,
        chains = chains,
        phases = phases,
        theta_names = private$theta_names,
        rungs = rung
      )
      
    },
    
    #' @description
    #' Get acceptance rates
    mc_acceptance_rate = function(phase = "sample"){
      stopifnot(phase %in% private$phases)
      swap_acceptance_counter <- private$chain_objects[[1]]$swap_acceptance_counter[[phase]]
      iteration_counter <- private$iteration_counter[[phase]]
      if(private$swap == 2){
        iteration_counter = iteration_counter / 2
      }
      swap_acceptance_counter / iteration_counter
    },
    
    #' @description
    #' Get DIC estimate
    dic = function(){
      output_df <- list_r_bind(private$output_df)
      estimate_dic(output_df)
    },
    
    #' @description
    #' Get effective sample size estimates
    ess = function(){
      output_df <- list_r_bind(private$output_df)
      estimate_ess(output_df, private$theta_names)
    },
    
    #' @description
    #' Get rhat estimates
    rhat = function(){
      output_df <- list_r_bind(private$output_df)
      iteration_counter <- private$iteration_counter
      iteration_counter <- list_c_bind(iteration_counter)
      estimate_rhat(output_df, private$theta_names, private$chains, iteration_counter)
    },
    
    #' @description
    #' Get run-time data
    timing = function(){
      duration <- chain_element(private$chain_objects, "duration")
      iteration_counter <- private$iteration_counter
      seconds <- list_c_bind(duration)
      iterations <- list_c_bind(iteration_counter)
      estimate_timing(
        seconds = seconds,
        iterations = iterations,
        phases = private$phases,
        chains = private$chains)
    },
    
    ### Plots ###
    
    #' @description Produce a parameter plot of the named parameter,
    #'   including the raw trace, the posterior histogram and an autocorrelation
    #'   plot. The combined plot or a list of plot elements can be returned.
    #'
    #' @param par The name of the parameter
    #' @param lag Maximum lag. Must be an integer between 1 and 500.
    #' @param downsample boolean. Whether to downsample chain to make plotting more
    #'   efficient.
    #' @param phase Optional selection of phases, can be from: tune, burn and sample
    #' @param chain Optional selection of chains
    #' @param return_elements boolean. If \code{TRUE} a list of plotting objects 
    #'   are returned without displaying.
    plot_par = function(par, lag = 20, downsample = TRUE, phase = NULL, chain = NULL, return_elements = FALSE){
      create_par_plot(
        par = par,
        output_df = private$output_df,
        lag = lag,
        downsample = downsample,
        phase = phase,
        chain = chain,
        return_elements = return_elements
      )
    },
    
    #' @description Plots the correlation between two parameters
    #'
    #' @param parx Name of parameter 1
    #' @param pary Name of parameter 2
    #' @param downsample boolean. Whether to downsample chain to make plotting more
    #'   efficient.
    #' @param phase Optional selection of phases, can be from: tune, burn and sample
    #' @param chain Optional selection of chains
    plot_cor = function(parx, pary, downsample = TRUE, phase = "sample", chain = NULL){
      create_cor_plot(
        parx = parx,
        pary = pary,
        output_df = private$output_df,
        downsample = downsample,
        phase = phase,
        chain = chain
      )
    },
    
    #' @description Plots posterior 95\% credible intervals over specified set of
    #'   parameters (defauls to all parameters).
    #'
    #' @param pars Vector of parameter names
    #' @param phase Optional selection of phases, can be from: tune, burn and sample
    #' @param chain Optional selection of chains
    #' @param param_names Optional vector of prameter names for plotting labels
    plot_credible = function(pars = NULL, phase = "sample", chain = NULL, param_names = NULL){
      if(is.null(pars)){
        pars = private$theta_names
      }
      if(is.null(param_names)){
        param_names <- pars
      }
      create_credible_plot(
        output_df = private$output_df,
        pars = pars,
        chain = chain,
        phase = phase,
        param_names = param_names
      )
    },
    
    #' @description Produces a matrix showing the correlation between all parameters
    #'   from posterior draws.
    #'   
    #' @param pars Vector of parameter names
    #' @param phase Optional selection of phases, can be from: tune, burn and sample
    #' @param chain Optional selection of chains
    #' @param param_names Optional vector of prameter names for plotting labels
    plot_cor_mat = function(pars, phase = "sample", chain = NULL, param_names = NULL){
      if(is.null(pars)){
        pars = private$theta_names
      }
      if(is.null(param_names)){
        param_names <- pars
      }
      create_cor_mat_plot(
        output_df = private$output_df,
        pars = pars,
        chain = chain,
        phase = phase,
        param_names = param_names)
    },
    
    #' @description Plots the acceptance rate between rungs
    #'
    #' @param phase Optional selection of phases, can be from: tune, burn or sample
    #' @param x_axis_type how to format the x-axis. 1 = integer rungs, 2 = values of
    #'   the thermodynamic power.
    plot_mc_acceptance_rate = function(phase = "sample", x_axis_type = 1){
      if(private$rungs == 1){
        stop("Not available for a single rung")
      }
      ar <- self$mc_acceptance_rate()
      ar <- ar[which(phase == private$phases),]
      create_mc_acceptance_plot(
        rungs = private$rungs,
        beta = private$beta,
        ar = ar,
        x_axis_type = x_axis_type
      )
    },
    
    plot_tuning_rejection_rate = function(){
      create_rejection_rate_plot(private$tune_beta, private$tune_beta_mid, private$tune_rejection_rate)
    },
    
    plot_local_communication_barrier = function(){
      create_local_communication_barrier_plot(private$tune_beta, private$tune_beta_mid, private$tune_rejection_rate)
    }
  )
)

#' @import R6
NULL
#' @import dust
NULL
#' @import parallel
NULL
#' @import cpp11
NULL
#' @import decor
NULL
