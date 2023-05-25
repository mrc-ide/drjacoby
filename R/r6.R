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
    
    target_acceptance = NULL,
    rungs = 1,
    beta = 1,
    swap = 0L,
    
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
      private$output_df <- vector("list", private$chains)
      
      # Chain-specific elements
      theta = initial(df_params, chains = private$chains, rungs = private$rungs)
      
      proposal_sd = lapply(1:private$chains, function(x){
        matrix(0.1, nrow = private$rungs, ncol = private$n_par)
      })
        
      acceptance_counter = lapply(1:private$chains, function(x){
        list(
          Tune = matrix(0L, nrow = private$rungs, ncol = private$n_par),
          Burn = matrix(0L, nrow = private$rungs, ncol = private$n_par),
          Sample = matrix(0L, nrow = private$rungs, ncol = private$n_par)
        )
      })
      
      iteration_counter = lapply(1:private$chains, function(x){
        matrix(0, nrow = 3, ncol = 1, dimnames = list(
          c("Tune", "Burn", "Sample"),
          paste0("Chain_", x)
        )
        )
      })
      
      swap_acceptance_counter = lapply(1:private$chains, function(x){
        rep(0L, private$rungs)
      })
      
      duration = lapply(1:private$chains, function(x){
        matrix(0, nrow = 3, ncol = 1, dimnames = list(
          c("Tune", "Burn", "Sample"),
          paste0("Chain_", x)
        )
        )
      })
      
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
          iteration_counter = iteration_counter[[i]],
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
      cat("  Chains: ", private$chains, "\n", sep = "")
      cat("  Rungs: ", private$rungs, "\n", sep = "")
      cat("  Tuning iterations: ", private$iteration_counter[1, 1], "\n", sep = "")
      cat("  Burn-in iterations: ", private$iteration_counter[2, 1], "\n", sep = "")
      cat("  Sampling iterations: ", private$iteration_counter[3, 1], "\n", sep = "")
      cat("  Parameters: ", private$n_par, "\n", sep = "")
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
    #' @param beta Initial beta schedule
    #' @param max_rungs The maximum number of rungs
    #' @param target_acceptance Target acceptance rate
    #' @param silent print progress (boolean)
    tune = function(iterations, swap = 1L, beta = seq(1, 0, -0.1), max_rungs = 100, target_acceptance = 0.44, silent = FALSE){
      private$tune_called <- TRUE
      phase <- 1
      # Infer rung number
      private$rungs = length(beta)
      private$beta = beta
      private$swap = swap
      # Update default initial values for new number of rungs:
      theta = initial(private$df_params, chains = private$chains, rungs = private$rungs)
      
      proposal_sd = lapply(1:private$chains, function(x){
        matrix(0.1, nrow = private$rungs, ncol = private$n_par)
      })
      
      acceptance_counter = lapply(1:private$chains, function(x){
        list(
          Tune = matrix(0L, nrow = private$rungs, ncol = private$n_par),
          Burn = matrix(0L, nrow = private$rungs, ncol = private$n_par),
          Sample = matrix(0L, nrow = private$rungs, ncol = private$n_par)
        )
      })
      
      swap_acceptance_counter = lapply(1:private$chains, function(x){
        rep(0L, private$rungs)
      })
      
      for(i in 1:private$chains){
        private$chain_objects[[i]]$theta = theta[[i]]
        private$chain_objects[[i]]$proposal_sd = proposal_sd[[i]]
        private$chain_objects[[i]]$acceptance_counter = acceptance_counter[[i]]
        private$chain_objects[[i]]$swap_acceptance_counter = swap_acceptance_counter[[i]]
      }
      # As well as updating internal elements as in burnin,
      # the following will all need to be updated if this function is called:
      ## private$rungs
      ## private$beta
      ## private$swap = FALSE
      ## private$swap_acceptance_counter
      ## private$proposal_sd - to increase ncol to match rungs
      ## private$acceptance_counter  - to increase ncol to match rungs
      ## private$swap_acceptance_counter - to increase length to match rungs
      ## private$iteration_counter - keep track of total iterations for tuning phase
      
      # Will also need to add in below new versions of:
      ## plot_mc_acceptance() 
      ## mc_acceptance() 
      # With internal functions
      ## create_mc_acceptance_plot()
      ## estimates_mc_acceptance()
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
      
      # Define the apply function
      is_parallel <- !is.null(cl)
      apply_func <- noclusterApply
      if(is_parallel){
        apply_func <- parallel::clusterApply
      }
      
      # Set burn parameters
      private$burn_called <- TRUE
      private$target_acceptance <- target_acceptance
      burnin <- TRUE
      phase <- 2
      
      # Run chains
      chain_output <- apply_func(cl, 
                                 private$chain_objects, run_mcmc,
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
                                 beta = private$beta,
                                 blocks = private$blocks,
                                 n_unique_blocks = private$n_unique_blocks
      )

      # Update R6 object with mcmc outputs
      for(i in 1:private$chains){
        # Check for error return
        if("error" %in% names(chain_output[[i]])){
          self$error_debug = chain_output[[i]]
          stop("Error in mcmc, check $error_debug")
        }
        # Update user output
        private$output_df[[i]] <- append_output(
          current = private$output_df[[i]],
          new = chain_output[[i]]$output,
          phase = "burn",
          theta_names = private$theta_names,
          chain = i
        )
        # Update internal states
        private$chain_objects[[i]]$duration[phase] = private$chain_objects[[i]]$duration[phase] + chain_output[[i]]$dur
        private$chain_objects[[i]]$iteration_counter[phase] = private$chain_objects[[i]]$iteration_counter[phase] + iterations
        private$chain_objects[[i]]$proposal_sd = chain_output[[i]]$proposal_sd
        private$chain_objects[[i]]$acceptance_counter[[phase]] = chain_output[[i]]$acceptance
        private$chain_objects[[i]]$theta = chain_output[[i]]$theta
        private$chain_objects[[i]]$swap_acceptance_counter = chain_output[[i]]$swap_acceptance
        if (is_parallel) {
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
      
      # Define the apply function
      is_parallel <- !is.null(cl)
      apply_func <- noclusterApply
      if(is_parallel){
        apply_func <- parallel::clusterApply
      }
      
      # Set sample parameters
      private$sample_called <- TRUE
      burnin = FALSE
      phase <- 3
      # Run chains
      chain_output <- apply_func(cl, 
                                 private$chain_objects, run_mcmc,
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
                                 beta = private$beta,
                                 blocks = private$blocks,
                                 n_unique_blocks = private$n_unique_blocks
      )
      
      # Update R6 object with mcmc outputs
      for(i in 1:private$chains){
        # Check for error return
        if("error" %in% names(chain_output[[i]])){
          self$error_debug = chain_output[[i]]
          stop("Error in mcmc, check $error_debug")
        }
        # Update user output
        private$output_df[[i]] <- append_output(
          current = private$output_df[[i]],
          new = chain_output[[i]]$output,
          phase = "sample",
          theta_names = private$theta_names,
          chain = i
        )
        # Update internal states
        private$chain_objects[[i]]$duration[phase] = private$chain_objects[[i]]$duration[phase] + chain_output[[i]]$dur
        private$chain_objects[[i]]$iteration_counter[phase] = private$chain_objects[[i]]$iteration_counter[phase] + iterations
        private$chain_objects[[i]]$theta = chain_output[[i]]$theta
        private$chain_objects[[i]]$acceptance_counter[[phase]] = chain_output[[i]]$acceptance
        private$chain_objects[[i]]$swap_acceptance_counter = chain_output[[i]]$swap_acceptance
        if (is_parallel) {
          private$chain_objects[[i]]$rng_ptr <- chain_output[[i]]$rng_ptr
        }
      }
    },
    
    ### Output ###
    #' @description
    #' Get mcmc output data.frame
    #' @param chain option chain(s) selection
    #' @param phase optional phase selection
    output = function(chain = NULL, phase = NULL){
      data <- list_r_bind(private$output_df)
      if(!is.null(phase)){
        return(data[data$phase %in% phase, ])
      }
      if(!is.null(chain)){
        return(data[data$chain %in% chain, ])
      }
      return(data)
    },
    
    ### Diagnostics ###
    #' @description
    #' Get acceptance rates
    acceptance_rate = function(){
      acceptance_counter <- chain_element(private$chain_objects, "acceptance_counter")
      iteration_counter <- chain_element(private$chain_objects, "iteration_counter")
      iteration_counter <- list_c_bind(iteration_counter)
      estimate_acceptance_rate(
        acceptance_counter,
        iteration_counter,
        private$chains,
        private$theta_names
      )
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
      iteration_counter <- chain_element(private$chain_objects, "iteration_counter")
      iteration_counter <- list_c_bind(iteration_counter)
      estimate_rhat(output_df, private$theta_names, private$chains, iteration_counter)
    },
    
    #' @description
    #' Get run-time data
    timing = function(){
      duration <- chain_element(private$chain_objects, "duration")
      iteration_counter <- chain_element(private$chain_objects, "iteration_counter")
      seconds <- list_c_bind(duration)
      iterations <- list_c_bind(iteration_counter)
      estimate_timing(seconds, iterations)
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
    #' @param phase Optional selection of phases, can be from: all, tune, burn and sample
    #' @param chain Optional selection of chains
    #' @param return_elements boolean. If \code{TRUE} a list of plotting objects 
    #'   are returned without displaying.
    plot_par = function(par, lag = 20, downsample = TRUE, phase = "sample", chain = NULL, return_elements = FALSE){
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
    #' @param phase Optional selection of phases, can be from: all, tune, burn and sample
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
    #' @param phase Optional selection of phases, can be from: all, tune, burn and sample
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
    #' @param phase Optional selection of phases, can be from: all, tune, burn and sample
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
