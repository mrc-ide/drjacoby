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
    
    theta = NULL,
    theta_names = NULL,
    theta_min =  NULL,
    theta_max =  NULL,
    theta_transform_type = NULL,
    infer_parameter = NULL,
    
    chains = NULL,
    
    blocks = NULL,
    n_unique_blocks = NULL,
    
    proposal_sd = NULL,
    target_acceptance = NULL,
    acceptance_counter = NULL,
    
    rungs = 1,
    beta = 1,
    swap = 0L,
    swap_acceptance_counter = rep(0L, 1),
    
    iteration_counter = NULL,
    duration = NULL,
    
    output_df = NULL,
    rng_ptr = NULL
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
    #' @param seed Seed for reproducible random number generation. Must be an integer
    #' @return A new `dj` object.
    initialize = function(data, df_params, loglikelihood, logprior, chains = 1L, misc = list(), seed = NULL){
      
      Sys.setenv(RSTUDIO = "1") # To get progress bar to show on Windows (progress package outstanding issue)
      stopifnot(is.list(data))
      stopifnot(is.list(df_params))
      check_params(df_params)
      stopifnot(is.list(misc))
      stopifnot(is.integer(chains))
      stopifnot(chains >= 1)
      
      private$data = data
      private$df_params = df_params
      private$misc = misc
      private$loglikelihood = loglikelihood
      private$logprior = logprior
      
      private$theta = set_init(df_params, chains)
      private$theta_names = unlist(df_params$name)
      private$theta_min = unlist(df_params$min)
      private$theta_max = unlist(df_params$max)
      private$theta_transform_type = get_transform_type(private$theta_min,  private$theta_max)
      private$proposal_sd = matrix(0.1, nrow = length(private$theta_names), ncol = private$rungs)
      private$infer_parameter = as.integer(!(private$theta_min == private$theta_max))
      
      private$chains = chains
      private$iteration_counter = matrix(
        0,
        nrow = 3,
        ncol = private$chains,
        dimnames = list(
          c("Tune", "Burn", "Sample"),
          paste0("Chain_", 1:private$chains)
        ))
      
      private$duration = matrix(
        0,
        nrow = 3,
        ncol = private$chains,
        dimnames = list(
          c("Tune", "Burn", "Sample"),
          paste0("Chain_", 1:private$chains)
        ))
      
      private$blocks = set_blocks(df_params)
      private$n_unique_blocks = length(unique(unlist(private$blocks)))
      
      private$acceptance_counter = lapply(1:private$chains, function(x){
        matrix(0L, nrow = length(private$theta_names), ncol = private$rungs)
      })
      private$output_df <- vector("list", private$chains)
      private$rng_ptr = dust::dust_rng_pointer$new(seed = seed, n_streams = private$chains)
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
      cat("  Parameters: ", length(private$theta_names), "\n", sep = "")
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
    #' @param beta Initial beta schedule
    #' @param max_rungs The maximum number of rungs
    #' @param target_acceptance Target acceptance rate
    #' @param silent print progress (boolean)
    tune = function(iterations, beta = seq(1, 0, -0.1), max_rungs = 100, target_acceptance = 0.44, silent = FALSE){
      private$tune_called <- TRUE
      
      # As well as updating internal elements as in burnin,
      # the following will all need to be updated if this function is called:
      ## private$rungs
      ## private$beta
      ## private$swap = FALSE
      ## private$swap_acceptance_counter
      ## private$proposal_sd - to increase ncol to match rungs
      ## private$acceptance_counter  - to increase ncol to match rungs
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
    burn = function(iterations, target_acceptance = 0.44, silent = FALSE){
      stopifnot(is.integer(iterations))
      if(private$sample_called){
        stop("Cannot call burn after sample has been called")
      }
      
      private$burn_called <- TRUE
      private$target_acceptance <- target_acceptance
      
      burnin = TRUE
      for(chain in 1:private$chains){
        raw_output <- mcmc(
          chain,
          burnin,
          iterations,
          silent,
          private$theta[[chain]],
          private$theta_names,
          private$theta_transform_type,
          private$theta_min,
          private$theta_max,
          private$infer_parameter,
          private$data,
          private$loglikelihood,
          private$logprior,
          private$misc,
          private$proposal_sd,
          private$acceptance_counter[[chain]],
          private$target_acceptance,
          private$swap,
          private$beta,
          private$swap_acceptance_counter,
          private$blocks,
          private$n_unique_blocks,
          private$iteration_counter[2, chain],
          private$rng_ptr
        )
        # Check for error return
        if("error" %in% names(raw_output)){
          self$error_debug = raw_output
          stop("Error in mcmc, check $error_debug")
        }
        # Update user output
        private$output_df[[chain]] <- append_output(
          current = private$output_df[[chain]],
          new = raw_output$output,
          phase = "burn",
          theta_names = private$theta_names,
          chain = chain
        )
        # Update internal states
        private$duration[2, chain] = private$duration[2, chain] + raw_output$dur
        private$iteration_counter[2, chain] = private$iteration_counter[2, chain] + iterations
        private$proposal_sd = raw_output$proposal_sd
        private$acceptance_counter[[chain]] = matrix(as.integer(raw_output$acceptance), nrow = length(private$theta_names), ncol = private$rungs)
        private$theta[[chain]] = unlist(private$output_df[[chain]][nrow(private$output_df[[chain]]),private$theta_names])
      }
    },
    
    ### Sampling ###
    #' @description
    #' Run sampling in. Runs the sampling phase of the MCMC.
    #' @param iterations Number of sampling iterations to run
    #' @param silent print progress (boolean)
    sample = function(iterations, silent = FALSE){
      stopifnot(is.integer(iterations))
      
      private$sample_called <- TRUE
      burnin = FALSE
      for(chain in 1:private$chains){
        raw_output <- mcmc(
          chain,
          burnin,
          iterations,
          silent,
          private$theta[[chain]],
          private$theta_names,
          private$theta_transform_type,
          private$theta_min,
          private$theta_max,
          private$infer_parameter,
          private$data,
          private$loglikelihood,
          private$logprior,
          private$misc,
          private$proposal_sd,
          private$acceptance_counter[[chain]],
          private$target_acceptance,
          private$swap,
          private$beta,
          private$swap_acceptance_counter,
          private$blocks,
          private$n_unique_blocks,
          private$iteration_counter[3, chain],
          private$rng_ptr
        )
        # Update user output
        private$output_df[[chain]] <- append_output(
          current = private$output_df[[chain]],
          new = raw_output$output,
          phase = "sample",
          theta_names = private$theta_names,
          chain = chain
        )
        # Update internal states
        private$duration[3, chain] = private$duration[3, chain] + raw_output$dur
        private$iteration_counter[3, chain] = private$iteration_counter[3, chain] + iterations
        private$theta[[chain]] = unlist(private$output_df[[chain]][nrow(private$output_df[[chain]]),private$theta_names])
      }
    },
    
    ### Output ###
    #' @description
    #' Get mcmc output data.frame
    #' @param phase optional phase selection
    output = function(phase = NULL){
      output_df <- dplyr::bind_rows(private$output_df)
      if(!is.null(phase)){
        return(output_df[output_df$phase == phase, ])
      }
      return(output_df)
    },
    
    ### Diagnostics ###
    #' @description
    #' Get acceptance rates
    acceptance_rate = function(){
      estimate_acceptance_rate(
        private$acceptance_counter,
        private$iteration_counter,
        private$chains,
        private$theta_names
      )
    },
    
    #' @description
    #' Get DIC estimate
    dic = function(){
      output_df <- dplyr::bind_rows(private$output_df)
      estimate_dic(output_df)
    },
    
    #' @description
    #' Get effective sample size estimates
    ess = function(){
      output_df <- dplyr::bind_rows(private$output_df)
      estimate_ess(output_df, private$theta_names)
    },
    
    #' @description
    #' Get rhat estimates
    rhat = function(){
      output_df <- dplyr::bind_rows(private$output_df)
      estimate_rhat(output_df, private$theta_names, private$chains, private$iteration_counter)
    },
    
    #' @description
    #' Get run-time data
    timing = function(){
      return(private$duration)
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
