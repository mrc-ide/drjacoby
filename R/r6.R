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
    error_debug = NULL,
    
    ### Initialisation ###
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
    
    #### Public functions ###
    tune = function(iterations){
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
      
    },
    
    burn = function(iterations, target_acceptance = 0.44, silent = FALSE){
      stopifnot(is.integer(iterations))
      
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
        if("error" %in% names(raw_output)){
          self$error_debug = raw_output
          message("Error in mcmc, check $error_debug")
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
        private$acceptance_counter[[chain]] = private$acceptance_counter[[chain]] + as.integer(raw_output$acceptance)
        private$theta[[chain]] = private$output_df[[chain]][nrow(private$output_df[[chain]]),private$theta_names]
      }
    },
    
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
        private$theta[[chain]] = private$output_df[[chain]][nrow(private$output_df[[chain]]),private$theta_names]
      }
    },
    
    output = function(phase = NULL){
      output_df <- dplyr::bind_rows(private$output_df)
      if(!is.null(phase)){
        return(output_df[output_df$phase == phase, ])
      }
      return(output_df)
    },
    
    acceptance_rate = function(){
      estimate_acceptance_rate(
        private$acceptance_counter,
        private$iteration_counter,
        private$chains,
        private$theta_names
      )
    },
    
    dic = function(){
      output_df <- dplyr::bind_rows(private$output_df)
      estimate_dic(output_df)
    },
    
    ess = function(){
      output_df <- dplyr::bind_rows(private$output_df)
      estimate_ess(output_df, private$theta_names)
    },
    
    rhat = function(){
      output_df <- dplyr::bind_rows(private$output_df)
      estimate_rhat(output_df, private$theta_names, private$chains, private$iteration_counter)
    },
    
    timing = function(){
      return(private$duration)
    },
    
    plot_par = function(par, lag = 20, downsample = TRUE, phase = "sample", chain = NULL){
      create_par_plot(
        par = par,
        output_df = private$output_df,
        lag = lag,
        downsample = downsample,
        phase = phase,
        chain = chain
      )
    }
  )
)
