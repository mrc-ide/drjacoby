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
    
    tune_iterations = 0L,
    burn_iterations = 0L,
    sample_iterations = 0L,
    
    output_df = NULL,
    rng_ptr = NULL
  ),
  
  public = list(
    ### Public variables ###
    error_debug = NULL,
    
    ### Initialisation ###
    initialize = function(data, df_params, loglikelihood, logprior, chains = 1L, misc = list(), seed = NULL){
      stopifnot(is.list(data))
      stopifnot(is.list(df_params))
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
      private$proposal_sd = matrix(0.1, nrow = length(private$theta), ncol = private$rungs)
      private$infer_parameter = as.integer(!(private$theta_min == private$theta_max))
      
      private$chains = chains
      
      private$blocks = set_blocks(df_params)
      private$n_unique_blocks = length(unique(unlist(private$blocks)))
      
      private$acceptance_counter = matrix(0L, nrow = length(private$theta), ncol = private$rungs)
      
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
        ## private$tune_iterations - keep track of total iterations for tuning phase
      
    },
    
    burn = function(iterations, target_acceptance = 0.44, silent = FALSE){
      stopifnot(is.integer(iterations))
      
      private$burn_called <- TRUE
      private$target_acceptance <- target_acceptance
      
      burnin = TRUE
      chain = 1L ## TODO loop over chains
      
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
        private$acceptance_counter,
        private$target_acceptance,
        private$swap,
        private$beta,
        private$swap_acceptance_counter,
        private$blocks,
        private$n_unique_blocks,
        private$burn_iterations,
        private$rng_ptr
      )
      if("error" %in% names(raw_output)){
        self$error_debug = raw_output
        message("Error in mcmc, check $error_debug")
      }
      # Update user output
      private$output_df <- append_output(
        current = private$output_df,
        new = raw_output$output,
        phase = "burn",
        theta_names = private$theta_names
      )
      # Update internal states
      private$burn_iterations = private$burn_iterations + iterations
      private$proposal_sd = raw_output$proposal_sd
      private$acceptance_counter = as.integer(raw_output$acceptance)
      private$theta[[chain]] = private$output_df[nrow(private$output_df),private$theta_names]
    },
    
    sample = function(iterations, silent = FALSE){
      stopifnot(is.integer(iterations))
      
      private$sample_called <- TRUE
      
      burnin = FALSE
      chain = 1L ## TODO loop over chains
      
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
        private$acceptance_counter,
        private$target_acceptance,
        private$swap,
        private$beta,
        private$swap_acceptance_counter,
        private$blocks,
        private$n_unique_blocks,
        private$sample_iterations,
        private$rng_ptr
      )
      # Update user output
      private$output_df <- append_output(
        current = private$output_df,
        new = raw_output$output,
        phase = "sample",
        theta_names = private$theta_names
      )
      # Update internal states
      private$sample_iterations = private$sample_iterations + iterations
      private$theta[[chain]] = private$output_df[nrow(private$output_df),private$theta_names]
      
    },
    
    output = function(phase = NULL){
      if(!is.null(phase)){
        return(private$output_df[private$output_df$phase == phase, ])
      }
      return(private$output_df)
    },
    
    acceptance_rate = function(){
      private$acceptance_counter / private$burn_iterations
    },
    
    dic = function(){
      estimate_dic(private$output_df)
    },
    
    ess = function(){
      estimate_ess(private$output_df, private$theta_names)
    }
  )
)
