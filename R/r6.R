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
    
    chains = NULL,
    
    blocks = NULL,
    n_unique_blocks = NULL,
    
    proposal_sd = NULL,
    target_acceptance = NULL,
    
    rungs = 1,
    n_rungs = 1,
    beta = 1,
    
    output = NULL
  ),
  
  public = list(
    ### Public variables ###

    ### Initialisation ###
    initialize = function(data, df_params, loglikelihood, logprior, chains = 1L, misc = list()){
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
      private$proposal_sd = rep(1, length(private$theta))
      
      private$chains = chains
      
      private$blocks = set_blocks(df_params)
      private$n_unique_blocks = length(unique(unlist(private$blocks)))
    },
    
    #### Public functions ###
    tune = function(iterations){
      private$tune_called <- TRUE
      # For Bob
    },
    
    burn = function(iterations, target_acceptance = 0.44){
      private$burn_called <- TRUE
      private$target_acceptance <- target_acceptance
      
      burnin = TRUE
      raw_output <- mcmc(
        chain,
        burnin,
        iterations,
        ## iteration_counter_init,
        ## silent,
        private$theta,
        private$theta_names,
        ##get_transform_type(unlist(df_params$min), unlist(df_params$max)),
        private$theta_min,
        private$theta_max,
        ##infer_parameter,
        private$data,
        private$loglikelihood,
        private$logprior,
        private$misc,
        private$proposal_sd,
        # acceptance_init,
        private$target_acceptance,
        #swap,
        #beta_init,
        #swap_acceptance_init,
        private$blocks,
        private$n_unique_blocks
      )
      
    },
    
    sample = function(iterations){
      private$sample_called <- TRUE
      
    },
    
    get_output = function(phase = NULL){
      if(!is.null(phase)){
        return(private$output[private$output$phase == phase, ])
      }
      return(private$output)
    }
  )
)
