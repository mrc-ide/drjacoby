run_mcmc <- function(chain_input, burnin, phase, iterations, silent, theta_names, theta_transform_type,
                     theta_min, theta_max, infer_parameter, data, loglikelihood,
                     logprior, misc, target_acceptance, swap,
                     beta, blocks, n_unique_blocks){
  raw_output <- mcmc(
    chain_input$chain,
    burnin,
    iterations,
    silent,
    chain_input$theta,
    theta_names,
    theta_transform_type,
    theta_min,
    theta_max,
    infer_parameter,
    data,
    loglikelihood,
    logprior,
    misc,
    chain_input$proposal_sd,
    chain_input$acceptance_counter[[phase]],
    target_acceptance,
    swap,
    beta,
    chain_input$swap_acceptance_counter,
    blocks,
    n_unique_blocks,
    chain_input$iteration_counter[phase],
    chain_input$rng_ptr
  )
  chain_input$rng_ptr$sync()
  return(raw_output)
}

# lapply set up with same argument as parallel::clusterApply
noclusterApply <- function(cl, x, fun, ...){
  lapply(X = x, FUN = fun, ...)
}
