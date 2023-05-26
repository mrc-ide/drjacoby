#' Append output from the cpp mcmc to R6 object
#'
#' @param current Any current output (from previous steps)
#' @param new New output from cpp mcmc call
#' @param phase Phase: tune, burn or sample
#' @param theta_names Vector of parameter names
#' @param chain Chain
append_output <- function(current, new, phase, theta_names, chain){
  new <- as.data.frame(new)
  colnames(new) <- c("iteration", theta_names, "logprior", "loglikelihood")
  new$phase <- phase
  new$chain <- chain
  
  if(!is.null(current)){
    # Append iteration count
    new$iteration <- new$iteration + max(current$iteration)
    new <- rbind(current, new)
  }
  return(new)
}

#' Extract names element from nested list
#'
#' @param x Nested list
#' @param name name
chain_element <- function(x, name){
  lapply(x, `[[`, name)  
}

#' Column bind all elements of list
#'
#' @param x List
list_c_bind <- function(x){
  if(!is.data.frame(x)){
    x <- do.call("cbind", x)
  }
  return(x)
}

#' Row bind all elements of list
#'
#' @param x List
list_r_bind <- function(x){
  if(!is.data.frame(x)){
    x <- do.call("rbind", x)
  }
  return(x)
}

create_swap_acceptance_counter = function(chains, rungs){
  lapply(1:chains, function(x){
    matrix(0L, nrow = 3, ncol = max(rungs - 1, 1))
  })
}

create_iteration_counter = function(chains){
  iteration_counter = lapply(1:chains, function(x){
    rep(0L, 3, ncol = 1)
  })
}

create_duration_log = function(chains){
  lapply(1:chains, function(x){
    matrix(0, nrow = 3, ncol = 1)
  })
}

create_acceptance_counter <- function(chains, rungs, n_par){
  lapply(1:chains, function(x){
    array(0L, dim = c(3, rungs, n_par)
    )
  })
}

create_proposal_sd_log <- function(chains, rungs, n_par, init_sd = 0.1){
  lapply(1:chains, function(x){
    matrix(init_sd, nrow = rungs, ncol = n_par)
  })
}
