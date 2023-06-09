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
  out <- lapply(x, `[[`, name)  
  names(out) <- paste0("Chain_", 1:length(out))
  return(out)
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

#' Create nested list of phases in chains
#'
#' @param chains Numebr of chains
#' @param base element at base of list
create_chain_phase_list <- function(chains, base){
  lapply(1:chains, function(x){
    list(
      tune = base,
      burn = base,
      sample = base
    )
  })
}

#' Create list containing proposal sd matrix for each chain
#'
#' @param chains Number of chains
#' @param rungs Number of rungs
#' @param n_par Number of parameters
#' @param init_sd Initial value
create_proposal_sd <- function(chains, rungs, n_par, init_sd = 0.1){
  lapply(1:chains, function(x){
    matrix(init_sd, nrow = rungs, ncol = n_par)
  })
}

#' Get mid points of beta between rungs
#'
#' @param beta Beta
beta_mid <- function(beta){
  beta[- 1] - diff(beta) / 2
}

check_chain <- function(chain, chains){
  present <- chain %in% 1:chains
  if(!all(present)){
    stop("Requested chains not all present in output")
  }
}

check_rung <- function(rung, rungs){
  present <- sapply(rungs, function(x){
    all(rung %in% x)
  })
  if(!all(present)){
    stop("Requested rungs not all present in output")
  }
}
