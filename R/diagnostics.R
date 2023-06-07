estimate_dic <- function(output){
  deviance <- -2 * output[output$phase == "sample", "loglikelihood"]
  dic  <- mean(deviance) + 0.5 * stats::var(deviance)
  return(dic)
}

estimate_ess <- function(output, parameter_names){
  input <- output[output$phase == "sample", parameter_names, drop = FALSE]
  ess <- apply(input, 2, coda::effectiveSize)
  return(ess)
}

#' Gelman-Rubin statistic
#' 
#' Estimate sthe Gelman-Rubin (rhat) convergence statistic for a single parameter
#' across multiple chains. Basic method, assuming all chains are of equal length
#'  
#' @references Gelman, A., and D. B. Rubin. 1992. 
#' Inference from Iterative Simulation Using Multiple Sequences. 
#' Statistical Science 7: 457–511.
#' @references \url{https://astrostatistics.psu.edu/RLectures/diagnosticsMCMC.pdf}
#'
#' @param par_matrix Matrix (iterations x chains)
#'
#' @return Gelman-Rubin statistic
gelman_rubin <- function(par_matrix){
  # Mean over all samples
  all_mean <- mean(par_matrix)
  
  # Mean of each chain
  chain_mean <- apply(par_matrix, 2, mean)
  
  chains <- ncol(par_matrix)
  samples <- nrow(par_matrix)
  
  # Variance of each chain
  chain_var <- apply(par_matrix, 2, stats::var)
  W <- (1 / chains) * sum(chain_var)
  B <- samples / (chains - 1) * sum((chain_mean - all_mean)^2)
  V <- (1 - 1 / samples) * W + (1 / samples) * B
  round(sqrt(V / W), 4)
}

estimate_rhat <- function(output, pars, chains){
  par_matrices <- lapply(pars, function(x){
    matrix(output[[x]], ncol = chains)
  })
  out <- sapply(par_matrices, gelman_rubin)
  names(out) <- pars
  return(out)
}

estimate_acceptance_rate <- function(
    acceptance_counter,
    iteration_counter,
    chains,
    phases,
    theta_names,
    rungs){
  
  ar <- list()
  counter <- 1
  for(i in chains){
    for(p in phases){
      d <- round(acceptance_counter[[i]][[p]] / iteration_counter[[p]], 3)
      colnames(d) <- theta_names
      ar[[counter]] <- 
        cbind(
          data.frame(
            chain = i,
            phase = p,
            rung = 1:nrow(d)
          ) ,
          d
        )
      counter <- counter + 1
    }
  }
  ar <- list_r_bind(ar)
  
  if(!is.null(rungs)){
    ar <- ar[ar$rung %in% rungs,]
  }
  
  return(ar)
}

estimate_timing <- function(duration, iteration_counter, phases, chains){
  seconds <- list_c_bind(lapply(duration, list_r_bind))
  seconds <- rbind(seconds, colSums(seconds))
  rownames(seconds) <- c(phases, "Total")
  colnames(seconds) <- paste0("Chain_", 1:chains)
  iterations_per_second <- seconds
  for(i in 1:3){
    iterations_per_second[i,] <- iteration_counter[[i]] / iterations_per_second[i,]
  }
  iterations_per_second[4,] <- sum(unlist(iteration_counter)) / iterations_per_second[4,]
  iterations_per_second <- round(iterations_per_second)
  
  return(
    list(
      seconds = seconds,
      iterations_per_second = iterations_per_second)
  )
}

estimate_mc_acceptance_rate <- function(swap_acceptance_counter, iteration_counter, phases, beta){
  mar <- list()
  counter <- 1
  for(p in phases){
    mar[[counter]] <- data.frame(
      phase = p,
      beta_mid = round(beta_mid(beta[[p]]), 3),
      coupling_acceptance = round(swap_acceptance_counter[[p]] / iteration_counter[[p]], 3)
    )
    counter <- counter + 1
  }
  mar <- list_r_bind(mar)
  mar
}


#' @title Estimate autocorrelation
#'
#' @param x samples
#' @param lag lag
#' @export
acf_data <- function(x, lag){
  stats::acf(x, plot = FALSE, lag.max = lag)$acf
}
