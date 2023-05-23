#' @title Null log likelihood function
#' @param params parmeter vector
#' @param data data list
#' @param misc misc list
#' @export
loglike_null <- function(params, data, misc){
  return(0)
}

#' @title Null log prior function
#' @inheritParams loglike_null
#' @export
logprior_null <- function(params, misc){
  return(0)
}

#' @title Simple normal log likelihood function
#' @inheritParams loglike_null
#' @export
loglike_normal <-  function(params, data, misc) {
  ret <- sum(
    stats::dnorm(
      x = data[["x"]],
      mean = params["mu"],
      sd = params["sigma"],
      log = TRUE
    )
  )
  return(ret)
}

#' @title Double well log likelihood function
#' @inheritParams loglike_null
#' @export
loglike_double_well <-  function(params, data, misc) {
  mu <- params["mu"]
  gamma <- params["gamma"]
  ret <- -gamma * (mu * mu - 1.0) * (mu * mu - 1.0)
  
  return(ret)
}

#' @title Log prior function for return prior check
#' @inheritParams loglike_null
#' @export
logprior_return <- function(params, misc){
 ret <- stats::dnorm(x = params["real_line"], mean = 0.0, sd = 1.0, log = TRUE) +
   stats::dgamma(x = -params["neg_line"], shape = 5.0, scale = 5.0, log = TRUE) +
   stats::dgamma(x = params["pos_line"], shape = 5.0, scale = 5.0, log = TRUE) +
   stats::dbeta(x = params["unit_interval"], shape1 = 3.0, shape2 = 3.0, log = TRUE)
 return(ret)
}

#' @title Blocked log likelihood function
#' @inheritParams loglike_null
#' @export
loglike_block <- function(params, data, misc){
  block <- misc[["block"]]
  if(block == 6){
    out <- sum(stats::dnorm(params[1:5], mean = 0, sd = 1, log = TRUE))
  } else {
    x <- data[[block]]
    out <- sum(stats::dnorm(x, mean = params[block], sd = 1, log = TRUE))
  }
  return(out)
}

