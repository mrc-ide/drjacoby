#' @title Null log likelihood function
#' @export
loglike_null <- function(params, data, misc){
  return(0)
}

#' @title Null log prior function
#' @export
logprior_null <- function(params, misc){
  return(0)
}

#' @title Simple normal log likelihood function
#' @export
loglike_normal <-  function(params, data, misc) {
  ret <- sum(
    stats::dnorm(
      x = data[["x"]],
      mean = params["mu"],
      sd = 1,
      log = TRUE
    )
  )
  return(ret)
}

#' @title Double well log likelihood function
#' @export
loglike_double_well <-  function(params, data, misc) {
  mu <- params["mu"]
  gamma <- params["gamma"]
  ret <- -gamma * (mu * mu - 1.0) * (mu * mu - 1.0)
  
  return(ret)
}

#' @title Log prior function for return prior check
#' @export
logprior_return <- function(params, misc){
 ret <- dnorm(x = params["real_line"], mean = 0.0, sd = 1.0, log = TRUE) +
    dgamma(x = -params["neg_line"], shape = 5.0, scale = 5.0, log = TRUE) +
    dgamma(x = params["pos_line"], shape = 5.0, scale = 5.0, log = TRUE) +
    dbeta(x = params["unit_interval"], shape1 = 3.0, shape2 = 3.0, log = TRUE)
 return(ret)
}
