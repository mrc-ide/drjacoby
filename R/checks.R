loglike_normal <-  function(params, data, misc) {
  ret <- sum(
    stats::dnorm(
      x = data[["x"]],
      mean = params["mu"],
      sd = 1,
      log = TRUE
    )
  )
  return(ret);
}

logprior_null <- function(params, misc){
  return(0)
}
