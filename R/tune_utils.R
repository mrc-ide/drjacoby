propose_new_beta <- function(n, beta_mid, rejection_rate, lambda){
  lambda <- sum(rejection_rate)
  li <- approxfun(
    y = c(0, rev(beta_mid), 1),
    x = c(0, cumsum(rev(rejection_rate)), lambda),
    ties = "ordered"
  )
  rev(li(seq(0, lambda, length.out = n)))
} 
