set.seed(1)

n <- 50
iterations <- 100
data <- list(x = rnorm(n))

df_params <- data.frame() |>
  add_parameter(name = "mu", min = -1, max = 1) |>
  add_parameter(name = "sigma", min = 0, max = 1)

## Want to test all combinations of:
  # Chains 1 or multiple
  # Rungs on or off
  # Function R/cpp11/Rcpp
  # x <- expand.grid(chains = c("One", "Multiple"), rungs = c("Yes", "No"), f = c("R", "cpp11", "Rcpp"))
  # x <- x[!(x$chain == "Multiple" & x$rungs == "Yes"),]

test_that("Chains = 1 Rungs = No, Functions = R", {
  mcmc <- dj$new( 
    data = data,
    df_params = df_params,
    loglikelihood = loglike_normal,
    logprior = logprior_null,
    chains = 1
  )
  mcmc$burn(iterations = iterations, silent = TRUE)
  mcmc$burn(iterations = iterations, silent = TRUE)
  mcmc$sample(iterations = iterations, silent = TRUE)
  mcmc$sample(iterations = iterations, silent = TRUE)
  
  # Output
  output <- mcmc$output()
  expect_type(output, "list")
  expect_equal(dim(output), c(iterations * 4, 7))
  expect_named(output, c("iteration", "mu", "sigma", "logprior", "loglikelihood", "phase", "chain"))
  expect_equal(sum(is.na(output)), 0)
  
  output <- mcmc$output(phase = "burn")
  expect_type(output, "list")
  expect_equal(dim(output), c(iterations * 2, 7))
  expect_named(output, c("iteration", "mu", "sigma", "logprior", "loglikelihood", "phase", "chain"))
  expect_equal(sum(is.na(output)), 0)
  
  output <- mcmc$output(phase = c("burn", "sample"))
  expect_type(output, "list")
  expect_equal(dim(output), c(iterations * 4, 7))
  expect_named(output, c("iteration", "mu", "sigma", "logprior", "loglikelihood", "phase", "chain"))
  expect_equal(sum(is.na(output)), 0)
  
  expect_error(
    mcmc$output(chain = 2),
    "Requested chains not all present in output",
    fixed = TRUE
  )
  
  # Acceptance rate
  ar <- mcmc$acceptance_rate()
  expect_type(ar, "list")
  expect_equal(dim(ar), c(1, 5))
  expect_named(ar, c("chain", "phase", "rung", "mu", "sigma"))
  expect_equal(sum(is.na(ar)), 0)
  
  ar <- mcmc$acceptance_rate(phase = "burn")
  expect_type(ar, "list")
  expect_equal(dim(ar), c(1, 5))
  expect_named(ar, c("chain", "phase", "rung", "mu", "sigma"))
  expect_equal(sum(is.na(ar)), 0)
  
  ar <- mcmc$acceptance_rate(phase = c("burn", "sample"))
  expect_type(ar, "list")
  expect_equal(dim(ar), c(2, 5))
  expect_named(ar, c("chain", "phase", "rung", "mu", "sigma"))
  expect_equal(sum(is.na(ar)), 0)
  
  expect_error(
    mcmc$acceptance_rate(chain = 2),
    "Requested chains not all present in output",
    fixed = TRUE
  )
  
  expect_error(
    mcmc$acceptance_rate(chain = 1:2),
    "Requested chains not all present in output",
    fixed = TRUE
  )
  
  expect_error(
    mcmc$acceptance_rate(rung = 2),
    "Requested rungs not all present in output",
    fixed = TRUE
  )
  
  expect_error(
    mcmc$acceptance_rate(rung = 1:2),
    "Requested rungs not all present in output",
    fixed = TRUE
  )

  # MC acceptance rate
  expect_error(
    mcmc$mc_acceptance_rate(),
    "private$tune_called is not TRUE",
    fixed = TRUE
  )
  
  expect_error(
    mcmc$mc_acceptance_rate(phase = "burn"),
    "private$tune_called is not TRUE",
    fixed = TRUE
  )
  
  expect_error(
    mcmc$mc_acceptance_rate(phase = c("burn", "sample")),
    "private$tune_called is not TRUE",
    fixed = TRUE
  )
  
  # DIC
  dic <- mcmc$dic()
  expect_type(dic, "double")
  expect_length(dic, 1)
  expect_equal(sum(is.na(dic)), 0)
  
  # ESS
  ess <- mcmc$ess()
  expect_type(ess, "double")
  expect_length(ess, 2)
  expect_equal(sum(is.na(ess)), 0)
  
  # Rhat
  expect_error(
    mcmc$rhat(),
    "Rhat not applicable for a single chain",
    fixed = TRUE
  )
  
  # Timing
  timing <- mcmc$timing()
  expect_type(timing, "list")
  expect_length(timing, 2)
  expect_named(timing, c("seconds", "iterations_per_second"))
  expect_equal(sum(is.na(unlist(timing))), 0)
  
  # Plot parameters
  p <- mcmc$plot_par(par = "mu")
  expect_type(p, "list")
  expect_equal(is(p), "gg")
  
  p <- mcmc$plot_par("mu", return_elements = TRUE)
  expect_type(p, "list")
  for(i in 1:3){
    expect_equal(is(p[[i]]), "gg")
  }
  
  expect_error(
    mcmc$plot_par("foo", return_elements = TRUE),
    "Parameter not present",
    fixed = TRUE
  )
})
