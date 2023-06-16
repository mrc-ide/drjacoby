set.seed(1)

rlang::is_installed("Rcpp")
rlang::is_installed("rogress")

n <- 50
iterations <- 100
data <- list(x = rnorm(n))

df_params <- data.frame() |>
  add_parameter(name = "mu", min = -1, max = 1) |>
  add_parameter(name = "sigma", min = 0, max = 1)

fpath <- system.file("extdata", "rcpp_likelihood_and_prior_functions.cpp", package = "drjacoby")
Rcpp::sourceCpp(fpath)

test_that("Chains = 1 Rungs = No, Functions = cpp11", {
  set.seed(1)
  mcmc <- dj$new( 
    data = data,
    df_params = df_params,
    loglikelihood = loglike_normal_rcpp,
    logprior = logprior_null_rcpp,
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
  
  expect_error(
    mcmc$output(phase = "foo"),
    "Requested phases not all present in output",
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
  
  expect_error(
    mcmc$acceptance_rate(phase = "foo"),
    "Requested phases not all present in output",
    fixed = TRUE
  )
  
  # MC acceptance rate
  expect_error(
    mcmc$mc_acceptance_rate(),
    "Not available for a single rung",
    fixed = TRUE
  )
  
  expect_error(
    mcmc$mc_acceptance_rate(phase = "burn"),
    "Not available for a single rung",
    fixed = TRUE
  )
  
  expect_error(
    mcmc$mc_acceptance_rate(phase = c("burn", "sample")),
    "Not available for a single rung",
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
  
  # Beta
  beta <- mcmc$get_beta()
  expect_type(beta, "list")
  expect_equal(length(beta), 3)
  expect_equal(sum(is.na(unlist(beta))), 0)
  
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
    mcmc$plot_par("foo"),
    "Requested parameters not all present in output",
    fixed = TRUE
  )
  
  expect_error(
    mcmc$plot_par(par = "mu", chain = 2),
    "Requested chains not all present in output",
    fixed = TRUE
  )
  
  # Plot correlation
  p <- mcmc$plot_cor(parx = "mu", pary = "sigma")
  expect_type(p, "list")
  expect_equal(is(p), "gg")
  
  expect_error(
    mcmc$plot_cor(parx = "foo", pary = "sigma"),
    "Requested parameters not all present in output",
    fixed = TRUE
  )
  
  expect_error(
    mcmc$plot_cor(parx = "mu", pary = "foo"),
    "Requested parameters not all present in output",
    fixed = TRUE
  )
  
  expect_error(
    mcmc$plot_cor(parx = "mu", pary = "sigma", chain = 2),
    "Requested chains not all present in output",
    fixed = TRUE
  )
  
  # Plot credible
  p <- mcmc$plot_credible()
  expect_type(p, "list")
  expect_equal(is(p), "gg")
  
  p <- mcmc$plot_credible(pars = "mu")
  expect_type(p, "list")
  expect_equal(is(p), "gg")
  
  expect_error(
    mcmc$plot_credible(pars = "foo"),
    "Requested parameters not all present in output",
    fixed = TRUE
  )
  
  expect_error(
    mcmc$plot_credible(pars = c("mu", "foo")),
    "Requested parameters not all present in output",
    fixed = TRUE
  )
  
  p <- mcmc$plot_credible(phase = "burn")
  expect_type(p, "list")
  expect_equal(is(p), "gg")
  
  expect_error(
    mcmc$plot_credible(phase = "foo"),
    "Requested phases not all present in output",
    fixed = TRUE
  )
  
  expect_error(
    mcmc$plot_credible(pars = "mu", param_names = c("foo", "bar")),
    "length(param_names) == length(pars) is not TRUE",
    fixed = TRUE
  )
  
  # Plot cor mat
  p <- mcmc$plot_cor_mat()
  expect_type(p, "list")
  expect_equal(is(p), "gg")
  
  p <- mcmc$plot_cor_mat(pars = c("mu", "sigma"))
  expect_type(p, "list")
  expect_equal(is(p), "gg")
  
  expect_error(
    mcmc$plot_cor_mat(pars = "mu"),
    "length(pars) > 1 is not TRUE",
    fixed = TRUE
  )
  
  expect_error(
    mcmc$plot_cor_mat(pars = c("foo", "bar")),
    "Requested parameters not all present in output",
    fixed = TRUE
  )
  
  expect_error(
    mcmc$plot_cor_mat(pars = c("mu", "foo")),
    "Requested parameters not all present in output",
    fixed = TRUE
  )
  
  p <- mcmc$plot_cor_mat(phase = "burn")
  expect_type(p, "list")
  expect_equal(is(p), "gg")
  
  expect_error(
    mcmc$plot_cor_mat(phase = "foo"),
    "Requested phases not all present in output",
    fixed = TRUE
  )
  
  expect_error(
    mcmc$plot_cor_mat(pars = "mu", param_names = c("foo", "bar")),
    "length(param_names) == length(pars) is not TRUE",
    fixed = TRUE
  )
  
  # plot mc acceptance_rate
  expect_error(
    mcmc$plot_mc_acceptance_rate(),
    "Not available for a single rung",
    fixed = TRUE
  )
  
  # plot communication barrier
  expect_error(
    mcmc$plot_local_communication_barrier(),
    "Not available for a single rung",
    fixed = TRUE
  )
})

test_that("Chains > 1 Rungs = No, Functions = cpp11", {
  set.seed(1)
  mcmc <- dj$new( 
    data = data,
    df_params = df_params,
    loglikelihood = loglike_normal_rcpp,
    logprior = logprior_null_rcpp,
    chains = 3
  )
  mcmc$burn(iterations = iterations, silent = TRUE)
  mcmc$burn(iterations = iterations, silent = TRUE)
  mcmc$sample(iterations = iterations, silent = TRUE)
  mcmc$sample(iterations = iterations, silent = TRUE)
  
  # Output
  output <- mcmc$output()
  expect_type(output, "list")
  expect_equal(dim(output), c(iterations * 4 * 3, 7))
  expect_named(output, c("iteration", "mu", "sigma", "logprior", "loglikelihood", "phase", "chain"))
  expect_equal(sum(is.na(output)), 0)
  
  output <- mcmc$output(phase = "burn")
  expect_type(output, "list")
  expect_equal(dim(output), c(iterations * 2 * 3, 7))
  expect_named(output, c("iteration", "mu", "sigma", "logprior", "loglikelihood", "phase", "chain"))
  expect_equal(sum(is.na(output)), 0)
  
  output <- mcmc$output(phase = c("burn", "sample"))
  expect_type(output, "list")
  expect_equal(dim(output), c(iterations * 4 * 3, 7))
  expect_named(output, c("iteration", "mu", "sigma", "logprior", "loglikelihood", "phase", "chain"))
  expect_equal(sum(is.na(output)), 0)
  
  output <- mcmc$output(chain = 2)
  expect_type(output, "list")
  expect_equal(dim(output), c(iterations * 4 * 1, 7))
  expect_named(output, c("iteration", "mu", "sigma", "logprior", "loglikelihood", "phase", "chain"))
  expect_equal(sum(is.na(output)), 0)
  
  expect_error(
    output <- mcmc$output(chain = 4),
    "Requested chains not all present in output",
    fixed = TRUE
  )
  
  expect_error(
    mcmc$output(phase = "foo"),
    "Requested phases not all present in output",
    fixed = TRUE
  )
  
  # Acceptance rate
  ar <- mcmc$acceptance_rate()
  expect_type(ar, "list")
  expect_equal(dim(ar), c(3, 5))
  expect_named(ar, c("chain", "phase", "rung", "mu", "sigma"))
  expect_equal(sum(is.na(ar)), 0)
  
  ar <- mcmc$acceptance_rate(phase = "burn")
  expect_type(ar, "list")
  expect_equal(dim(ar), c(3, 5))
  expect_named(ar, c("chain", "phase", "rung", "mu", "sigma"))
  expect_equal(sum(is.na(ar)), 0)
  
  ar <- mcmc$acceptance_rate(phase = c("burn", "sample"))
  expect_type(ar, "list")
  expect_equal(dim(ar), c(6, 5))
  expect_named(ar, c("chain", "phase", "rung", "mu", "sigma"))
  expect_equal(sum(is.na(ar)), 0)
  
  ar <- mcmc$acceptance_rate(chain = 2)
  expect_type(ar, "list")
  expect_equal(dim(ar), c(1, 5))
  expect_named(ar, c("chain", "phase", "rung", "mu", "sigma"))
  expect_equal(sum(is.na(ar)), 0)
  
  expect_error(
    mcmc$acceptance_rate(chain = 4),
    "Requested chains not all present in output",
    fixed = TRUE
  )
  
  expect_error(
    mcmc$acceptance_rate(chain = 1:5),
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
  
  expect_error(
    mcmc$acceptance_rate(phase = "foo"),
    "Requested phases not all present in output",
    fixed = TRUE
  )
  
  # MC acceptance rate
  expect_error(
    mcmc$mc_acceptance_rate(),
    "Not available for a single rung",
    fixed = TRUE
  )
  
  expect_error(
    mcmc$mc_acceptance_rate(phase = "burn"),
    "Not available for a single rung",
    fixed = TRUE
  )
  
  expect_error(
    mcmc$mc_acceptance_rate(phase = c("burn", "sample")),
    "Not available for a single rung",
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
  rhat <- mcmc$rhat()
  expect_type(rhat, "double")
  expect_length(rhat, 2)
  expect_named(rhat, c("mu", "sigma"))
  expect_equal(sum(is.na(rhat)), 0)
  
  # Timing
  timing <- mcmc$timing()
  expect_type(timing, "list")
  expect_length(timing, 2)
  expect_named(timing, c("seconds", "iterations_per_second"))
  expect_equal(sum(is.na(unlist(timing))), 0)
  
  # Beta
  beta <- mcmc$get_beta()
  expect_type(beta, "list")
  expect_equal(length(beta), 3)
  expect_equal(sum(is.na(unlist(beta))), 0)
  
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
    mcmc$plot_par("foo"),
    "Requested parameters not all present in output",
    fixed = TRUE
  )
  
  p <- mcmc$plot_par(par = "mu", chain = 2)
  expect_type(p, "list")
  expect_equal(is(p), "gg")
  
  expect_error(
    mcmc$plot_par(par = "mu", chain = 4),
    "Requested chains not all present in output",
    fixed = TRUE
  )
  
  # Plot correlation
  p <- mcmc$plot_cor(parx = "mu", pary = "sigma")
  expect_type(p, "list")
  expect_equal(is(p), "gg")
  
  expect_error(
    mcmc$plot_cor(parx = "foo", pary = "sigma"),
    "Requested parameters not all present in output",
    fixed = TRUE
  )
  
  expect_error(
    mcmc$plot_cor(parx = "mu", pary = "foo"),
    "Requested parameters not all present in output",
    fixed = TRUE
  )
  
  p <- mcmc$plot_cor(parx = "mu", pary = "sigma", chain = 2)
  expect_type(p, "list")
  expect_equal(is(p), "gg")
  
  expect_error(
    mcmc$plot_cor(parx = "mu", pary = "sigma", chain = 4),
    "Requested chains not all present in output",
    fixed = TRUE
  )
  
  # Plot credible
  p <- mcmc$plot_credible()
  expect_type(p, "list")
  expect_equal(is(p), "gg")
  
  p <- mcmc$plot_credible(pars = "mu")
  expect_type(p, "list")
  expect_equal(is(p), "gg")
  
  expect_error(
    mcmc$plot_credible(pars = "foo"),
    "Requested parameters not all present in output",
    fixed = TRUE
  )
  
  expect_error(
    mcmc$plot_credible(pars = c("mu", "foo")),
    "Requested parameters not all present in output",
    fixed = TRUE
  )
  
  p <- mcmc$plot_credible(phase = "burn")
  expect_type(p, "list")
  expect_equal(is(p), "gg")
  
  expect_error(
    mcmc$plot_credible(phase = "foo"),
    "Requested phases not all present in output",
    fixed = TRUE
  )
  
  expect_error(
    mcmc$plot_credible(pars = "mu", param_names = c("foo", "bar")),
    "length(param_names) == length(pars) is not TRUE",
    fixed = TRUE
  )
  
  # Plot cor mat
  p <- mcmc$plot_cor_mat()
  expect_type(p, "list")
  expect_equal(is(p), "gg")
  
  p <- mcmc$plot_cor_mat(pars = c("mu", "sigma"))
  expect_type(p, "list")
  expect_equal(is(p), "gg")
  
  expect_error(
    mcmc$plot_cor_mat(pars = "mu"),
    "length(pars) > 1 is not TRUE",
    fixed = TRUE
  )
  
  expect_error(
    mcmc$plot_cor_mat(pars = c("foo", "bar")),
    "Requested parameters not all present in output",
    fixed = TRUE
  )
  
  expect_error(
    mcmc$plot_cor_mat(pars = c("mu", "foo")),
    "Requested parameters not all present in output",
    fixed = TRUE
  )
  
  p <- mcmc$plot_cor_mat(phase = "burn")
  expect_type(p, "list")
  expect_equal(is(p), "gg")
  
  expect_error(
    mcmc$plot_cor_mat(phase = "foo"),
    "Requested phases not all present in output",
    fixed = TRUE
  )
  
  expect_error(
    mcmc$plot_cor_mat(pars = "mu", param_names = c("foo", "bar")),
    "length(param_names) == length(pars) is not TRUE",
    fixed = TRUE
  )
  
  # plot mc acceptance_rate
  expect_error(
    mcmc$plot_mc_acceptance_rate(),
    "Not available for a single rung",
    fixed = TRUE
  )
  
  # plot communication barrier
  expect_error(
    mcmc$plot_local_communication_barrier(),
    "Not available for a single rung",
    fixed = TRUE
  )
})

test_that("Chains = 1 Rungs = Yes, Functions = cpp11", {
  set.seed(1)
  mcmc <- dj$new( 
    data = data,
    df_params = df_params,
    loglikelihood = loglike_normal_rcpp,
    logprior = logprior_null_rcpp,
    chains = 1
  )
  mcmc$tune(iterations = iterations, silent = TRUE)
  mcmc$burn(iterations = iterations, silent = TRUE)
  mcmc$burn(iterations = iterations, silent = TRUE)
  mcmc$sample(iterations = iterations, silent = TRUE)
  mcmc$sample(iterations = iterations, silent = TRUE)
  
  # Output
  output <- mcmc$output()
  expect_type(output, "list")
  expect_equal(dim(output), c(iterations * 4 * 1, 7))
  expect_named(output, c("iteration", "mu", "sigma", "logprior", "loglikelihood", "phase", "chain"))
  expect_equal(sum(is.na(output)), 0)
  
  output <- mcmc$output(phase = "burn")
  expect_type(output, "list")
  expect_equal(dim(output), c(iterations * 2 * 1, 7))
  expect_named(output, c("iteration", "mu", "sigma", "logprior", "loglikelihood", "phase", "chain"))
  expect_equal(sum(is.na(output)), 0)
  
  output <- mcmc$output(phase = c("burn", "sample"))
  expect_type(output, "list")
  expect_equal(dim(output), c(iterations * 4 * 1, 7))
  expect_named(output, c("iteration", "mu", "sigma", "logprior", "loglikelihood", "phase", "chain"))
  expect_equal(sum(is.na(output)), 0)
  
  expect_error(
    mcmc$output(chain = 2),
    "Requested chains not all present in output",
    fixed = TRUE
  )
  
  expect_error(
    mcmc$output(phase = "foo"),
    "Requested phases not all present in output",
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
    ar <- mcmc$acceptance_rate(chain = 1:5),
    "Requested chains not all present in output",
    fixed = TRUE
  )
  
  ar <- mcmc$acceptance_rate(rung = 2)
  expect_type(ar, "list")
  expect_equal(dim(ar), c(1, 5))
  expect_named(ar, c("chain", "phase", "rung", "mu", "sigma"))
  expect_equal(sum(is.na(ar)), 0)
  
  ar <- mcmc$acceptance_rate(rung = 1:2)
  expect_type(ar, "list")
  expect_equal(dim(ar), c(2, 5))
  expect_named(ar, c("chain", "phase", "rung", "mu", "sigma"))
  expect_equal(sum(is.na(ar)), 0)
  
  expect_error(
    mcmc$acceptance_rate(phase = "foo"),
    "Requested phases not all present in output",
    fixed = TRUE
  )
  
  # MC acceptance rate
  mc_ar <- mcmc$mc_acceptance_rate()
  expect_type(mc_ar, "list")
  expect_equal(dim(mc_ar), c(59, 3))
  expect_named(mc_ar, c("phase", "beta_mid", "coupling_acceptance"))
  expect_equal(sum(is.na(mc_ar)), 0)
  
  mc_ar <- mcmc$mc_acceptance_rate(phase = "burn")
  expect_type(mc_ar, "list")
  expect_equal(dim(mc_ar), c(5, 3))
  expect_named(mc_ar, c("phase", "beta_mid", "coupling_acceptance"))
  expect_equal(sum(is.na(mc_ar)), 0)
  
  mc_ar <- mcmc$mc_acceptance_rate(phase = c("burn", "sample"))
  expect_type(mc_ar, "list")
  expect_equal(dim(mc_ar), c(10, 3))
  expect_named(mc_ar, c("phase", "beta_mid", "coupling_acceptance"))
  expect_equal(sum(is.na(mc_ar)), 0)
  
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
  
  # Beta
  beta <- mcmc$get_beta()
  expect_type(beta, "list")
  expect_equal(length(beta), 3)
  expect_equal(sum(is.na(unlist(beta))), 0)
  
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
    mcmc$plot_par("foo"),
    "Requested parameters not all present in output",
    fixed = TRUE
  )
  
  expect_error(
    mcmc$plot_par(par = "mu", chain = 2),
    "Requested chains not all present in output",
    fixed = TRUE
  )
  
  # Plot correlation
  p <- mcmc$plot_cor(parx = "mu", pary = "sigma")
  expect_type(p, "list")
  expect_equal(is(p), "gg")
  
  expect_error(
    mcmc$plot_cor(parx = "foo", pary = "sigma"),
    "Requested parameters not all present in output",
    fixed = TRUE
  )
  
  expect_error(
    mcmc$plot_cor(parx = "mu", pary = "foo"),
    "Requested parameters not all present in output",
    fixed = TRUE
  )
  
  expect_error(
    mcmc$plot_cor(parx = "mu", pary = "sigma", chain = 2),
    "Requested chains not all present in output",
    fixed = TRUE
  )
  
  # Plot credible
  p <- mcmc$plot_credible()
  expect_type(p, "list")
  expect_equal(is(p), "gg")
  
  p <- mcmc$plot_credible(pars = "mu")
  expect_type(p, "list")
  expect_equal(is(p), "gg")
  
  expect_error(
    mcmc$plot_credible(pars = "foo"),
    "Requested parameters not all present in output",
    fixed = TRUE
  )
  
  expect_error(
    mcmc$plot_credible(pars = c("mu", "foo")),
    "Requested parameters not all present in output",
    fixed = TRUE
  )
  
  p <- mcmc$plot_credible(phase = "burn")
  expect_type(p, "list")
  expect_equal(is(p), "gg")
  
  expect_error(
    mcmc$plot_credible(phase = "foo"),
    "Requested phases not all present in output",
    fixed = TRUE
  )
  
  expect_error(
    mcmc$plot_credible(pars = "mu", param_names = c("foo", "bar")),
    "length(param_names) == length(pars) is not TRUE",
    fixed = TRUE
  )
  
  # Plot cor mat
  p <- mcmc$plot_cor_mat()
  expect_type(p, "list")
  expect_equal(is(p), "gg")
  
  p <- mcmc$plot_cor_mat(pars = c("mu", "sigma"))
  expect_type(p, "list")
  expect_equal(is(p), "gg")
  
  expect_error(
    mcmc$plot_cor_mat(pars = "mu"),
    "length(pars) > 1 is not TRUE",
    fixed = TRUE
  )
  
  expect_error(
    mcmc$plot_cor_mat(pars = c("foo", "bar")),
    "Requested parameters not all present in output",
    fixed = TRUE
  )
  
  expect_error(
    mcmc$plot_cor_mat(pars = c("mu", "foo")),
    "Requested parameters not all present in output",
    fixed = TRUE
  )
  
  p <- mcmc$plot_cor_mat(phase = "burn")
  expect_type(p, "list")
  expect_equal(is(p), "gg")
  
  expect_error(
    mcmc$plot_cor_mat(phase = "foo"),
    "Requested phases not all present in output",
    fixed = TRUE
  )
  
  expect_error(
    mcmc$plot_cor_mat(pars = "mu", param_names = c("foo", "bar")),
    "length(param_names) == length(pars) is not TRUE",
    fixed = TRUE
  )
  
  # plot mc acceptance_rate
  p <- mcmc$plot_mc_acceptance_rate()
  expect_type(p, "list")
  expect_equal(is(p), "gg")
  
  
  # plot communication barrier
  p <- mcmc$plot_local_communication_barrier()
  expect_type(p, "list")
  expect_equal(is(p), "gg")
})

