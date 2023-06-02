
test_that("plots data subsetting works", {
  phases <- c("tune", "burn", "sample")
  chains <- 1:2
  data = data.frame(
    x = 1:6,
    chain = rep(chains, each = 3),
    phase = rep(phases, 2)
  )
  
  expect_equal(
    plot_data_subset(output_df = data, phase = NULL, chain = NULL),
    data
  )
  for(p in phases){
    expect_equal(
      plot_data_subset(output_df = data, phase = p, chain = NULL),
      data[data$phase == p, ]
    )
  }
  for(c in chains){
    expect_equal(
      plot_data_subset(output_df = data, phase = NULL, chain = c),
      data[data$chain == c, ]
    )
  }
  for(p in phases){
    for(c in chains){
      expect_equal(
        plot_data_subset(output_df = data, phase = p, chain = c),
        data[data$phase == p & data$chain == c, ]
      )
    }
  }
})

test_that("quantile 95 works", {
  qs <- quantile(1:10, c(0.025, 0.5, 0.975))
  names(qs) <- c("Q2.5", "Q50", "Q97.5")
  
  expect_equal(
    quantile_95(1:10), 
    qs
  )
})

test_that("par_plot works", {
  par <- "x"
  phases <- c("tune", "burn", "sample")
  chains <- 1:2
  
  data = data.frame(
    iteration = rep(1:3000, 2),
    x = runif(6000),
    chain = rep(chains, each = 3000),
    phase = rep(phases, 2000)
  )
  lag <- 20
  downsample <- TRUE
  phase <- "sample"
  chain <- NULL
  return_elements <- FALSE
  
  # Return combined plot
  p1 <- create_par_plot(
    par = par,
    output_df = data,
    lag = lag,
    downsample = downsample,
    phase = phase,
    chain = chain,
    return_elements = return_elements
  )
  expect_is(p1, "gg")
  expect_type(p1, "list")
  
  # Return plot elements 
  p2 <- create_par_plot(
    par = par,
    output_df = data,
    lag = lag,
    downsample = downsample,
    phase = phase,
    chain = chain,
    return_elements = TRUE
  )
  expect_is(p2, "list")
  expect_type(p2, "list")
  expect_length(p2, 3)
  expect_is(p2[[1]], "gg")
  expect_is(p2[[2]], "gg")
  expect_is(p2[[3]], "gg")
})


test_that("cor_plot works", {
  phases <- c("tune", "burn", "sample")
  chains <- 1:2
  
  data = data.frame(
    iteration = rep(1:3000, 2),
    x = runif(6000),
    y = runif(6000),
    chain = rep(chains, each = 3000),
    phase = rep(phases, 2000)
  )
  downsample <- TRUE
  phase <- "sample"
  chain <- NULL
  
  p1 <- create_cor_plot(
    parx = "x",
    pary = "y",
    output_df = data,
    downsample = downsample,
    phase = phase,
    chain = chain
  )
  expect_is(p1, "gg")
})

test_that("create_credible works", {
  phases <- c("tune", "burn", "sample")
  chains <- 1:2
  
  data = data.frame(
    iteration = rep(1:3000, 2),
    x = runif(6000),
    y = runif(6000),
    chain = rep(chains, each = 3000),
    phase = rep(phases, 2000)
  )
  phase <- "sample"
  chain <- NULL
  pars <- c("x", "y")
  param_names = pars
  
  p1 <- create_credible_plot(
    pars = c("x", "y"),
    output_df = data,
    phase = phase,
    chain = chain,
    param_names = param_names
  )
  expect_is(p1, "gg")
})

test_that("create_cor_mat works", {
  phases <- c("tune", "burn", "sample")
  chains <- 1:2
  
  data = data.frame(
    iteration = rep(1:3000, 2),
    x = runif(6000),
    y = runif(6000),
    chain = rep(chains, each = 3000),
    phase = rep(phases, 2000)
  )
  phase <- "sample"
  chain <- NULL
  pars <- c("x", "y")
  param_names = pars
  
  p1 <- create_cor_mat_plot(
    output_df = data,
    pars = pars,
    chain = chain,
    phase = phase,
    param_names = param_names
  )
  expect_is(p1, "gg")
})
  
