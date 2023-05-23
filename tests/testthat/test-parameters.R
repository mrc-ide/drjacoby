test_that("define_params() working as expected", {
  
  # expect error if bad name
  expect_error(define_params(name = "mu", min = -10, max = 10, foo = 5))
  
  # expect error if duplicated names
  expect_error(define_params(name = "mu", min = -10, max = 10,
                             name = "mu", min = -10, max = 10))
  
  # expect error if missing some args
  expect_error(define_params(name = "mu1", min = -10, max = 10, init = 5,
                             name = "mu2", min = -10, max = 10))
  
  # expect blocks in correct format
  df_params <- define_params(name = "mu", min = -10, max = 10, block = 1,
                             name = "sigma", min = 0, max = Inf, block = 1)
  expect_equal(df_params$block, list(block = 1, block = 1))
  
  # correctly defined dataframe
  df_params <- define_params(name = "mu", min = -10, max = 10, init = 5,
                             name = "sigma", min = 0, max = Inf, init = 1)
  
  # same parameters dataframe using base R method
  df_params_base <- data.frame(name = c("mu", "sigma"),
                               min = c(-10, 0),
                               max = c(10, Inf))
  df_params_base$init <- list(init = 5, init = 1)
  
  # check identical
  expect_identical(df_params, df_params_base)
  
})
