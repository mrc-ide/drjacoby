test_that("gelman_rubin works", {
  
  # create dummy posterior draws matrix
  m <- cbind(chain = rep(1:2, each = 1e2), mu = rnorm(2e2))
  
  # expect error when running with single chain or sample
  expect_error(gelman_rubin(m, 1, 10))
  expect_error(gelman_rubin(m, 2, 1))
  
  # expect no warnings when run with correct parameters
  expect_silent(gelman_rubin(m, 2, 10))
})
