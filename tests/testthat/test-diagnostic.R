#------------------------------------------------
test_that("gelman_rubin works", {
  
  # create dummy posterior draws matrix
  m <- cbind(chain = rep(1:2, each = 1e2), mu = rnorm(2e2))
  
  # expect error when running with single chain or sample
  expect_error(gelman_rubin(m, 1, 10))
  expect_error(gelman_rubin(m, 2, 1))
  
  # expect no warnings when run with correct parameters
  expect_silent(gelman_rubin(m, 2, 10))
})

#------------------------------------------------
test_that("test_convergence works", {
  
  set.seed(1)
  
  # FALSE if single value
  x <- rnorm(1e2)
  expect_false(test_convergence(x, n = 1))
  
  # FALSE if ESS gives warning
  x <- rep(1e100, 100)
  #coda::effectiveSize(x)
  expect_false(test_convergence(x, n = 10))
  
  # FALSE if ESS too small
  x <- rep(0, 100)
  #coda::effectiveSize(x)
  expect_false(test_convergence(x, n = 100))
  
  # expect TRUE when run with appropriate values
  x <- rnorm(1e2)
  expect_true(test_convergence(x, n = 100))
  
})

#------------------------------------------------
test_that("geweke_pvalue works", {
  
  # NaN if coda's geweke.diag() cannot calculate valid output
  expect_true(is.nan(geweke_pvalue(rep(0, 10))))
  
  # expect no warnings when run with appropriate values
  expect_silent(geweke_pvalue(rnorm(10)))
  
})

