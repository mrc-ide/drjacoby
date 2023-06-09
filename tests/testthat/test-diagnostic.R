test_that("gelman_rubin works", {
  
  set.seed(1)
  m <- matrix(runif(300), ncol = 3)
  expect_lt(gelman_rubin(m), 1.1)
  m <- matrix(c(
    runif(100, 0, 1),
    runif(100, 10, 11),
    runif(100, 20, 21)
  ), ncol = 3)
  expect_gt(gelman_rubin(m), 1.1)
})
