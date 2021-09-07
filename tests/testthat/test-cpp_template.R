test_that("template works", {
  expect_error(cpp_template(NULL))
  expect_error(cpp_template(1))
  expect_error(cpp_template(save_as = "foo"))
  expect_error(cpp_template(save_as = "foo.csv"))
})

