test_that("cpp file scraping  works", {
  sourceCpp("test_input_files/sourced_likelihood.cpp")
  expect_equal(loglike(c(0, 2), c(0,0, -1, 1)),  -6.698343, tol = 0.00001)
  
  t1 <- cpp_function_get("test_input_files/sourced_likelihood.cpp")
  expect_null(drjacoby::check_likelihood_compilation(t1))
  expect_error(cpp_function_get("badname.cqq"), "The file at address should be a .cpp file")
  
  expect_error(cpp_function_get("test_input_files/sourced_likelihood_2.cpp"),
               "The string: '// drjacoby_function_start' was found more than once in the specified file")
  expect_error(cpp_function_get("test_input_files/sourced_likelihood_3.cpp"),
               "The string: '// drjacoby_function_start' could not be found in the specified file")
  expect_error(cpp_function_get("test_input_files/sourced_likelihood_4.cpp"),
               "The string: '// drjacoby_function_end' could not be found in the specified file. This
         must be included on the final line to flag the function end")
})
