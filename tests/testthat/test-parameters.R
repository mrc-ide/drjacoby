test_that("add_parameter() working as expected", {
  expect_error(
    add_parameter(1),
    "is.data.frame(x) is not TRUE",
    fixed = TRUE
  )
  
  expect_error(
    add_parameter(data.frame(), name = 1),
    "is.character(name) is not TRUE",
    fixed = TRUE
  )
  
  expect_error(
    add_parameter(data.frame(), name = "a") |>
      add_parameter(name = "a"),
    "!name %in% x$name is not TRUE",
    fixed = TRUE
  )
  
  expect_error(
    add_parameter(data.frame(), name = "a", min = "a"),
    "is.numeric(min) is not TRUE",
    fixed = TRUE
  )
  
  expect_error(
    add_parameter(data.frame(), name = "a", max = "a"),
    "is.numeric(max) is not TRUE",
    fixed = TRUE
  )
  
  expect_error(
    add_parameter(data.frame(), name = "a", min = 1, max = 0),
    "max >= min is not TRUE",
    fixed = TRUE
  )
  
  expect_error(
    add_parameter(data.frame(), name = "a", initial_values = "a"),
    "is.numeric(initial_values) is not TRUE",
    fixed = TRUE
  )
  
  expect_error(
    add_parameter(data.frame(), name = "a", min = 0, initial_values = -1),
    "all(initial_values >= min) is not TRUE",
    fixed = TRUE
  )
  
  expect_error(
    add_parameter(data.frame(), name = "a", max = 1, initial_values = 2),
    "all(initial_values <= max) is not TRUE",
    fixed = TRUE
  )
  
  expect_error(
    add_parameter(data.frame(), name = "a", blocks = "a"),
    "is.numeric(blocks) is not TRUE",
    fixed = TRUE
  )
  
  p <- data.frame() |>
    add_parameter(name = "a", min = 0, max = 1, initial_values = 0.5)

  expect_equal(p$name, "a")
  expect_equal(p$min, 0)
  expect_equal(p$max, 1)
  expect_equal(unlist(p$initial_values), 0.5)
  expect_equal(unlist(p$blocks), 1)
  expect_equal(p$transform_type, 3)
  expect_equal(p$infer_parameter, TRUE)
  
  p <- p |>
    add_parameter(name = "b", min = 10, max = 10, initial_values = 10)
  
  expect_equal(p$name, c("a", "b"))
  expect_equal(p$min, c(0, 10))
  expect_equal(p$max, c(1, 10))
  expect_equal(unlist(p$initial_values), c(0.5, 10))
  expect_equal(unlist(p$blocks), c(1, 1))
  expect_equal(p$transform_type, c(3, 3))
  expect_equal(p$infer_parameter, c(TRUE, FALSE))
  
})
