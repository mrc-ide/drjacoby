
#### HELPER FUNCTIONS ####################################################################

#------------------------------------------------
# for single value, return value as string. For vector of values return string
# of comma-separated values enclosed in curly brackets
#' @noRd
nice_format <- function(x) {
  if (is.null(x)) {
    return("")
  }
  if (length(x)==1) {
    ret <- as.character(x)
  } else {
    ret <- paste0("{", paste(x, collapse = ", "), "}")
  }
  return(ret)
}

#### BASIC OBJECT TYPES ####################################################################
# messages should be NULL with the default speficied within function

#------------------------------------------------
# x is NULL
#' @noRd
assert_null <- function(x, message = NULL,
                        name = paste(deparse(substitute(x)), collapse = "")) {
  # default message
  if (is.null(message)) {
    message <- sprintf("%s must be null", name)
  }
  
  if (!is.null(x)) {
    stop(message, call. = FALSE)
  }
  return(TRUE)
}

#------------------------------------------------
# x is not NULL
#' @noRd
assert_non_null <- function(x, message = NULL,
                            name = paste(deparse(substitute(x)), collapse = "")) {
  # default message
  if (is.null(message)) {
    message <- sprintf("%s cannot be null", name)
  }

  if (is.null(x)) {
    stop(message, call. = FALSE)
  }
  return(TRUE)
}

#------------------------------------------------
# x is NA
#' @noRd
assert_NA <- function(x, message = NULL,
                      name = paste(deparse(substitute(x)), collapse = "")) {
  # default message
  if (is.null(message)) {
    message <- sprintf("%s must be NA", name)
  }
  
  assert_non_null(x)
  if (!all(is.na(x))) {
    stop(message, call. = FALSE)
  }
  return(TRUE)
}

#------------------------------------------------
# x is not NA
#' @noRd
assert_non_NA <- function(x, message = NULL,
                      name = paste(deparse(substitute(x)), collapse = "")) {
  # default message
  if (is.null(message)) {
    message <- sprintf("%s cannot be NA", name)
  }
  
  if (any(is.na(x))) {
    stop(message, call. = FALSE)
  }
  return(TRUE)
}

#------------------------------------------------
# x is atomic
#' @noRd
assert_atomic <- function(x, message = NULL,
                          name = paste(deparse(substitute(x)), collapse = "")) {
  # default message
  if (is.null(message)) {
    message <- sprintf("%s must be atomic (see ?is.atomic)", name)
  }
  
  if (!is.atomic(x)) {
    stop(message, call. = FALSE)
  }
  return(TRUE)
}

#------------------------------------------------
# x is character string
#' @noRd
assert_string <- function(x, message = NULL,
                          name = paste(deparse(substitute(x)), collapse = "")) {
  # default message
  if (is.null(message)) {
    message <- sprintf("%s must be a character string", name)
  }

  if (!is.character(x)) {
    stop(message, call. = FALSE)
  }
  return(TRUE)
}

#------------------------------------------------
# x is logical
#' @noRd
assert_logical <- function(x, message = NULL,
                           name = paste(deparse(substitute(x)), collapse = "")) {
  # default message
  if (is.null(message)) {
    message <- sprintf("%s must be logical", name)
  }

  if (!is.logical(x)) {
    stop(message, call. = FALSE)
  }
  return(TRUE)
}

#------------------------------------------------
# x is numeric
#' @noRd
assert_numeric <- function(x, message = NULL,
                           name = paste(deparse(substitute(x)), collapse = "")) {
  # default message
  if (is.null(message)) {
    message <- sprintf("%s must be numeric", name)
  }

  if (!is.numeric(x)) {
    stop(message, call. = FALSE)
  }
  return(TRUE)
}

#------------------------------------------------
# x is integer
#' @noRd
assert_int <- function(x, message = NULL,
                       name = paste(deparse(substitute(x)), collapse = "")) {
  # default message
  if (is.null(message)) {
    message <- sprintf("%s must be integer valued", name)
  }

  assert_numeric(x, name = name)
  if (!isTRUE(all.equal(x, as.integer(x), check.attributes = FALSE))) {
    stop(message, call. = FALSE)
  }
  return(TRUE)
}

#------------------------------------------------
# x is positive (with or without zero allowed)
#' @noRd
assert_pos <- function(x, zero_allowed = TRUE, message = NULL,
                       name = paste(deparse(substitute(x)), collapse = "")) {
  # default message
  if (is.null(message)) {
    if (zero_allowed) {
      message <- sprintf("%s must be greater than or equal to zero", name)
    } else {
      message <- sprintf("%s must be greater than zero", name)
    }
  }

  assert_numeric(x, name = name)
  if (zero_allowed) {
    if (!all(x >= 0)) {
      stop(message, call. = FALSE)
    }
  } else {
    if (!all(x > 0)) {
      stop(message, call. = FALSE)
    }
  }
  return(TRUE)
}

#------------------------------------------------
# x is a vector (and is not a list or another recursive type)
#' @noRd
assert_vector <- function(x, message = NULL,
                          name = paste(deparse(substitute(x)), collapse = "")) {
  # default message
  if (is.null(message)) {
    message <- sprintf("%s must be a non-recursive vector", name)
  }

  if (!is.vector(x) || is.recursive(x)) {
    stop(message, call. = FALSE)
  }
  return(TRUE)
}

#------------------------------------------------
# x is a matrix
#' @noRd
assert_matrix <- function(x, message = NULL,
                          name = paste(deparse(substitute(x)), collapse = "")) {
  # default message
  if (is.null(message)) {
    message <- sprintf("%s must be a matrix", name)
  }

  if (!is.matrix(x)) {
    stop(message, call. = FALSE)
  }
  return(TRUE)
}

#------------------------------------------------
# x is a list
assert_list <- function(x, message = NULL,
                        name = paste(deparse(substitute(x)), collapse = "")) {
  # default message
  if (is.null(message)) {
    message <- sprintf("%s must be a list", name)
  }

  if (!is.list(x)) {
    stop(message, call. = FALSE)
  }
  return(TRUE)
}

#------------------------------------------------
# x is named
assert_named <- function(x, message = NULL,
                         name = paste(deparse(substitute(x)), collapse = "")) {
  # default message
  if (is.null(message)) {
    message <- sprintf("%s must be named", name)
  }

  if (is.null(names(x))) {
    stop(message, call. = FALSE)
  }
  return(TRUE)
}

#------------------------------------------------
# x is a data frame
assert_dataframe <- function(x, message = NULL,
                             name = paste(deparse(substitute(x)), collapse = "")) {
  # default message
  if (is.null(message)) {
    message <- sprintf("%s must be a data frame", name)
  }
  
  if (!is.data.frame(x)) {
    stop(message, call. = FALSE)
  }
  return(TRUE)
}

#------------------------------------------------
# x inherits from class c
assert_class <- function(x, c, message = NULL,
                         name = paste(deparse(substitute(x)), collapse = "")) {
  # default message
  if (is.null(message)) {
    message <- sprintf("%s must inherit from class '%s'", name, c)
  }
  
  if (!inherits(x, c)) {
    stop(message, call. = FALSE)
  }
  return(TRUE)
}


#### COMPOUND OBJECT TYPES ####################################################################
# messages should be NULL, allowing them to be propagated through to more
# elementary levels

#------------------------------------------------
# x is atomic and single valued (has length 1)
#' @noRd
assert_single <- function(x, message = NULL,
                          name = paste(deparse(substitute(x)), collapse = "")) {
  
  # default message
  if (is.null(message)) {
    message <- sprintf("%s must be a single value (atomic and of length 1)", name)
  }
  
  assert_non_null(x, name = name, message = message)
  assert_atomic(x, name = name, message = message)
  assert_length(x, 1, name = name, message = message)
  return(TRUE)
}

#------------------------------------------------
# x is single character string
#' @noRd
assert_single_string <- function(x, message = NULL,
                                 name = paste(deparse(substitute(x)), collapse = "")) {
  
  # default message
  if (is.null(message)) {
    message <- sprintf("%s must be a single character string", name)
  }
  
  assert_length(x, n = 1, name = name, message = message)
  assert_string(x, name = name, message = message)
  return(TRUE)
}

#------------------------------------------------
# x is single logical
#' @noRd
assert_single_logical <- function(x, message = NULL,
                                  name = paste(deparse(substitute(x)), collapse = "")) {
  
  # default message
  if (is.null(message)) {
    message <- sprintf("%s must be a single logical value", name)
  }
  
  assert_length(x, n = 1, name = name, message = message)
  assert_logical(x, name = name, message = message)
  return(TRUE)
}

#------------------------------------------------
# x is single numeric
#' @noRd
assert_single_numeric <- function(x, message = NULL,
                                  name = paste(deparse(substitute(x)), collapse = "")) {
  
  # default message
  if (is.null(message)) {
    message <- sprintf("%s must be a single numeric value", name)
  }
  
  assert_length(x, n = 1, name = name, message = message)
  assert_numeric(x, name = name, message = message)
  return(TRUE)
}

#------------------------------------------------
# x is single integer
#' @noRd
assert_single_int <- function(x, message = NULL,
                              name = paste(deparse(substitute(x)), collapse = "")) {
  
  # default message
  if (is.null(message)) {
    message <- sprintf("%s must be a single integer", name)
  }
  
  assert_length(x, n = 1, name = name, message = message)
  assert_int(x, name = name, message = message)
  return(TRUE)
}

#------------------------------------------------
# x is single positive (with or without zero allowed)
#' @noRd
assert_single_pos <- function(x, zero_allowed = TRUE, message = NULL,
                              name = paste(deparse(substitute(x)), collapse = "")) {
  
  # default message
  if (is.null(message)) {
    if (zero_allowed) {
      message <- sprintf("%s must be a single positive value or zero", name)
    } else {
      message <- sprintf("%s must be a single positive value", name)
    }
  }
  
  assert_length(x, n = 1, name = name, message = message)
  assert_pos(x, zero_allowed = zero_allowed, name = name, message = message)
  return(TRUE)
}

#------------------------------------------------
# x is positive integer (with or without zero allowed)
#' @noRd
assert_pos_int <- function(x, zero_allowed = TRUE, message = NULL,
                           name = paste(deparse(substitute(x)), collapse = "")) {
  
  # default message
  if (is.null(message)) {
    if (zero_allowed) {
      message <- sprintf("%s must be a positive integer or zero", name)
    } else {
      message <- sprintf("%s must be a positive integer", name)
    }
  }
  
  assert_int(x, name = name, message = message)
  assert_pos(x, zero_allowed = zero_allowed, name = name, message = message)
  return(TRUE)
}

#------------------------------------------------
# x is single positive integer (with or without zero allowed)
#' @noRd
assert_single_pos_int <- function(x, zero_allowed = TRUE, message = NULL,
                                  name = paste(deparse(substitute(x)), collapse = "")) {
  
  # default message
  if (is.null(message)) {
    if (zero_allowed) {
      message <- sprintf("%s must be a single positive integer or zero", name)
    } else {
      message <- sprintf("%s must be a single positive integer", name)
    }
  }
  
  assert_length(x, n = 1, name = name, message = message)
  assert_pos_int(x, zero_allowed = zero_allowed, name = name, message = message)
  return(TRUE)
}

#------------------------------------------------
# x is single value bounded between limits
#' @noRd
assert_single_bounded <- function(x, left = 0, right = 1, inclusive_left = TRUE, inclusive_right = TRUE,
                                  message = NULL,
                                  name = paste(deparse(substitute(x)), collapse = "")) {
  
  # default message
  if (is.null(message)) {
    message <- sprintf("%s must be a single value bounded between %s and %s", name, left, right)
  }
  
  assert_length(x, n = 1, name = name, message = message)
  assert_bounded(x, left = left, right = right,
                 inclusive_left = inclusive_left, inclusive_right = inclusive_right,
                 message = message,
                 name = name)
  return(TRUE)
}

#------------------------------------------------
# x is a numeric vector
#' @noRd
assert_vector_numeric <- function(x, message = NULL,
                                  name = paste(deparse(substitute(x)), collapse = "")) {
  
  # default message
  if (is.null(message)) {
    message <- sprintf("%s must be a numeric vector", name)
  }
  
  assert_numeric(x, name = name, message = message)
  assert_vector(x, name = name, message = message)
  return(TRUE)
}

#------------------------------------------------
# x is a vector of integers
#' @noRd
assert_vector_int <- function(x, message = NULL,
                              name = paste(deparse(substitute(x)), collapse = "")) {
  
  # default message
  if (is.null(message)) {
    message <- sprintf("%s must be a vector of integers", name)
  }
  
  assert_int(x, name = name, message = message)
  assert_vector(x, name = name, message = message)
  return(TRUE)
}

#------------------------------------------------
# x is a vector of positive values
#' @noRd
assert_vector_pos <- function(x, zero_allowed = TRUE,
                              message = NULL,
                              name = paste(deparse(substitute(x)), collapse = "")) {
  
  # default message
  if (is.null(message)) {
    if (zero_allowed) {
      message <- sprintf("%s must be a vector of positive values or zero", name)
    } else {
      message <- sprintf("%s must be a vector of positive values", name)
    }
  }
  
  assert_pos(x, zero_allowed = zero_allowed, name = name, message = message)
  assert_vector(x, name = name, message = message)
  return(TRUE)
}

#------------------------------------------------
# x is a vector of positive integers
#' @noRd
assert_vector_pos_int <- function(x, zero_allowed = TRUE,
                                  message = NULL,
                                  name = paste(deparse(substitute(x)), collapse = "")) {
  
  # default message
  if (is.null(message)) {
    if (zero_allowed) {
      message <- sprintf("%s must be a vector of positive integers or zero", name)
    } else {
      message <- sprintf("%s must be a vector of positive integers", name)
    }
  }
  
  assert_pos(x, zero_allowed = zero_allowed, name = name, message = message)
  assert_vector(x, name = name, message = message)
  assert_int(x, name = name, message = message)
  return(TRUE)
}

#------------------------------------------------
# x is a vector of bounded values
#' @noRd
assert_vector_bounded <- function(x, left = 0, right = 1, inclusive_left = TRUE, inclusive_right = TRUE,
                                  message = NULL,
                                  name = paste(deparse(substitute(x)), collapse = "")) {
  
  # default message
  if (is.null(message)) {
    message <- sprintf("%s must be a vector of values bounded between %s and %s", name, left, right)
  }
  
  assert_vector(x, name = name, message = message)
  assert_bounded(x, left = left, right = right,
                 inclusive_left = inclusive_left, inclusive_right = inclusive_right,
                 message = message,
                 name = name)
  return(TRUE)
}

#------------------------------------------------
# x is a vector of strings
#' @noRd
assert_vector_string <- function(x, message = NULL,
                                 name = paste(deparse(substitute(x)), collapse = "")) {
  
  # default message
  if (is.null(message)) {
    message <- sprintf("%s must be a vector of character strings", name)
  }
  
  assert_vector(x, name = name, message = message)
  assert_string(x, name = name, message = message)
  return(TRUE)
}

#------------------------------------------------
# x is a matrix of numeric values
#' @noRd
assert_matrix_numeric <- function(x, message = NULL,
                                  name = paste(deparse(substitute(x)), collapse = "")) {
  
  # default message
  if (is.null(message)) {
    message <- sprintf("%s must be a matrix of numeric values", name)
  }
  
  assert_matrix(x, name = name, message = message)
  assert_numeric(x, name = name, message = message)
  return(TRUE)
}

#------------------------------------------------
# x is a named list
assert_list_named <- function(x, message = NULL,
                              name = paste(deparse(substitute(x)), collapse = "")) {
  
  # default message
  if (is.null(message)) {
    message <- sprintf("%s must be a named list", name)
  }
  
  assert_named(x, name = name, message = message)
  assert_list(x, name = name, message = message)
  return(TRUE)
}

#------------------------------------------------
# x is a plotting limit, i.e. contains two increasing values
assert_limit <- function(x, message = NULL,
                         name = paste(deparse(substitute(x)), collapse = "")) {
  
  # default message
  if (is.null(message)) {
    message <- sprintf("%s must be vector of two increasing values", name)
  }
  
  assert_vector(x, name = name, message = message)
  assert_length(x, 2, name = name, message = message)
  assert_numeric(x, name = name, message = message)
  assert_increasing(x, name = name, message = message)
  return(TRUE)
}


#### VALUE COMPARISONS ####################################################################

#------------------------------------------------
# x and y are equal in all matched comparisons. x and y can be any type
#' @noRd
assert_eq <- function(x, y, message = NULL,
                      name_x = paste(deparse(substitute(x)), collapse = ""), name_y = nice_format(y)) {
  # default message
  if (is.null(message)) {
    message <- sprintf("%s must equal %s", name_x, name_y)
  }
  
  assert_non_null(x, name = name_x, message = message)
  assert_non_null(y, name = name_y, message = message)
  assert_same_length(x, y, name_x = name_x, name_y = name_y, message = message)
  if (!isTRUE(all.equal(x, y, check.attributes = FALSE))) {
    stop(message, call. = FALSE)
  }
  return(TRUE)
}

#------------------------------------------------
# x and y are unequal in all matched comparisons. x and y can be any type
#' @noRd
assert_neq <- function(x, y, message = NULL,
                       name_x = paste(deparse(substitute(x)), collapse = ""), name_y = nice_format(y)) {
  # default message
  if (is.null(message)) {
    message <- sprintf("%s cannot equal %s", name_x, name_y)
  }

  assert_non_null(x, name = name_x, message = message)
  assert_non_null(y, name = name_y, message = message)
  assert_same_length(x, y, name_x = name_x, name_y = name_y, message = message)
  if (any(x == y)) {
    stop(message, call. = FALSE)
  }
  return(TRUE)
}

#------------------------------------------------
# x is greater than y in all matched comparisons
#' @noRd
assert_gr <- function(x, y, message = NULL,
                      name_x = paste(deparse(substitute(x)), collapse = ""),
                      name_y = nice_format(y)) {
  # default message
  if (is.null(message)) {
    message <- sprintf("%s must be greater than %s", name_x, name_y)
  }

  assert_numeric(x, name = name_x, message = message)
  assert_numeric(y, name = name_y, message = message)
  assert_in(length(y), c(1, length(x)), message = message)
  if (!all(x > y)) {
    stop(message, call. = FALSE)
  }
  return(TRUE)
}

#------------------------------------------------
# x is greater than or equal to y in all matched comparisons
#' @noRd
assert_greq <- function(x, y, message = NULL,
                        name_x = paste(deparse(substitute(x)), collapse = ""),
                        name_y = nice_format(y)) {
  # default message
  if (is.null(message)) {
    message <- sprintf("%s must be greater than or equal to %s", name_x, name_y)
  }

  assert_numeric(x, name = name_x, message = message)
  assert_numeric(y, name = name_y, message = message)
  assert_in(length(y), c(1, length(x)), message = message)
  if (!all(x >= y)) {
    stop(message, call. = FALSE)
  }
  return(TRUE)
}

#------------------------------------------------
# x is less than y in all matched comparisons
#' @noRd
assert_le <- function(x, y, message = NULL,
                      name_x = paste(deparse(substitute(x)), collapse = ""),
                      name_y = nice_format(y)) {
  # default message
  if (is.null(message)) {
    message <- sprintf("%s must be less than %s", name_x, name_y)
  }

  assert_numeric(x, name = name_x, message = message)
  assert_numeric(y, name = name_y, message = message)
  assert_in(length(y), c(1, length(x)), message = message)
  if (!all(x < y)) {
    stop(message, call. = FALSE)
  }
  return(TRUE)
}

#------------------------------------------------
# x is less than or equal to y in all matched comparisons
#' @noRd
assert_leq <- function(x, y, message = NULL,
                       name_x = paste(deparse(substitute(x)), collapse = ""),
                       name_y = nice_format(y)) {
  # default message
  if (is.null(message)) {
    message <- sprintf("%s must be less than or equal to %s", name_x, name_y)
  }

  assert_numeric(x, name = name_x, message = message)
  assert_numeric(y, name = name_y, message = message)
  assert_in(length(y), c(1, length(x)), message = message)
  if (!all(x <= y)) {
    stop(message, call. = FALSE)
  }
  return(TRUE)
}

#------------------------------------------------
# x is between bounds (inclusive or exclusive)
#' @noRd
assert_bounded <- function(x, left = 0, right = 1, inclusive_left = TRUE, inclusive_right = TRUE,
                           message = NULL,
                           name = paste(deparse(substitute(x)), collapse = "")) {

  if (inclusive_left) {
    assert_greq(x, left, name_x = name, message = message)
  } else {
    assert_gr(x, left, name_x = name, message = message)
  }
  if (inclusive_right) {
    assert_leq(x, right, name_x = name, message = message)
  } else {
    assert_le(x, right, name_x = name, message = message)
  }
  return(TRUE)
}

#------------------------------------------------
# all x are in y
#' @noRd
assert_in <- function(x, y, message = NULL,
                      name_x = paste(deparse(substitute(x)), collapse = ""), name_y = nice_format(y)) {
  # default message
  if (is.null(message)) {
    message <- sprintf("all %s must be in %s", name_x, name_y)
  }

  assert_non_null(x, name = name_x, message = message)
  assert_non_null(y, name = name_y, message = message)
  if (!all(x %in% y)) {
    stop(message, call. = FALSE)
  }
  return(TRUE)
}

#------------------------------------------------
# none of x are in y
#' @noRd
assert_not_in <- function(x, y, message = NULL,
                          name_x = paste(deparse(substitute(x)), collapse = ""), name_y = nice_format(y)) {
  # default message
  if (is.null(message)) {
    message <- sprintf("none of %s can be in %s", name_x, name_y)
  }

  assert_non_null(x, name = name_x, message = message)
  assert_non_null(y, name = name_y, message = message)
  if (any(x %in% y)) {
    stop(message, call. = FALSE)
  }
  return(TRUE)
}


#### DIMENSIONS ####################################################################

#------------------------------------------------
# length(x) equals n
#' @noRd
assert_length <- function(x, n, message = NULL,
                          name = paste(deparse(substitute(x)), collapse = "")) {
  # default message
  if (is.null(message)) {
    message <- sprintf("%s must be of length %s", name, n)
  }

  assert_pos_int(n, name = name, message = message)
  if (length(x) != n[1]) {
    stop(message, call. = FALSE)
  }
  return(TRUE)
}

#------------------------------------------------
# x and y are same length
#' @noRd
assert_same_length <- function(x, y, message =  NULL,
                               name_x = paste(deparse(substitute(x)), collapse = ""),
                               name_y = paste(deparse(substitute(y)))) {
  # default message
  if (is.null(message)) {
    message <- sprintf("%s and %s must be the same length", name_x, name_y)
  }

  if (length(x) != length(y)) {
    stop(message, call. = FALSE)
  }
  return(TRUE)
}

#------------------------------------------------
# multiple objects all same length
#' @noRd
assert_same_length_multiple <- function(...) {
  l <- mapply(length, list(...))
  if (length(unique(l)) != 1) {
    l_names <- sapply(match.call(expand.dots = FALSE)$..., deparse)
    stop(sprintf("variables %s must be the same length", nice_format(l_names)), call. = FALSE)
  }
  return(TRUE)
}

#------------------------------------------------
# x is two-dimensional
#' @noRd
assert_2d <- function(x, message = NULL,
                      name = paste(deparse(substitute(x)), collapse = "")) {
  # default message
  if (is.null(message)) {
    message <- sprintf("%s must be two-dimensional", name)
  }

  is_2d <- FALSE
  if (!is.null(dim(x))) {
    is_2d <- (length(dim(x)) == 2)
  }
  if (!is_2d) {
    stop(message, call. = FALSE)
  }
  return(TRUE)
}

#------------------------------------------------
# nrow(x) equals n
#' @noRd
assert_nrow <- function(x, n, message = NULL,
                        name = paste(deparse(substitute(x)), collapse = "")) {
  # default message
  if (is.null(message)) {
    message <- sprintf("%s must have %s rows", name, n)
  }
  
  assert_2d(x, name = name, message = message)
  if (nrow(x) != n) {
    stop(message, call. = FALSE)
  }
  return(TRUE)
}

#------------------------------------------------
# ncol(x) equals n
#' @noRd
assert_ncol <- function(x, n, message = NULL,
                        name = paste(deparse(substitute(x)), collapse = "")) {
  # default message
  if (is.null(message)) {
    message <- sprintf("%s must have %s cols", name, n)
  }
  
  assert_2d(x, name = name, message = message)
  if (ncol(x) != n) {
    stop(message, call. = FALSE)
  }
  return(TRUE)
}

#------------------------------------------------
# dim(x) equals y
#' @noRd
assert_dim <- function(x, y, message = NULL,
                       name = paste(deparse(substitute(x)), collapse = "")) {
  # default message
  if (is.null(message)) {
    message <- sprintf("%s must have %s rows and %s columns", name, y[1], y[2])
  }
  
  assert_2d(x, name = name, message = message)
  assert_pos_int(y, name = name, message = message)
  assert_length(y, 2, name = name, message = message)
  if (nrow(x) != y[1] | ncol(x) != y[2]) {
    stop(message, call. = FALSE)
  }
  return(TRUE)
}

#------------------------------------------------
# x is square matrix
#' @noRd
assert_square_matrix <- function(x, message = NULL,
                                 name = paste(deparse(substitute(x)), collapse = "")) {
  # default message
  if (is.null(message)) {
    message <- sprintf("%s must be a square matrix", name)
  }
  
  assert_matrix(x, name = name, message = message)
  if (nrow(x) != ncol(x)) {
    stop(message, call. = FALSE)
  }
  return(TRUE)
}

#------------------------------------------------
# is symmetric matrix
#' @noRd
assert_symmetric_matrix <- function(x, message = NULL,
                                    name = paste(deparse(substitute(x)), collapse = "")) {
  # default message
  if (is.null(message)) {
    message <- sprintf("%s must be a symmetric matrix", name)
  }

  assert_square_matrix(x, name = name, message = message)
  if (!isSymmetric(x)) {
    stop(message, call. = FALSE)
  }
  return(TRUE)
}

#### MISC ####################################################################

#------------------------------------------------
# x contains no duplicates
#' @noRd
assert_noduplicates <- function(x, message = NULL,
                                name = paste(deparse(substitute(x)), collapse = "")) {
  # default message
  if (is.null(message)) {
    message <- sprintf("%s must contain no duplicates", name)
  }

  if (any(duplicated(x))) {
    stop(message, call. = FALSE)
  }
  return(TRUE)
}

#------------------------------------------------
# file exists at chosen path
#' @noRd
assert_file_exists <- function(x, message = NULL,
                               name = paste(deparse(substitute(x)), collapse = "")) {
  # default message
  if (is.null(message)) {
    message <- sprintf("%s must be a string defining a file that exists at that path", name)
  }
  
  assert_single_string(x, name = name, message = message)
  if (!file.exists(x)) {
    stop(message, call. = FALSE)
  }
  return(TRUE)
}

#------------------------------------------------
# x is increasing
#' @noRd
assert_increasing <- function(x, message = NULL,
                              name = paste(deparse(substitute(x)), collapse = "")) {
  # default message
  if (is.null(message)) {
    message <- sprintf("%s must be increasing", name)
  }
  
  assert_non_null(x, name = name, message = message)
  assert_numeric(x, name = name, message = message)
  if (!isTRUE(all.equal(x, sort(x)))) {
    stop(message, call. = FALSE)
  }
  return(TRUE)
}

#------------------------------------------------
# x is decreasing
#' @noRd
assert_decreasing <- function(x, message = NULL,
                              name = paste(deparse(substitute(x)), collapse = "")) {
  # default message
  if (is.null(message)) {
    message <- sprintf("%s must be decreasing", name)
  }
  
  assert_non_null(x, name = name, message = message)
  assert_numeric(x, name = name, message = message)
  if (!isTRUE(all.equal(x, sort(x, decreasing = TRUE)))) {
    stop(message, call. = FALSE)
  }
  return(TRUE)
}

