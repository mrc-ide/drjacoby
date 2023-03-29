#------------------------------------------------
#' @title Define parameters dataframe
#'
#' @description Provides a convenient way of defining parameters in the format
#'   required by \code{run_mcmc()}. Each parameter must have the following three
#'   elements, defined in order:
#'   \itemize{
#'     \item \code{name} - the parameter name.
#'     \item \code{min} - the minimum value of the parameter. \code{-Inf} is
#'     allowed.
#'     \item \code{max} - the maximum value of the parameter. \code{Inf} is
#'     allowed.
#'   }
#'   There following arguments are also optional:
#'   \itemize{
#'     \item \code{init} - the initial value of the parameter. If running
#'     multiple chains a vector of initial values can be used to specify distinct
#'     values for each chain.
#'     \item \code{block} - which likelihood block(s) this parameter belongs to.
#'     See vignettes for instructions on using likelihood blocks.
#'   }
#'
#' @param ... a series of named input arguments.
#'   
#'
#' @export
#' @examples
#' define_params(name = "mu", min = -10, max = 10, init = 0,
#'               name = "sigma", min = 0, max = 5, init = c(1, 2))
#'               
#' define_params(name = "mu1", min = -10, max = 10, init = 0, block = 1,
#'               name = "mu2", min = -10, max = 10, init = 0, block = 2,
#'               name = "sigma", min = 0, max = 5, init = 1, block = c(1, 2))

define_params <- function(...) {
  x <- list(...)
  
  # check input format of arguments
  assert_gr(length(x), 0, message = "input cannot be empty")
  assert_in(names(x), c("name", "min", "max", "init", "block"))
  use_init <- ("init" %in% names(x))
  use_block <- ("block" %in% names(x))
  n_cols <- 3 + use_init + use_block
  if ((length(x) %% n_cols) != 0) {
    stop("must have the same number of inputs per parameter")
  }
  n_param <- length(x) / n_cols
  arg_names <-c("name", "min", "max")
  if (use_init) {
    arg_names <- c(arg_names, "init")
  }
  if (use_block) {
    arg_names <- c(arg_names, "block")
  }
  assert_eq(names(x), rep(arg_names, n_param))
  
  # create params dataframe
  v <- n_cols*(0:(n_param - 1))
  ret <- data.frame(name = unlist(x[1 + v]),
                    min = unlist(x[2 + v]),
                    max = unlist(x[3 + v]))
  if (use_init) {
    ret$init <- x[which(arg_names == "init") + v]
  }
  if (use_block) {
    ret$block <- x[which(arg_names == "block") + v]
  }
  
  # run checks and standardise format
  ret <- check_params(ret)
  
  return(ret)
}

#------------------------------------------------
# Check that params dataframe is formatted correctly, and return in standardised
# format (init and block coerced to list)
#' @noRd
check_params <- function(x) {
  
  # check dataframe has correct elements
  assert_dataframe(x)
  assert_in(c("name", "min", "max"), names(x),
            message = "df_params must contain the columns 'name', 'min', 'max'")
  if (any(duplicated(x$name))) {
    stop("parameter names must be unique")
  }
  use_init <- ("init" %in% names(x))
  use_block <- ("block" %in% names(x))
  
  # coerce init and block to list
  if (use_init) {
    if (!is.list(x$init)) {
      x$init <- as.list(x$init)
    }
  }
  if (use_block) {
    if (!is.list(x$block)) {
      x$block <- as.list(x$block)
    }
  }
  
  # check each row in turn
  for (i in seq_len(nrow(x))) {
    
    # check format
    assert_single_string(x$name[i], message = "parameter names must be character strings")
    assert_single_numeric(x$min[i], message = "min values must be single values")
    assert_single_numeric(x$max[i], message = "min values must be single values")
    if (use_init) {
      assert_vector_numeric(x$init[[i]], message = "init values must be numeric")
    }
    if (use_block) {
      assert_vector_numeric(x$block[[i]], message = "block values must be numeric")
    }
    
    # check order
    assert_leq(x$min[i], x$max[i], message = "min values must be less than or equal to max values")
    if (use_init) {
      assert_greq(x$init[[i]], x$min[i], message = "init values must be greater than or equal to min values")
      assert_leq(x$init[[i]], x$max[i], message = "init values must be less than or equal to max values")
    }
  }
  
  return(x)
}
