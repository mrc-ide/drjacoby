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
  stopifnot(length(x) > 0)
  if(!all(c("name", "min", "max") %in% names(x))){
    stop("Columns: name, min, max, must all be present in df_params")
  }
  
  
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
  stopifnot(identical(names(x), rep(arg_names, n_param)))
  
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
  stopifnot(is.data.frame(x))
  if(!all(c("name", "min", "max") %in% names(x))){
    stop("df_params must contain the columns 'name', 'min', 'max'")
  }
  if (any(duplicated(x$name))) {
    stop("parameter names must be unique")
  }
  use_init <- ("init" %in% names(x))
  use_block <- ("block" %in% names(x))
  
  # coerce init and block to list
  if(use_init) {
    if (!is.list(x$init)) {
      x$init <- as.list(x$init)
    }
  }
  if(use_block) {
    if (!is.list(x$block)) {
      x$block <- as.list(x$block)
    }
  }
  
  # check each row in turn
  for(i in seq_len(nrow(x))) {
    # check format
    stopifnot(is.character(x$name[i]))
    stopifnot(is.numeric(x$min[i]))
    stopifnot(is.numeric(x$max[i]))
    
    if(use_init) {
      stopifnot(is.numeric(x$init[[i]]))
    }
    if(use_block) {
      stopifnot(is.numeric(x$block[[i]]))
    }
    stopifnot(
      x$min[i] <= x$max[i]
    )
    stopifnot(
      all(x$init[[i]] >= x$min[i])
    )
    stopifnot(
      all(x$init[[i]] <= x$max[i])
    )
  }
  
  return(x)
}

get_transform_type <- function(theta_min, theta_max){
  as.integer(2 * is.finite(theta_min) + is.finite(theta_max))
}

# define default init values
set_init <- function(df_params, chains){
  use_init <- ("init" %in% names(df_params))
  if(use_init){
    check_init(df_params, chains)
  } else {
    df_params$init <- get_init(df_params, chains)
  }
  init <- lapply(1:chains, function(x){
    sapply(df_params$init, '[', x)
  }
  )
  return(init)
}

check_init <- function(df_params, chains){
  for (i in 1:nrow(df_params)) {
    stopifnot(length(df_params$init[[i]]) == chains)
  }
}

get_init <- function(df_params, chains){
  init_list <- list()
  for (i in 1:nrow(df_params)) {
    transform_type <- get_transform_type(df_params[i,]$min, df_params[i,]$max)
    p <- runif(chains)
    if (transform_type == 0) {
      init_list[[i]] <- log(p) - log(1 - p)
    } else if (transform_type == 1) {
      init_list[[i]] <- log(p) + df_params$max[i]
    } else if (transform_type == 2) {
      init_list[[i]] <- df_params$min[i] - log(p)
    } else if (transform_type == 3) {
      init_list[[i]] <- df_params$min[i] + (df_params$max[i] - df_params$min[i])*p
    }
  }
  return(init_list)
}

set_blocks <- function(df_params){
  if(!"block" %in% names(df_params)){
    df_params$block <- 1
  }
  blocks_list <- lapply(df_params$block, as.integer)
  return(blocks_list)
}
