add_parameter <- function(x, name, min = -Inf, max = Inf, initial_values = NULL, blocks = 1L){
  # Input checks
  stopifnot(is.data.frame(x))
  stopifnot(is.character(name))
  stopifnot(!name %in% x$name)
  stopifnot(is.numeric(min))
  stopifnot(is.numeric(max))
  stopifnot(max >= min)
  if(!is.null(initial_values)){
    stopifnot(is.numeric(initial_values))
    stopifnot(all(initial_values >= min))
    stopifnot(all(initial_values <= max))
  }
  stopifnot(is.numeric(blocks))
  blocks <- as.integer(blocks)
  
  # Append new parameter
  x |>
    rbind(
      data.frame(
        name = name,
        min = min,
        max = max,
        initial_values = I(list(initial_values)),
        blocks = I(list(blocks)),
        transform_type = get_transform_type(min, max),
        infer_parameter = !min == max
      )
    )
}

get_transform_type <- function(theta_min, theta_max){
  as.integer(2 * is.finite(theta_min) + is.finite(theta_max))
}

draw_initial_values <- function(min, max, transform_type, chains){
    p <- stats::runif(chains)
    if (transform_type == 0) {
      initial_values <- log(p) - log(1 - p)
    } else if (transform_type == 1) {
      initial_values <- log(p) + max
    } else if (transform_type == 2) {
      initial_values <- min - log(p)
    } else if (transform_type == 3) {
      initial_values <- min + (max - min) * p
    }
    initial_values <- list(initial_values)
  return(initial_values)
}

add_initial_values <- function(df_params, chains){
  for(i in 1:nrow(df_params)){
    if(is.null(df_params$initial_values[[i]])){
      df_params$initial_values[i] <- draw_initial_values(
        min = df_params$min[i],
        max = df_params$max[i],
        transform_type = df_params$transform_type[i],
        chains = chains
      )
    } else {
      if(length(df_params$initial_values[[i]]) != chains){
        stop("Length of user provided initial values must match number of chains")
      }
    }
  }
  return(df_params)
}

initial_theta <- function(df_params, chains, rungs){
  # Convert list of vectors into list of matrices, with nrow = n rungs
  initial_values <- list()
  for(i in 1:chains){
    chain_init <- sapply(df_params$initial_values, "[", i)
    initial_values[[i]]<- matrix(chain_init, ncol = length(chain_init), nrow = rungs, byrow = TRUE)
  }
  return(initial_values)
}
