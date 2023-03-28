# define default init values
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

check_init <- function(df_params, chains){
  for (i in 1:nrow(df_params)) {
    if (length(df_params$init[[i]]) != 1) {
      assert_length(df_params$init[[i]], chains, message = paste0("must define one df_params$init value per parameter, ",
                                                                  "or alternatively a list of values one for each chain"))
    }
  }
}
