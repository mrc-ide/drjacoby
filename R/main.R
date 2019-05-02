#' @useDynLib drjacoby, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

#------------------------------------------------
#' @title Dummy function
#'
#' @description Simple test function that demonstrates some of the features of
#'   this package.
#'
#' @details Takes a vector of values, returns the square.
#'
#' @param x vector of values
#'
#' @export
#' @examples
#' # Find square of first 100 values
#' dummy1(1:100)

dummy1 <- function(x = 1:5) {
  
  # print message to console
  message("running R dummy1 function")
  
  # get arguments in list form
  args <- list(x = x)
  
  # run C++ function with these arguments
  output_raw <- dummy1_cpp(args)
  
  # some optional processing of output
  message("processing output")
  ret <- output_raw$x_squared
  
  # return
  return(ret)
}

#------------------------------------------------
#' @title Run MCMC
#'
#' @description Run MCMC.
#'
#' @param df_params a dataframe of parameters.
#'
#' @export

run_mcmc <- function(df_params) {
  
  # check inputs
  
  # calculate transformation type for each parameter
  # 1 = [0,Inf] -> log
  # 2 = 
  df_params$trans_type <- mapply(function(x,y) {
    if (x == 0) {
      if (y == 1) {
        
      }
    }
  }, df_params$min, df_params$max)
  
  ret <- df_params
  
  # return
  return(ret)
}
