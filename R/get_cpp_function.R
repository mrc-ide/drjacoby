#' Get .cpp function
#' 
#' Scrapes an appropriately tagged .cpp function for use in drJacoby. This is convenient as it
#' facilitates development and testing of a .cpp likelihood or prior via standard Rcpp::sourceCpp.
#' Once ready the same .cpp file can be imported for use in drJacoby using \code{cpp_function_get()}.
#' To allow this the user must add two tags. The first, \code{// drjacoby_function_start}, should
#' be on the line immediately before the start of the function and the second, \code{// drjacoby_function_end},
#' on the line immediately after the function. Only one function should be tagged in this manner per file.
#'
#' @param address Address of .cpp file
#'
#' @return A drJacoby compatible string function
#' @export
cpp_function_get <- function(address){
  assert_string(address)
  if(!grepl(".cpp", address)){
    stop("The file at address should be a .cpp file")
  }
  
  # Read in file line by line
  function_text <- readr::read_lines(address)
  # Isolate where the function starts
  f_start <- which(grepl("// drjacoby_function_start", function_text)) + 1
  if(length(f_start) == 0){
    stop("The string: '// drjacoby_function_start' could not be found in the specified file")
  }
  if(length(f_start) > 1){
    stop("The string: '// drjacoby_function_start' was found more than once in the specified file")
  }
  
  # Isolate where function ends
  f_end <- which(grepl("// drjacoby_function_end", function_text))
  if(length(f_end) == 0){
    stop("The string: '// drjacoby_function_end' could not be found in the specified file. This
         must be included on the final line to flag the function end")
  }
  
  # Create string of function
  paste(function_text[f_start:(f_end - 1)], collapse =  "\n")
}
