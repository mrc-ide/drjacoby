#' Create template for cpp
#'
#' @param save_as of file to create, relative to root of active project.
#'
#' @export
cpp_template <- function (save_as){
  ext <- tools::file_ext(save_as)
  if(ext != ".cpp"){
    stop("File must be .cpp")
  }
  usethis::use_template(
    template = "cpp_template.cpp",
    save_as = save_as,
    open = rlang::is_interactive(),
    package = "drjacoby"
  )
}


