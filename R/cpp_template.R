#' Create template for cpp
#'
#' @param save_as of file to create, relative to root of active project.
#'
#' @return 
#' @export
cpp_template <- function (save_as){
  usethis::use_template(
    template = "cpp_template.cpp",
    save_as = save_as,
    open = rlang::is_interactive(),
    package = "drjacoby"
  )
}


