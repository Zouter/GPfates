#' @importFrom glue glue
.onLoad <- function(libname, pkgname) {
  if(!dir.exists(glue::glue("bash {find.package('GPfates')}/venv"))) {
    reinstall()
  }
}

#' @export
#' @importFrom glue glue
reinstall <- function() {
  system(glue::glue("bash {find.package('GPfates')}/make {find.package('GPfates')}/venv"))
}
