#' @importFrom glue glue
.onLoad <- function(libname, pkgname) {
  path <- find.package('GPfates')
  if(!dir.exists(glue::glue("{path}/venv"))) {
    reinstall()
  }
}

#' @importFrom glue glue
reinstall <- function() {
  path <- find.package('GPfates')
  system(glue::glue("bash {path}/make {path}/venv"))
}
