.onLoad <- function(libname, pkgname) {
  op <- options()
  op.Trans2Kegg <- list(
    Trans2Kegg.path = "~/R-dev",
    Trans2Kegg.install.args = "",
    Trans2Kegg.name = "Chuck Roesel",
    Trans2Kegg.desc.author = "Roesel Chuck <croesel@gmail.com> [aut, cre]",
    Trans2Kegg.desc.license = "What license is it under?",
    Trans2Kegg.desc.suggests = NULL,
    Trans2Kegg.desc = list()
  )
  toset <- !(names(op.Trans2Kegg) %in% names(op))
  if(any(toset)) options(op.Trans2Kegg[toset])

  invisible()
}

.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Welcome to my package")
}
