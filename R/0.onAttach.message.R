.onAttach <- function(libname = find.package("ENMwizard"), pkgname = "ENMwizard") {
  packageStartupMessage(("    Welcome to ENMwizard \n  To cite this package use: citation('ENMwizard'). \n  Please, remember to cite other packages that ENMwizard depends on: 'raster', 'dismo', 'ENMeval', 'spThin'"))
}
