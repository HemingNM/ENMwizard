# library(testthat)
# library(ENMwizard)

testthat::test_check("ENMwizard")


# devtools::use_testthat()

# devtools::create("ENMwizard")
# library(devtools)
# current.code <- devtools::as.package("/Volumes/Samsung SSD/neanderh/Documents/Dropbox/Artigos Projetos Trabalhos/1e. Pacotes R/ENMwizard")
# load_all(current.code)
# document(current.code)

## Wickham 2015 R packages (pg 98) Documenting package
# library(roxygen2)
# roxygen2::roxygenise()
devtools::document()
devtools::check()
devtools::build_win()


?load_all()
devtools::load_all()


Bvarieg.occ <- read.table(paste(system.file(package="dismo"), "/ex/bradypus.csv", sep=""), header=TRUE, sep=",")
colnames(Bvarieg.occ) <- c("SPEC", "LONG", "LAT")
spp.occ.list <- list(Bvarieg = Bvarieg.occ)
occ_polys <- f.poly.batch(spp.occ.list, o.path="occ_poly")
