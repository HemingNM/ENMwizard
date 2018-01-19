library(devtools)
# devtools::use_testthat()
# devtools::use_vignette("introduction")
# senha Git: GT!he36
## Wickham 2015 R packages (pg 98) Documenting package
# ?devtools::use_package()
# devtools::create()
devtools::use_build_ignore(c("documenting", "parallel.md", "readme.md"))
devtools::document()
devtools::check()
# devtools::load_all()
# roxygen2::roxygenise()

# for building windows package
devtools::build_win()




# # general function documentation
#
# #' Function Title (short description)
# #'
# #' General function description. A short paragraph (or more) describing what the function does.
# #' @param arg1 List of occurence data. See argument "occ_locs" in mxnt.cp.
# #' @inheritParams function1
# #' @return objects returned from function
# #' @examples
# #' plot(mxnt.mdls.preds.lst[[1]][[4]]) # MaxEnt predictions, based on the model selection criteria
# #' @export
#


# Bvarieg.occ <- read.table(paste(system.file(package="dismo"), "/ex/bradypus.csv", sep=""), header=TRUE, sep=",")
# colnames(Bvarieg.occ) <- c("SPEC", "LONG", "LAT")
# spp.occ.list <- list(Bvarieg = Bvarieg.occ)
# occ_polys <- f.poly.batch(spp.occ.list, o.path="occ_poly")
