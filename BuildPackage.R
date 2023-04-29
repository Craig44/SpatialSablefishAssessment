#'
#' Build the sablefish spatial stock assessment model.
#'
#'

library(roxygen2)
library(devtools)
library(testthat)

#remove.packages("SpatialSablefishAssessment")

## build package
## Note!! if you change source code (C++) and want to recompile
## go to src/ & src/TMB/ and delete all files with extension '.dll' and '.o',
## otherwise this compile function wont register a change in source code.
pkgbuild::compile_dll() # need to compile src first
devtools::document()
#devtools::check()
#rcmdcheck::rcmdcheck(args = "--no-manual", error_on = "error")
## Run unit tests for Package
result = devtools::test(stop_on_failure = F)

devtools::build()

test_active_file(file.path("tests","testthat","test-validate-and-production-compatible.R"))
## build Gitbook
bookdown::render_book(input = 'Gitbook')
