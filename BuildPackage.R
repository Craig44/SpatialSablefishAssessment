#'
#'
#' Build the sablefish spatial stock assessment model.
#'
#'

library(roxygen2)
library(devtools)
library(testthat)

remove.packages("SpatialSablefishAssessment")

## build package
pkgbuild::compile_dll() # need to compile src first
devtools::document()
#devtools::build()

## Run unit tests for Package
devtools::test()

## build Gitbook
bookdown::render_book(input = "Gitbook")
