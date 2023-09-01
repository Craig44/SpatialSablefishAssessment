#'
#' Build the sablefish spatial stock assessment model.
#'
#'

# delete previous package from lower right 'packages', click red x, then 'session' tab 'Restart R'...then run rest of script

library(roxygen2)
library(devtools)
library(testthat)

delete_previous_version = F ## sometimes you don't want to delete the previously installed version
install_using_Rstudio = T ## I prefer to use R-studio inbuilt panel to install.
## the r-stuido install is located i nthe 'Build' tab which should be next to the 'environment' & 'history' tab.

if(delete_previous_version)
  remove.packages("SpatialSablefishAssessment")
## sometimes if the package is in use, it won't delete. In this case sometimes
## I close all R environments and delete the package folder manually.
## You can find the location of the package directory using .libPaths()


## build package

######################################################################################
## Note!! if you change source code (C++) and want to recompile
## go to src/ & src/TMB/ and delete all files with extension '.dll' and '.o',
## otherwise this compile function wont register a change in source code.
######################################################################################

pkgbuild::compile_dll() # need to compile src first
devtools::document()

if(!install_using_Rstudio) {
  devtools::install() ## but I prefer to use R-studios "Build/Install/Button"
} else {
  cat("-----------\n")
  cat("This is a prompt to install using Rstudio's inbuild methods\n")
  cat("-----------\n")
  browser()
}

#########################################################################################################################
## If Install_Using_R_Studio==TRUE, you opt for the R-studios method and you need to do this manually at this point
## the code will pause and wait for you to do the install manually
## go to 'build' tab in upper right (next to 'environment'), click install button and everything should run and install the package
##########################################################################################################################


## Run unit tests for Package to check everything works as expected.
## make sure nothing is broken
## these won't prevent package installation/loading, but tell you if any new code is inconsistent with previous code (i.e., that something has changed)
## can add unit tests as update code by adding them to the 'tests' 'testthat' folder

devtools::test(stop_on_failure = F)

## Rebuild Gitbook to check you haven't broken anything
bookdown::render_book(input = 'Gitbook')
