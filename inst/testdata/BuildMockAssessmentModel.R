#'
#' Build MockAssessment model. This initially loads "SimdataADMB.RDS" and  "SimparsADMB.RDS"
#' These were created in the "Spatial_Sablefish_Dev" repo, specficially in the "R/ADMBComparison/SetUpAssessmentModelUsingRpackage"
#' This script will add additional switches and inputs as the model changes from this point.
#'
set.seed(234)

data = readRDS(file = file.path("inst", "testdata", "SimdataADMB.RDS"))
parameters = readRDS(file = file.path("inst", "testdata", "SimparsADMB.RDS"))
region_key = data.frame(area = "Alaska", TMB_ndx = 0)
## change these as the model changes


## save them, these will be called by unit-tests to validate the model
## hasn't changed
save(data, parameters, region_key, file = file.path("inst", "testdata", "MockAssessmentModel.RData"))
