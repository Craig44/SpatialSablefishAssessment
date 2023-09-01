#'
#' Test the current assessment has a specific likelihood value
#'

#'
#' test-current-assessment-likelihood
#'
test_that("test-current-assessment-likelihood", {
  ## Read in mock data
  load(system.file("testdata", "MockAssessmentModel.RData",package="SpatialSablefishAssessment"))

  test_model <- TMB::MakeADFun(data = data,
                               parameters = parameters,
                               DLL = "SpatialSablefishAssessment_TMBExports", silent  = T)

  test_report = test_model$report()

  ## check negative loglikelihood
  expect_equal(sum(test_report$nll_weighted), 38163.63, tolerance = 0.001)
})

