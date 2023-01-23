

#' test-Projection
#' @description this is just to test the projection component of TagIntegrateModel currently just checks the model doesn't crash
#'
test_that("test-Projection", {
  ## Read in mock data
  load(system.file("testdata", "MockProjectionData.RData",package="SpatialSablefishAssessment"))
  data$model = "TagIntegratedValidate"
  ## no Z or movement
  data$apply_fixed_movement = 1 ## no movement initial age-structure should be only a function of ageing and M

  expect_no_condition(test_model <- TMB::MakeADFun(data = data,
                               parameters = parameters,
                               DLL = "SpatialSablefishAssessment_TMBExports", silent  = T))

  test_report = test_model$report()
})

