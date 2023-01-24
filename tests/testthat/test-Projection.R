

#' test-Projection
#' @description this is just to test the projection component of TagIntegrateModel currently just checks the model doesn't crash
#'
test_that("test-Projection", {
  ## Read in mock data
  load(system.file("testdata", "MockProjectionData.RData",package="SpatialSablefishAssessment"))
  data$model = "TagIntegratedValidate"
  ## no Z or movement
  data$apply_fixed_movement = 1 ## no movement initial age-structure should be only a function of ageing and M

  ## make sure it doesn't crash
  expect_no_condition(test_model <- TMB::MakeADFun(data = data,
                               parameters = parameters,
                               DLL = "SpatialSablefishAssessment_TMBExports", silent  = T))

  proj_rep = test_model$report()
  ## check all the projection devs are being populated
  expect_true(all(proj_rep$recruitment_devs[1,(length(data$years) + 1):(length(data$years) + data$n_projections_years)] != 0.0))

  ## empircially resample recruits
  data$future_recruitment_type = 1
  parameters$trans_rec_dev = matrix(rnorm(length(data$years), 0, 1.2), nrow = 1)
  ## make sure it doesn't crash
  expect_no_condition(test_model <- TMB::MakeADFun(data = data,
                                                   parameters = parameters,
                                                   DLL = "SpatialSablefishAssessment_TMBExports", silent  = T))

  proj_rep = test_model$report()
  proj_rep$recruitment_devs
  expect_true(all(proj_rep$recruitment_devs[1,(length(data$years) + 1):(length(data$years) + data$n_projections_years)] != 0.0))
  expect_true(all(proj_rep$recruitment_devs[1,(length(data$years) + 1):(length(data$years) + data$n_projections_years)] %in% proj_rep$recruitment_devs[1,1:length(data$years)]))

  ## resample from the first three years
  data$year_ndx_for_empirical_resampling = c(0,2)


  expect_no_condition(test_model <- TMB::MakeADFun(data = data,
                                                   parameters = parameters,
                                                   DLL = "SpatialSablefishAssessment_TMBExports", silent  = T))

  proj_rep = test_model$report()
  proj_rep$recruitment_devs
  expect_true(all(proj_rep$recruitment_devs[1,(length(data$years) + 1):(length(data$years) + data$n_projections_years)] != 0.0))
  expect_true(all(proj_rep$recruitment_devs[1,(length(data$years) + 1):(length(data$years) + data$n_projections_years)] %in% proj_rep$recruitment_devs[1,1:3]))


  ## resample from the last three years
  data$year_ndx_for_empirical_resampling = c(8,10)


  expect_no_condition(test_model <- TMB::MakeADFun(data = data,
                                                   parameters = parameters,
                                                   DLL = "SpatialSablefishAssessment_TMBExports", silent  = T))

  proj_rep = test_model$report()
  proj_rep$recruitment_devs
  expect_true(all(proj_rep$recruitment_devs[1,(length(data$years) + 1):(length(data$years) + data$n_projections_years)] != 0.0))
  expect_true(all(proj_rep$recruitment_devs[1,(length(data$years) + 1):(length(data$years) + data$n_projections_years)] %in% proj_rep$recruitment_devs[1,9:11]))


  ## test production model projection ability
  data$model = "TagIntegrated"
  ## make sure it doesn't crash
  expect_no_condition(test_model <- TMB::MakeADFun(data = data,
                                                   parameters = parameters,
                                                   DLL = "SpatialSablefishAssessment_TMBExports", silent  = T))

})

