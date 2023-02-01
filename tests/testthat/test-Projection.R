

#' test-Projection-Recruitment
#' @description this is to test projected recruitment works as expected
#'
test_that("test-Projection-Recruitment", {
  ## Read in mock data
  load(system.file("testdata", "MockProjectionData.RData", package="SpatialSablefishAssessment"))
  ## no Z or movement
  data$apply_fixed_movement = 1 ## no movement initial age-structure should be only a function of ageing and M
  data$model = "TagIntegrated" # so we can run the validate function
  expect_true(validate_input_data_and_parameters(data, parameters))
  data$model = "TagIntegratedValidate"
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
  data$model = "TagIntegrated" # so we can run the validate function
  expect_true(validate_input_data_and_parameters(data, parameters))
  data$model = "TagIntegratedValidate"

  expect_no_condition(test_model <- TMB::MakeADFun(data = data,
                                                   parameters = parameters,
                                                   DLL = "SpatialSablefishAssessment_TMBExports", silent  = T))

  proj_rep = test_model$report()
  proj_rep$recruitment_devs
  expect_true(all(proj_rep$recruitment_devs[1,(length(data$years) + 1):(length(data$years) + data$n_projections_years)] != 0.0))
  expect_true(all(proj_rep$recruitment_devs[1,(length(data$years) + 1):(length(data$years) + data$n_projections_years)] %in% proj_rep$recruitment_devs[1,1:length(data$years)]))

  ## resample from the first three years
  data$year_ndx_for_empirical_resampling = c(0,2)

  data$model = "TagIntegrated" # so we can run the validate function
  expect_true(validate_input_data_and_parameters(data, parameters))
  data$model = "TagIntegratedValidate"

  expect_no_condition(test_model <- TMB::MakeADFun(data = data,
                                                   parameters = parameters,
                                                   DLL = "SpatialSablefishAssessment_TMBExports", silent  = T))

  proj_rep = test_model$report()
  proj_rep$recruitment_devs
  expect_true(all(proj_rep$recruitment_devs[1,(length(data$years) + 1):(length(data$years) + data$n_projections_years)] != 0.0))
  expect_true(all(proj_rep$recruitment_devs[1,(length(data$years) + 1):(length(data$years) + data$n_projections_years)] %in% proj_rep$recruitment_devs[1,1:3]))


  ## resample from the last three years
  data$year_ndx_for_empirical_resampling = c(8,10)

  data$model = "TagIntegrated" # so we can run the validate function
  expect_true(validate_input_data_and_parameters(data, parameters))
  data$model = "TagIntegratedValidate"

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



  ####### test constant recruitment
  data$future_recruitment_type = 2
  expect_no_condition(test_model <- TMB::MakeADFun(data = data,
                                                   parameters = parameters,
                                                   DLL = "SpatialSablefishAssessment_TMBExports", silent  = T))
  proj_rep = test_model$report()
  ## check multiplers are all = 1.0
  expect_equal(proj_rep$recruitment_multipliers[1,(length(data$years) + 1):(length(data$years) + data$n_projections_years)], rep(1, data$n_projections_years), tolerance = 0.0001)

})



#' test-Projection-Mortality
#' @description this is to test projected F methods work as expected
#'
test_that("test-Projection-Mortality", {
  ## Read in mock data
  load(system.file("testdata", "MockProjectionData.RData",package="SpatialSablefishAssessment"))
  data$model = "TagIntegratedValidate"
  ## no Z or movement
  data$apply_fixed_movement = 1 ## no movement initial age-structure should be only a function of ageing and M
  data$future_fishing_type = 0
  data$future_recruitment_type = 2
  ## make sure it doesn't crash
  expect_no_condition(test_model <- TMB::MakeADFun(data = data,
                                                   parameters = parameters,
                                                   DLL = "SpatialSablefishAssessment_TMBExports", silent  = T))

  proj_rep = test_model$report()
  ## check Future Fs are the same as the input
  expect_equal(data$future_fishing_inputs_fixed, proj_rep$annual_F_fixed[, (length(data$years) + 1):(length(data$years) + data$n_projections_years)], tolerance = 0.00001)
  expect_equal(data$future_fishing_inputs_trwl, proj_rep$annual_F_trwl[, (length(data$years) + 1):(length(data$years) + data$n_projections_years)], tolerance = 0.00001)


  data$future_fishing_type = 1
  data$future_fishing_inputs_trwl = array(mean(proj_rep$annual_trwl_catch_pred), dim = c(data$n_regions, data$n_projections_years))
  data$future_fishing_inputs_fixed = array(mean(proj_rep$annual_fixed_catch_pred), dim = c(data$n_regions, data$n_projections_years))

  ## make sure it doesn't crash
  expect_no_condition(test_model <- TMB::MakeADFun(data = data,
                                                   parameters = parameters,
                                                   DLL = "SpatialSablefishAssessment_TMBExports", silent  = T))

  proj_rep = test_model$report()
  ## check predicted catch is close to input catch
  expect_equal(data$future_fishing_inputs_fixed, proj_rep$annual_fixed_catch_pred[, (length(data$years) + 1):(length(data$years) + data$n_projections_years)], tolerance = 0.01)
  expect_equal(data$future_fishing_inputs_trwl, proj_rep$annual_trwl_catch_pred[, (length(data$years) + 1):(length(data$years) + data$n_projections_years)], tolerance = 0.01)

  ## check  Fs are greater than zero
  expect_true(all(proj_rep$annual_F_fixed[, (length(data$years) + 1):(length(data$years) + data$n_projections_years)] > 0))
  expect_true(all(proj_rep$annual_F_trwl[, (length(data$years) + 1):(length(data$years) + data$n_projections_years)] > 0))

})
