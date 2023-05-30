#'
#' These unit-tests will run some basic unit-tests on the catch at age and associated likelihoods
#'

#'
#' AgeBasedMovement
#'
test_that("AgeBasedMovementTagIntegratedModel", {
  ## Read in mock data
  load(system.file("testdata", "MockAgeBasedMovementModel.RData",package="SpatialSablefishAssessment"))
  data$model = "TagIntegratedAgeBasedMovement"
  ## no Z or movement
  data$apply_fixed_movement = 1
  data$age_based_movement = 1
  data$fixed_movement_matrix_young = data$fixed_movement_matrix
  data$fixed_movement_matrix_old = data$fixed_movement_matrix
  data$tag_likelihood = 2
  data$obs_tag_recovery = array(10, dim = c(length(data$ages), dim(data$obs_tag_recovery)))
  parameters$transformed_movement_pars_young = parameters$transformed_movement_pars
  parameters$transformed_movement_pars_old = parameters$transformed_movement_pars
  parameters$ln_a50_movement = log(6)
  parameters$ln_ato95_movement = log(1)

  test_model <- TMB::MakeADFun(data = data,
                               parameters = parameters,
                               DLL = "SpatialSablefishAssessment_TMBExports", silent  = T)

  test_report = test_model$report()

  ## check the age based ogive works as expected
  expect_equal(test_report$old_age_based_movement_ogive, logis(data$ages, exp(parameters$ln_a50_movement), exp(parameters$ln_ato95_movement)), tolerance = 0.001)

  ## check  both young and old matrices are being built as expected
  expect_equal(test_report$movement_matrix_young, test_report$movement_matrix_old)

})

#'
#' AgeBasedMovement-noagebasedMovement
#'
test_that("AgeBasedMovementTagIntegratedModel-noagebasedMovement-no-move", {
  ## Read in mock data
  load(system.file("testdata", "MockSablefishModel.RData",package="SpatialSablefishAssessment"))

  data$model = "TagIntegrated"
  data$apply_fixed_movement = 1

  data$tag_likelihood = 0 ## Poisson

  test_model_int <- TMB::MakeADFun(data = data,
                               parameters = parameters,
                               DLL = "SpatialSablefishAssessment_TMBExports", silent  = T)

  test_report_int = test_model_int$report()

  load(system.file("testdata", "MockAgeBasedMovementModel.RData",package="SpatialSablefishAssessment"))
  data$model = "TagIntegratedAgeBasedMovement"
  ## no Z or movement
  data$age_based_movement = 0
  data$fixed_movement_matrix_young = data$fixed_movement_matrix
  data$fixed_movement_matrix_old = data$fixed_movement_matrix
  data$obs_tag_recovery = array(10/length(data$ages), dim = c(length(data$ages), dim(data$obs_tag_recovery)))
  parameters$transformed_movement_pars_young = parameters$transformed_movement_pars
  parameters$transformed_movement_pars_old = parameters$transformed_movement_pars
  parameters$ln_a50_movement = log(6)
  parameters$ln_ato95_movement = log(1)

  test_model <- TMB::MakeADFun(data = data,
                               parameters = parameters,
                               DLL = "SpatialSablefishAssessment_TMBExports", silent  = T)

  test_report = test_model$report()

  #round(cbind(test_report$nll, test_report_int$nll), 2)

  expect_equal(test_report$nll[8], test_report_int$nll[8], tolerance = 0.0001)

  temp = test_report$pred_aggregated_tag_recovery[1,,,]
  expect_equal(sum(temp[data$tag_recovery_indicator == 1]), sum(test_report_int$pred_tag_recovery[data$tag_recovery_indicator == 1]), tolerance = 0.0001)


})

#'
#' AgeBasedMovement-noagebasedMovement
#'
test_that("AgeBasedMovementTagIntegratedModel-noagebasedMovement-move", {
  ## Read in mock data
  load(system.file("testdata", "MockSablefishModel.RData",package="SpatialSablefishAssessment"))

  data$model = "TagIntegrated"
  data$apply_fixed_movement = 0

  data$tag_likelihood = 0 ## Poisson

  test_model_int <- TMB::MakeADFun(data = data,
                                   parameters = parameters,
                                   DLL = "SpatialSablefishAssessment_TMBExports", silent  = T)

  test_report_int = test_model_int$report()

  load(system.file("testdata", "MockAgeBasedMovementModel.RData",package="SpatialSablefishAssessment"))
  data$model = "TagIntegratedAgeBasedMovement"
  ## no Z or movement
  data$age_based_movement = 0
  data$apply_fixed_movement = 0

  data$tag_likelihood = 0 ## Poisson
  data$fixed_movement_matrix_young = data$fixed_movement_matrix
  data$fixed_movement_matrix_old = data$fixed_movement_matrix
  data$obs_tag_recovery = array(10/length(data$ages), dim = c(length(data$ages), dim(data$obs_tag_recovery)))
  parameters$transformed_movement_pars_young = parameters$transformed_movement_pars
  parameters$transformed_movement_pars_old = parameters$transformed_movement_pars
  parameters$ln_a50_movement = log(6)
  parameters$ln_ato95_movement = log(1)

  test_model <- TMB::MakeADFun(data = data,
                               parameters = parameters,
                               DLL = "SpatialSablefishAssessment_TMBExports", silent  = T)

  test_report = test_model$report()

  #round(cbind(test_report$nll, test_report_int$nll), 2)

  expect_equal(test_report$nll[8], test_report_int$nll[8], tolerance = 0.0001)
  expect_equal(sum(test_report$nll), sum(test_report_int$nll), tolerance = 0.0001)
  temp = test_report$pred_aggregated_tag_recovery[1,,,]

  expect_equal(sum(temp[data$tag_recovery_indicator == 1]), sum(test_report_int$pred_tag_recovery[data$tag_recovery_indicator == 1]), tolerance = 0.0001)

})
