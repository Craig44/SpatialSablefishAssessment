#'
#' These unit-tests will run some basic unit-tests on the catch at age and associated likelihoods
#'

#'
#' AgeBasedMovement
#'
test_that("AgeBasedMovement", {
  ## Read in mock data
  load(system.file("testdata", "MockSablefishModel.RData",package="SpatialSablefishAssessment"))
  data$model = "TagIntegratedAgeBasedMovement"
  ## no Z or movement
  data$apply_fixed_movement = 1
  data$age_based_movement = 1
  data$fixed_movement_matrix_young = data$fixed_movement_matrix
  data$fixed_movement_matrix_old = data$fixed_movement_matrix


  parameters$transformed_movement_pars_young = parameters$transformed_movement_pars
  parameters$transformed_movement_pars_old = parameters$transformed_movement_pars
  parameters$ln_a50_movement = log(6)
  parameters$ln_ato95_movement = log(1)

  test_model <- TMB::MakeADFun(data = data,
                               parameters = parameters,
                               DLL = "SpatialSablefishAssessment_TMBExports", silent  = T)

  test_report = test_model$report()

  ## check the age based ogive works as expected
  expect_equal(test_report$old_age_based_movement_ogive, logis(data$ages, exp(parameters$ln_a50_movement), exp(parameters$ln_ato95_movement))/ max(logis(data$ages, exp(parameters$ln_a50_movement), exp(parameters$ln_ato95_movement))), tolerance = 0.001)

  ## check  both young and old matrices are being built as expected
  expect_equal(test_report$movement_matrix_young, test_report$movement_matrix_old)


})

