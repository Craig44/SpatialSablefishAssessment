#'
#' Check that the Validate model and production model are compatible
#'


#' compatibility-single-release-with-movement-and-Z
#' @description this will run a fairly complex tag module and check nll and predicted values
#' are the same between the validation TMB model and production TMB model
#'
test_that("compatibility-single-release-with-movement-and-Z", {
  ## Read in mock data
  load(system.file("testdata", "MockSablefishModel.RData",package="SpatialSablefishAssessment"))

  data$model = "TagIntegratedValidate"
  ## no Z or movement
  data$apply_Z_on_tagged_fish = 1
  data$apply_fixed_movement = 1
  data$apply_tag_reporting_rate = 1
  ## this assumes no movement
  data$fixed_movement_matrix  = data$movement_matrix
  data$apply_fishery_tag_reporting = 1 ## all tagged fish will be recovered not a function of F
  ## turn off tag shedding and initial mortality
  data$initial_tag_induced_mortality = rep(0.1, sum(data$tag_release_event_this_year))
  data$annual_tag_shedding_rate = 0.05
  parameters$logistic_tag_reporting_rate = matrix(logit(0.2), nrow = data$n_regions, ncol = sum(data$tag_recovery_indicator))


  validate_model <- TMB::MakeADFun(data = data,
                         parameters = parameters,
                         DLL = "SpatialSablefishAssessment_TMBExports", silent  = T)

  data$model = "TagIntegrated"
  production_model <- TMB::MakeADFun(data = data,
                                   parameters = parameters,
                                   DLL = "SpatialSablefishAssessment_TMBExports", silent  = T)

  validate_report = validate_model$report()
  production_report = production_model$report()

  validate_predicted_recoveries = plot_tag_recovery_obs(validate_report, region_key = region_key, sex = "both", release_ndx_to_plot = 1:5)
  validate_predicted_recoveries = validate_predicted_recoveries$data
  production_predicted_recoveries = plot_tag_recovery_obs(production_report, region_key = region_key, sex = "both", release_ndx_to_plot = 1:5)
  production_predicted_recoveries = production_predicted_recoveries$data
  ## check they are reporting the same number of recoveries
  expect_true(nrow(validate_predicted_recoveries) == nrow(production_predicted_recoveries))
  ## check the predicted recoveries are the same
  expect_true(all(validate_predicted_recoveries$predicted == production_predicted_recoveries$predicted))

  ## check the tagged partitions at the end of the model run are the same.
  expect_true(all(validate_report$tagged_natage_f == production_report$tagged_natage_f))
  expect_true(all(validate_report$tagged_natage_m == production_report$tagged_natage_m))

  ## check likelihood contribution is the same
  expect_true(sum(validate_report$nll) == sum(production_report$nll))

  ## check SSBs
  expect_true(all(validate_report$SSB_yr == production_report$SSB_yr))
})

