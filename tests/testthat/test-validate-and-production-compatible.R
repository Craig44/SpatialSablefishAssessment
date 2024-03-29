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

  ## turn on observations to compare likelihood evaluations
  ## assumes multinomial
  data$trwl_catchatlgth_indicator[data$trwl_catchatlgth_indicator == 0] = 1
  data$obs_trwl_catchatlgth = array(5, dim = dim(data$obs_trwl_catchatlgth))
  data$fixed_catchatlgth_indicator[data$fixed_catchatlgth_indicator == 0] = 1
  data$obs_fixed_catchatlgth = array(5, dim = dim(data$obs_fixed_catchatlgth))
  data$fixed_catchatage_indicator[data$fixed_catchatage_indicator == 0] = 1
  data$obs_fixed_catchatage = array(5, dim = dim(data$obs_fixed_catchatage))
  data$srv_catchatage_indicator[data$srv_catchatage_indicator == 0] = 1
  data$obs_srv_catchatage = array(5, dim = dim(data$obs_srv_catchatage))


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

  # check Bzero and Binit calcualtions are consistent when F_init = 0
  expect_equal(production_report$Bzero, production_report$Binit, tolerance = 0.0001)
  ## get tag recoveries
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
  ## partition - male
  expect_true(all(validate_report$natage_m == production_report$natage_m))
  ## partition - female
  expect_true(all(validate_report$natage_f == production_report$natage_f))

})

#' compatibility-single-release-with-movement-and-Z-dirichletmultinomial
#' @description this will run a fairly complex tag module and check nll and predicted values
#' are the same between the validation TMB model and production TMB model
#'
test_that("compatibility-single-release-with-movement-and-Z-dirichletmultinomial", {
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

  ## turn on observations to compare likelihood evaluations
  ## assumes multinomial
  data$trwl_catchatlgth_indicator[data$trwl_catchatlgth_indicator == 0] = 1
  data$obs_trwl_catchatlgth = array(5, dim = dim(data$obs_trwl_catchatlgth))
  data$fixed_catchatlgth_indicator[data$fixed_catchatlgth_indicator == 0] = 1
  data$obs_fixed_catchatlgth = array(6, dim = dim(data$obs_fixed_catchatlgth))
  data$fixed_catchatage_indicator[data$fixed_catchatage_indicator == 0] = 1
  data$obs_fixed_catchatage = array(7, dim = dim(data$obs_fixed_catchatage))
  data$srv_catchatage_indicator[data$srv_catchatage_indicator == 0] = 1
  data$obs_srv_catchatage = array(8, dim = dim(data$obs_srv_catchatage))

  data$trwl_catchatlgth_indicator[data$trwl_catchatlgth_indicator == 0] = 1
  data$fixed_catchatlgth_indicator[data$fixed_catchatlgth_indicator == 0] = 1
  data$fixed_catchatage_comp_likelihood = 1
  data$srv_catchatage_comp_likelihood = 1

  parameters$logistic_tag_reporting_rate = matrix(logit(0.2), nrow = data$n_regions, ncol = sum(data$tag_recovery_indicator))
  parameters$trans_fixed_catchatlgth_error = log(0.1)
  parameters$trans_theta_trwl_catchatlgth_error = log(0.2)
  parameters$trans_fixed_catchatage_error = log(10)
  parameters$trans_srv_catchatage = log(7.4)

  validate_model <- TMB::MakeADFun(data = data,
                                   parameters = parameters,
                                   DLL = "SpatialSablefishAssessment_TMBExports", silent  = T)

  data$model = "TagIntegrated"
  production_model <- TMB::MakeADFun(data = data,
                                     parameters = parameters,
                                     DLL = "SpatialSablefishAssessment_TMBExports", silent  = T)

  validate_report = validate_model$report()
  production_report = production_model$report()

  ## check likelihood contribution is the same
  expect_true(sum(validate_report$nll) == sum(production_report$nll))
})


#' compatibility-single-release-with-movement-and-Z-tag-likelihood-multinomial-release
#' @description this will run a fairly complex tag module and check nll and predicted values
#' are the same between the validation TMB model and production TMB model
#'
test_that("compatibility-single-release-with-movement-and-Z-tag-likelihood-multinomial-release", {
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

  data$tag_recovery_indicator_by_year = rep(0, length(data$years)) ## no tag releases
  data$tag_recovery_indicator = array(0, dim = c(length(data$years), data$n_regions))
  data$tag_recovery_indicator[1,] = 1
  data$obs_tag_recovery = array(10, dim = c(data$n_regions * data$n_years_to_retain_tagged_cohorts_for + 1, data$n_regions, length(data$years)))
  for(r in 1:data$n_regions)
    data$obs_tag_recovery[,r,1] = (data$obs_tag_recovery[,r,1] / sum(data$obs_tag_recovery[,r,1])) * (sum(data$male_tagged_cohorts_by_age[,r,1]) + sum(data$female_tagged_cohorts_by_age[,r,1]))

  data$tag_likelihood = 2

  ## turn on observations to compare likelihood evaluations
  ## assumes multinomial
  data$trwl_catchatlgth_indicator[data$trwl_catchatlgth_indicator == 0] = 1
  data$obs_trwl_catchatlgth = array(5, dim = dim(data$obs_trwl_catchatlgth))
  data$fixed_catchatlgth_indicator[data$fixed_catchatlgth_indicator == 0] = 1
  data$obs_fixed_catchatlgth = array(6, dim = dim(data$obs_fixed_catchatlgth))
  data$fixed_catchatage_indicator[data$fixed_catchatage_indicator == 0] = 1
  data$obs_fixed_catchatage = array(7, dim = dim(data$obs_fixed_catchatage))
  data$srv_catchatage_indicator[data$srv_catchatage_indicator == 0] = 1
  data$obs_srv_catchatage = array(8, dim = dim(data$obs_srv_catchatage))

  data$trwl_catchatlgth_indicator[data$trwl_catchatlgth_indicator == 0] = 1
  data$fixed_catchatlgth_indicator[data$fixed_catchatlgth_indicator == 0] = 1
  data$fixed_catchatage_comp_likelihood = 1
  data$srv_catchatage_comp_likelihood = 1

  parameters$logistic_tag_reporting_rate = matrix(logit(0.2), nrow = data$n_regions, ncol = sum(data$tag_recovery_indicator))
  parameters$trans_fixed_catchatlgth_error = log(0.1)
  parameters$trans_theta_trwl_catchatlgth_error = log(0.2)
  parameters$trans_fixed_catchatage_error = log(10)
  parameters$trans_srv_catchatage = log(7.4)

  validate_model <- TMB::MakeADFun(data = data,
                                   parameters = parameters,
                                   DLL = "SpatialSablefishAssessment_TMBExports", silent  = T)

  data$model = "TagIntegrated"
  production_model <- TMB::MakeADFun(data = data,
                                     parameters = parameters,
                                     DLL = "SpatialSablefishAssessment_TMBExports", silent  = T)

  validate_report = validate_model$report()
  production_report = production_model$report()

  ## check likelihood contribution is the same
  expect_true(sum(validate_report$nll) == sum(production_report$nll))
})
