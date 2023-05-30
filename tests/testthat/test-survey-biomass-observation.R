#'
#' These unit-tests will run some basic unit-tests survey biomass observation and associated likelihoods
#'

#'
#' abundance-test
#'
test_that("survey-abundance-test", {
  ## Read in mock data
  load(system.file("testdata", "MockSablefishModel.RData",package="SpatialSablefishAssessment"))
  data$model = "TagIntegratedValidate"
  ## no Z or movement
  data$apply_Z_on_tagged_fish = 0
  data$apply_fixed_movement = 1
  data$srv_bio_indicator[data$srv_bio_indicator == 0] = 1
  data$srv_obs_is_abundance = rep(1, data$n_surveys)
  mean_obs_val = 12
  data$obs_srv_bio = array(mean_obs_val, dim = c(data$n_regions, ncol(data$srv_bio_indicator), data$n_surveys))
  data$obs_srv_se = array(0.3 * mean_obs_val, dim = c(data$n_regions, ncol(data$srv_bio_indicator), data$n_surveys))
  data$srv_bio_likelihood = rep(0, data$n_surveys)
  data$q_is_nuisance = rep(0, data$n_surveys)

  ## make sure it doesn't crash
  expect_no_condition(test_model <- TMB::MakeADFun(data = data,
                                                   parameters = parameters,
                                                   DLL = "SpatialSablefishAssessment_TMBExports", silent  = T))
  test_report = test_model$report()

  expect_true(all(test_report$pred_srv_bio > 0))

  ## change likelihood and thus SE intepretation
  data$obs_srv_se = array(0.3, dim = c(data$n_regions, ncol(data$srv_bio_indicator), data$n_surveys))
  data$srv_bio_likelihood = rep(1, data$n_surveys)
  expect_no_condition(test_model <- TMB::MakeADFun(data = data,
                                                   parameters = parameters,
                                                   DLL = "SpatialSablefishAssessment_TMBExports", silent  = T))

  test_report = test_model$report()
  ## test against R's dlnorm function they should be the same
  expect_equal(test_report$nll[5], -1.0 * sum(dlnorm(test_report$obs_srv_bio, log(test_report$pred_srv_bio) - 0.5*test_report$obs_srv_se^2, test_report$obs_srv_se, log = T)), tolerance = 0.001)

  ## nuisance q calculation
  data$q_is_nuisance = rep(1, data$n_surveys)
  ## check it doesn't crash
  expect_no_condition(test_model <- TMB::MakeADFun(data = data,
                                                   parameters = parameters,
                                                   DLL = "SpatialSablefishAssessment_TMBExports", silent  = T))
  test_report = test_model$report()


  ## repeat but for weight instead of abundance
  data$srv_obs_is_abundance = rep(0, data$n_surveys)
  data$srv_bio_likelihood = rep(0, data$n_surveys)
  data$q_is_nuisance = rep(1, data$n_surveys)
  mean_obs_val = 45
  data$obs_srv_bio = array(mean_obs_val, dim = c(data$n_regions, ncol(data$srv_bio_indicator), data$n_surveys))
  data$obs_srv_se = array(0.3 * mean_obs_val, dim = c(data$n_regions, ncol(data$srv_bio_indicator), data$n_surveys))
  ## check it doesn't crash
  expect_no_condition(test_model <- TMB::MakeADFun(data = data,
                                                   parameters = parameters,
                                                   DLL = "SpatialSablefishAssessment_TMBExports", silent  = T))
  test_report = test_model$report()
  test_report = test_model$report()

  expect_true(all(test_report$pred_srv_dom_ll_bio > 0))

  ## change likelihood and thus SE intepretation
  data$obs_srv_se = array(0.3, dim = c(data$n_regions, ncol(data$srv_bio_indicator),data$n_surveys))
  data$srv_bio_likelihood = rep(1, data$n_surveys)
  expect_no_condition(test_model <- TMB::MakeADFun(data = data,
                                                   parameters = parameters,
                                                   DLL = "SpatialSablefishAssessment_TMBExports", silent  = T))

  test_report = test_model$report()
  ## test against R's dlnorm function they should be the same
  expect_equal(test_report$nll[5], -1.0 * sum(dlnorm(test_report$obs_srv_bio, log(test_report$pred_srv_bio) - 0.5*test_report$obs_srv_se^2, test_report$obs_srv_se, log = T)), tolerance = 0.001)

  ## nuisance q calculation
  data$q_is_nuisance = rep(1, data$n_surveys)
  ## check it doesn't crash
  expect_no_condition(test_model <- TMB::MakeADFun(data = data,
                                                   parameters = parameters,
                                                   DLL = "SpatialSablefishAssessment_TMBExports", silent  = T))
  test_report = test_model$report()

})

#'
#' abundance-test-production
#'
test_that("survey-abundance-production-test", {
  ## Read in mock data
  load(system.file("testdata", "MockSablefishModel.RData",package="SpatialSablefishAssessment"))
  data$model = "TagIntegrated"
  ## no Z or movement
  data$apply_Z_on_tagged_fish = 0
  data$apply_fixed_movement = 1
  data$srv_bio_indicator[data$srv_bio_indicator == 0] = 1
  data$srv_obs_is_abundance = rep(0, data$n_surveys)
  mean_obs_val = 12
  data$obs_srv_bio = array(mean_obs_val, dim = c(data$n_regions, ncol(data$srv_bio_indicator), data$n_surveys))
  data$obs_srv_se = array(0.3 * mean_obs_val, dim = c(data$n_regions, ncol(data$srv_bio_indicator), data$n_surveys))
  data$srv_dom_ll_bio_comp_likelihood = rep(0, data$n_surveys)
  data$q_is_nuisance = rep(0, data$n_surveys)

  ## make sure it doesn't crash
  expect_no_condition(test_model <- TMB::MakeADFun(data = data,
                                                   parameters = parameters,
                                                   DLL = "SpatialSablefishAssessment_TMBExports", silent  = T))
  test_report = test_model$report()

  expect_true(all(test_report$pred_srv_bio > 0))

  ## change likelihood and thus SE intepretation
  data$obs_srv_se = array(0.3, dim = c(data$n_regions, ncol(data$srv_bio_indicator), data$n_surveys))
  data$srv_bio_comp_likelihood = rep(1, data$n_surveys)
  expect_no_condition(test_model <- TMB::MakeADFun(data = data,
                                                   parameters = parameters,
                                                   DLL = "SpatialSablefishAssessment_TMBExports", silent  = T))

  test_report = test_model$report()
  ## test against R's dlnorm function they should be the same
  expect_equal(test_report$nll[5], -1.0 * sum(dlnorm(test_report$obs_srv_bio, log(test_report$pred_srv_bio) - 0.5*test_report$obs_srv_se^2, test_report$obs_srv_se, log = T)), tolerance = 0.001)

  ## nuisance q calculation
  data$q_is_nuisance = rep(1, data$n_surveys)
  ## check it doesn't crash
  expect_no_condition(test_model <- TMB::MakeADFun(data = data,
                                                   parameters = parameters,
                                                   DLL = "SpatialSablefishAssessment_TMBExports", silent  = T))
  test_report = test_model$report()


  ## repeat but for weight instead of abundance
  data$srv_obs_is_abundance = rep(0, data$n_surveys)
  data$srv_bio_comp_likelihood = rep(0, data$n_surveys)
  data$q_is_nuisance = rep(0, data$n_surveys)
  mean_obs_val = 45
  data$obs_srv_bio = array(mean_obs_val, dim = c(data$n_regions, ncol(data$srv_bio_indicator), data$n_surveys))
  data$obs_srv_se = array(0.3 * mean_obs_val, dim = c(data$n_regions, ncol(data$srv_bio_indicator), data$n_surveys))
  ## check it doesn't crash
  expect_no_condition(test_model <- TMB::MakeADFun(data = data,
                                                   parameters = parameters,
                                                   DLL = "SpatialSablefishAssessment_TMBExports", silent  = T))
  test_report = test_model$report()

  expect_true(all(test_report$pred_srv_bio > 0))

  ## change likelihood and thus SE intepretation
  data$obs_srv_se = array(0.3, dim = c(data$n_regions, ncol(data$srv_bio_indicator), data$n_surveys))
  data$srv_bio_comp_likelihood = rep(1, data$n_surveys)
  expect_no_condition(test_model <- TMB::MakeADFun(data = data,
                                                   parameters = parameters,
                                                   DLL = "SpatialSablefishAssessment_TMBExports", silent  = T))

  test_report = test_model$report()
  ## test against R's dlnorm function they should be the same
  expect_equal(test_report$nll[5], -1.0 * sum(dlnorm(test_report$obs_srv_bio, log(test_report$pred_srv_bio) - 0.5*test_report$obs_srv_se^2, test_report$obs_srv_se, log = T)), tolerance = 0.001)

  ## nuisance q calculation
  data$q_is_nuisance = rep(1, data$n_surveys)
  ## check it doesn't crash
  expect_no_condition(test_model <- TMB::MakeADFun(data = data,
                                                   parameters = parameters,
                                                   DLL = "SpatialSablefishAssessment_TMBExports", silent  = T))
  test_report = test_model$report()

})
