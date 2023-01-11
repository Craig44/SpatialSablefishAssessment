#'
#' These unit-tests will run some basic unit-tests on the catch at length using
#' simulated data
#'


#' single-release-no-movement-and-Z
#' @description this will track taged cohorts in a model with no movement or mortality
#' to check that tagged fish move as we expect. Because they have a slightly abstract
#' partition design this is a useful first validation
#'
test_that("catch-at-length-trwl", {
  ## Read in mock data
  load(system.file("testdata", "MockSablefishModel.RData",package="SpatialSablefishAssessment"))
  data$model = "TagIntegratedValidate"
  ## no Z or movement
  data$apply_Z_on_tagged_fish = 0
  data$apply_fixed_movement = 0
  data$trwl_catchatlgth_indicator[data$trwl_catchatlgth_indicator == 0] = 1
  data$obs_trwl_catchatlgth = array(5, dim = dim(data$obs_trwl_catchatlgth))
  data$fixed_catchatlgth_indicator[data$fixed_catchatlgth_indicator == 0] = 1
  data$obs_fixed_catchatlgth = array(5, dim = dim(data$obs_fixed_catchatlgth))
  ## multinomial
  test_model <- TMB::MakeADFun(data = data,
                               parameters = parameters,
                               DLL = "SpatialSablefishAssessment_TMBExports", silent  = T)

  test_report = test_model$report()

  ## Trawl
  F_pred_age = test_report$natage_f[,1,1] * test_report$F_trwl_f[,1,1] / test_report$Z_f[,1,1] * (1 - test_report$S_f[,1,1])
  M_pred_age = test_report$natage_m[,1,1] * test_report$F_trwl_m[,1,1] / test_report$Z_m[,1,1] * (1 - test_report$S_m[,1,1])

  ## check catch at age
  expect_equal(M_pred_age, test_report$catchatage_trwl_m[,1,1], tolerance = 0.001)
  expect_equal(F_pred_age, test_report$catchatage_trwl_f[,1,1], tolerance = 0.001)

  ## catch at length
  m_pred_at_length = M_pred_age %*% data$male_age_length_transition[,,1]
  f_pred_at_length = F_pred_age %*% data$female_age_length_transition[,,1]

  ##
  pred_full = c(m_pred_at_length, f_pred_at_length)
  ## test
  expect_equal(pred_full / sum(pred_full), test_report$pred_trwl_catchatlgth[,1,1], tolerance = 0.001)

  ## likelihood
  expected_nll = 0
  for(y in 1:(length(data$years))) {
    for(r in 1:data$n_regions) {
      expected_nll = expected_nll - sum(dmultinom_upd(x = test_report$obs_trwl_catchatlgth[,r,y], prob = test_report$pred_trwl_catchatlgth[,r,y], log = T))
    }
  }
  expect_equal(test_report$nll[2], expected_nll, tolerance = 0.1) ## with one log likelihood value

  ## Fixed gear

  F_pred_age = test_report$natage_f[,1,1] * test_report$F_fixed_f[,1,1] / test_report$Z_f[,1,1] * (1 - test_report$S_f[,1,1])
  M_pred_age = test_report$natage_m[,1,1] * test_report$F_fixed_m[,1,1] / test_report$Z_m[,1,1] * (1 - test_report$S_m[,1,1])

  ## check catch at age
  expect_equal(M_pred_age, test_report$catchatage_fixed_m[,1,1], tolerance = 0.001)
  expect_equal(F_pred_age, test_report$catchatage_fixed_f[,1,1], tolerance = 0.001)

  ## catch at length
  m_pred_at_length = M_pred_age %*% data$male_age_length_transition[,,1]
  f_pred_at_length = F_pred_age %*% data$female_age_length_transition[,,1]

  ##
  pred_full = c(m_pred_at_length, f_pred_at_length)
  ## test
  expect_equal(pred_full / sum(pred_full), test_report$pred_fixed_catchatlgth[,1,1], tolerance = 0.001)

  ## likelihood
  expected_nll = 0
  for(y in 1:(length(data$years))) {
    for(r in 1:data$n_regions) {
      expected_nll = expected_nll - sum(dmultinom_upd(x = test_report$obs_fixed_catchatlgth[,r,y], prob = test_report$pred_fixed_catchatlgth[,r,y], log = T))
    }
  }
  expect_equal(test_report$nll[3], expected_nll, tolerance = 0.1) ## with one log likelihood value


  ## Dirichlet-Multinomial
  data$fixed_catchatlgth_comp_likelihood = 1
  data$trwl_catchatlgth_comp_likelihood = 1
  test_model <- TMB::MakeADFun(data = data,
                               parameters = parameters,
                               DLL = "SpatialSablefishAssessment_TMBExports", silent  = T)

  test_report = test_model$report()
  ## Fixed gear
  expected_nll = 0
  for(y in 1:(length(data$years))) {
    for(r in 1:data$n_regions) {
      expected_nll = expected_nll - sum(ddirichmult(obs = test_report$obs_fixed_catchatlgth[,r,y] / sum(test_report$obs_fixed_catchatlgth[,r,y]), est = test_report$pred_fixed_catchatlgth[,r,y], beta = sum(test_report$obs_fixed_catchatlgth[,r,y]) * test_report$theta_fixed_catchatlgth, n = sum(test_report$obs_fixed_catchatlgth[,r,y]), log_it = T))
    }
  }
  expect_equal(test_report$nll[3], expected_nll, tolerance = 0.1) ## with one log likelihood value

  ## trawl gear
  expected_nll = 0
  for(y in 1:(length(data$years))) {
    for(r in 1:data$n_regions) {
      expected_nll = expected_nll - sum(ddirichmult(obs = test_report$obs_trwl_catchatlgth[,r,y] / sum(test_report$obs_trwl_catchatlgth[,r,y]), est = test_report$pred_trwl_catchatlgth[,r,y], beta = sum(test_report$obs_trwl_catchatlgth[,r,y]) * test_report$theta_trwl_catchatlgth, n = sum(test_report$obs_trwl_catchatlgth[,r,y]), log_it = T))
    }
  }
  expect_equal(test_report$nll[2], expected_nll, tolerance = 0.1) ## with one log likelihood value


  ## smaller theta
  parameters$trans_fixed_catchatlgth_error = log(0.1)
  parameters$trans_theta_trwl_catchatlgth_error = log(0.2)
  test_model <- TMB::MakeADFun(data = data,
                               parameters = parameters,
                               DLL = "SpatialSablefishAssessment_TMBExports", silent  = T)

  test_report = test_model$report()
  ## Fixed gear
  expected_nll = 0
  for(y in 1:(length(data$years))) {
    for(r in 1:data$n_regions) {
      expected_nll = expected_nll - sum(ddirichmult(obs = test_report$obs_fixed_catchatlgth[,r,y] / sum(test_report$obs_fixed_catchatlgth[,r,y]), est = test_report$pred_fixed_catchatlgth[,r,y], beta = sum(test_report$obs_fixed_catchatlgth[,r,y]) * test_report$theta_fixed_catchatlgth, n = sum(test_report$obs_fixed_catchatlgth[,r,y]), log_it = T))
    }
  }
  expect_equal(test_report$nll[3], expected_nll, tolerance = 0.1) ## with one log likelihood value

  ## trawl gear
  expected_nll = 0
  for(y in 1:(length(data$years))) {
    for(r in 1:data$n_regions) {
      expected_nll = expected_nll - sum(ddirichmult(obs = test_report$obs_trwl_catchatlgth[,r,y] / sum(test_report$obs_trwl_catchatlgth[,r,y]), est = test_report$pred_trwl_catchatlgth[,r,y], beta = sum(test_report$obs_trwl_catchatlgth[,r,y]) * test_report$theta_trwl_catchatlgth, n = sum(test_report$obs_trwl_catchatlgth[,r,y]), log_it = T))
    }
  }
  expect_equal(test_report$nll[2], expected_nll, tolerance = 0.1) ## with one log likelihood value

  ## Larger theta
  parameters$trans_fixed_catchatlgth_error = log(6.7)
  parameters$trans_theta_trwl_catchatlgth_error = log(8.3)
  test_model <- TMB::MakeADFun(data = data,
                               parameters = parameters,
                               DLL = "SpatialSablefishAssessment_TMBExports", silent  = T)

  test_report = test_model$report()
  ## Fixed gear
  expected_nll = 0
  for(y in 1:(length(data$years))) {
    for(r in 1:data$n_regions) {
      expected_nll = expected_nll - sum(ddirichmult(obs = test_report$obs_fixed_catchatlgth[,r,y] / sum(test_report$obs_fixed_catchatlgth[,r,y]), est = test_report$pred_fixed_catchatlgth[,r,y], beta = sum(test_report$obs_fixed_catchatlgth[,r,y]) * test_report$theta_fixed_catchatlgth, n = sum(test_report$obs_fixed_catchatlgth[,r,y]), log_it = T))
    }
  }
  expect_equal(test_report$nll[3], expected_nll, tolerance = 0.1) ## with one log likelihood value

  ## trawl gear
  expected_nll = 0
  for(y in 1:(length(data$years))) {
    for(r in 1:data$n_regions) {
      expected_nll = expected_nll - sum(ddirichmult(obs = test_report$obs_trwl_catchatlgth[,r,y] / sum(test_report$obs_trwl_catchatlgth[,r,y]), est = test_report$pred_trwl_catchatlgth[,r,y], beta = sum(test_report$obs_trwl_catchatlgth[,r,y]) * test_report$theta_trwl_catchatlgth, n = sum(test_report$obs_trwl_catchatlgth[,r,y]), log_it = T))
    }
  }
  expect_equal(test_report$nll[2], expected_nll, tolerance = 0.1) ## with one log likelihood value
})



