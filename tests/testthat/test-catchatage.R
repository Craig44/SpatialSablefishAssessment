#'
#' These unit-tests will run some basic unit-tests on the catch at age and associated likelihoods
#'

#'
#' catch-at-age-fixed-multinomial
#'
test_that("catch-at-age-fixed-multinomial", {
  ## Read in mock data
  load(system.file("testdata", "MockSablefishModel.RData",package="SpatialSablefishAssessment"))
  data$model = "TagIntegratedValidate"
  ## no Z or movement
  data$apply_Z_on_tagged_fish = 0
  data$apply_fixed_movement = 1
  data$fixed_catchatage_indicator[data$fixed_catchatage_indicator == 0] = 1
  data$obs_fixed_catchatage = array(5, dim = dim(data$obs_fixed_catchatage))
  data$fixed_catchatage_comp_likelihood = 0
  data$srv_dom_ll_catchatage_indicator[data$srv_dom_ll_catchatage_indicator == 0] = 1
  data$obs_srv_dom_ll_catchatage = array(10, dim = dim(data$obs_srv_dom_ll_catchatage))
  data$srv_dom_ll_catchatage_comp_likelihood = 0
  test_model <- TMB::MakeADFun(data = data,
                               parameters = parameters,
                               DLL = "SpatialSablefishAssessment_TMBExports", silent  = T)

  test_report = test_model$report()

  F_pred_age = test_report$natage_f[,1,1] * test_report$F_fixed_f[,1,1] / test_report$Z_f[,1,1] * (1 - test_report$S_f[,1,1])
  M_pred_age = test_report$natage_m[,1,1] * test_report$F_fixed_m[,1,1] / test_report$Z_m[,1,1] * (1 - test_report$S_m[,1,1])

  ## check catch at age
  expect_equal(M_pred_age, test_report$catchatage_fixed_m[,1,1], tolerance = 0.001)
  expect_equal(F_pred_age, test_report$catchatage_fixed_f[,1,1], tolerance = 0.001)

  ## fixed gear likelihood
  expected_nll = 0
  for(y in 1:(length(data$years))) {
    for(r in 1:data$n_regions) {
      expected_nll = expected_nll - sum(dmultinom_upd(x = test_report$obs_fixed_catchatage[,r,y], prob = test_report$pred_fixed_catchatage[,r,y], log = T))
    }
  }
  expect_equal(test_report$nll[1], expected_nll, tolerance = 0.1) ## with one log likelihood value

  ## Survey likelihood
  expected_nll = 0
  for(y in 1:(length(data$years))) {
    for(r in 1:data$n_regions) {
      expected_nll = expected_nll - sum(dmultinom_upd(x = test_report$obs_srv_dom_ll_catchatage[,r,y], prob = test_report$pred_srv_dom_ll_catchatage[,r,y], log = T))
    }
  }
  expect_equal(test_report$nll[4], expected_nll, tolerance = 0.1) ## with one log likelihood value

})


#'
#' catch-at-age-fixed-dirichlet-multinomial
#'
test_that("catch-at-age-fixed-dirichlet-multinomial", {
  ## Read in mock data
  load(system.file("testdata", "MockSablefishModel.RData",package="SpatialSablefishAssessment"))
  data$model = "TagIntegratedValidate"
  ## no Z or movement
  data$apply_Z_on_tagged_fish = 0
  data$apply_fixed_movement = 1
  data$fixed_catchatage_indicator[data$fixed_catchatage_indicator == 0] = 1
  data$obs_fixed_catchatage = array(5, dim = dim(data$obs_fixed_catchatage))
  data$fixed_catchatage_comp_likelihood = 1
  data$srv_dom_ll_catchatage_indicator[data$srv_dom_ll_catchatage_indicator == 0] = 1
  data$obs_srv_dom_ll_catchatage = array(10, dim = dim(data$obs_srv_dom_ll_catchatage))
  data$srv_dom_ll_catchatage_comp_likelihood = 1
  test_model <- TMB::MakeADFun(data = data,
                               parameters = parameters,
                               DLL = "SpatialSablefishAssessment_TMBExports", silent  = T)

  test_report = test_model$report()

  F_pred_age = test_report$natage_f[,1,1] * test_report$F_fixed_f[,1,1] / test_report$Z_f[,1,1] * (1 - test_report$S_f[,1,1])
  M_pred_age = test_report$natage_m[,1,1] * test_report$F_fixed_m[,1,1] / test_report$Z_m[,1,1] * (1 - test_report$S_m[,1,1])

  ## check catch at age
  expect_equal(M_pred_age, test_report$catchatage_fixed_m[,1,1], tolerance = 0.001)
  expect_equal(F_pred_age, test_report$catchatage_fixed_f[,1,1], tolerance = 0.001)


  ## likelihood
  expected_nll = 0
  for(y in 1:(length(data$years))) {
    for(r in 1:data$n_regions) {
      expected_nll = expected_nll - sum(ddirichmult(obs = test_report$obs_fixed_catchatage[,r,y] / sum(test_report$obs_fixed_catchatage[,r,y]), est = test_report$pred_fixed_catchatage[,r,y], beta = sum(test_report$obs_fixed_catchatage[,r,y]) * test_report$theta_fixed_catchatage, n = sum(test_report$obs_fixed_catchatage[,r,y]), log_it = T))
    }
  }
  expect_equal(test_report$nll[1], expected_nll, tolerance = 0.1) ## with one log likelihood value

  ## Survey likelihood
  expected_nll = 0
  for(y in 1:(length(data$years))) {
    for(r in 1:data$n_regions) {
      expected_nll = expected_nll - sum(ddirichmult(obs = test_report$obs_srv_dom_ll_catchatage[,r,y]/ sum(test_report$obs_srv_dom_ll_catchatage[,r,y]), est = test_report$pred_srv_dom_ll_catchatage[,r,y], beta = sum(test_report$obs_srv_dom_ll_catchatage[,r,y]) * test_report$theta_srv_dom_ll_catchatage, n = sum(test_report$obs_srv_dom_ll_catchatage[,r,y]), log_it = T))
    }
  }
  expect_equal(test_report$nll[4], expected_nll, tolerance = 0.1) ## with one log likelihood value


  ## small theta
  parameters$trans_fixed_catchatage_error = log(0.1)
  parameters$trans_srv_dom_ll_catchatage = log(0.2)
  test_model <- TMB::MakeADFun(data = data,
                               parameters = parameters,
                               DLL = "SpatialSablefishAssessment_TMBExports", silent  = T)

  test_report = test_model$report()
  ## likelihood
  expected_nll = 0
  for(y in 1:(length(data$years))) {
    for(r in 1:data$n_regions) {
      expected_nll = expected_nll - sum(ddirichmult(obs = test_report$obs_fixed_catchatage[,r,y] / sum(test_report$obs_fixed_catchatage[,r,y]), est = test_report$pred_fixed_catchatage[,r,y], beta = sum(test_report$obs_fixed_catchatage[,r,y]) * test_report$theta_fixed_catchatage, n = sum(test_report$obs_fixed_catchatage[,r,y]), log_it = T))
    }
  }
  expect_equal(test_report$nll[1], expected_nll, tolerance = 0.1) ## with one log likelihood value
  ## Survey likelihood
  expected_nll = 0
  for(y in 1:(length(data$years))) {
    for(r in 1:data$n_regions) {
      expected_nll = expected_nll - sum(ddirichmult(obs = test_report$obs_srv_dom_ll_catchatage[,r,y]/ sum(test_report$obs_srv_dom_ll_catchatage[,r,y]), est = test_report$pred_srv_dom_ll_catchatage[,r,y], beta = sum(test_report$obs_srv_dom_ll_catchatage[,r,y]) * test_report$theta_srv_dom_ll_catchatage, n = sum(test_report$obs_srv_dom_ll_catchatage[,r,y]), log_it = T))
    }
  }
  expect_equal(test_report$nll[4], expected_nll, tolerance = 0.1) ## with one log likelihood value

  ## large Theta
  parameters$trans_fixed_catchatage_error = log(10)
  parameters$trans_srv_dom_ll_catchatage = log(7.4)

  test_model <- TMB::MakeADFun(data = data,
                               parameters = parameters,
                               DLL = "SpatialSablefishAssessment_TMBExports", silent  = T)

  test_report = test_model$report()
  ## likelihood
  expected_nll = 0
  for(y in 1:(length(data$years))) {
    for(r in 1:data$n_regions) {
      expected_nll = expected_nll - sum(ddirichmult(obs = test_report$obs_fixed_catchatage[,r,y] / sum(test_report$obs_fixed_catchatage[,r,y]), est = test_report$pred_fixed_catchatage[,r,y], beta = sum(test_report$obs_fixed_catchatage[,r,y]) * test_report$theta_fixed_catchatage, n = sum(test_report$obs_fixed_catchatage[,r,y]), log_it = T))
    }
  }
  expect_equal(test_report$nll[1], expected_nll, tolerance = 0.1) ## with one log likelihood value
  ## Survey likelihood
  expected_nll = 0
  for(y in 1:(length(data$years))) {
    for(r in 1:data$n_regions) {
      expected_nll = expected_nll - sum(ddirichmult(obs = test_report$obs_srv_dom_ll_catchatage[,r,y]/ sum(test_report$obs_srv_dom_ll_catchatage[,r,y]), est = test_report$pred_srv_dom_ll_catchatage[,r,y], beta = sum(test_report$obs_srv_dom_ll_catchatage[,r,y]) * test_report$theta_srv_dom_ll_catchatage, n = sum(test_report$obs_srv_dom_ll_catchatage[,r,y]), log_it = T))
    }
  }
  expect_equal(test_report$nll[4], expected_nll, tolerance = 0.1) ## with one log likelihood value

})

