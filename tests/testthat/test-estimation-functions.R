#' test-fix-pars
#' @description tests fix pars function
#'
test_that("test-fix-pars", {
  load(system.file("testdata", "MockSablefishModel.RData",package="SpatialSablefishAssessment"))

  ## add an additional time-block
  ## test fix-pars on 3D arrays
  data$srv_sel_type = matrix(0,nrow = 2, ncol = data$n_surveys)
  data$srv_sel_by_year_indicator[5:length(data$srv_sel_by_year_indicator[,1]),1] = 1
  parameters$ln_srv_sel_pars = aperm(array(parameters$ln_srv_sel_pars, dim = c(2,2,2,2)))
  vals_to_exclude = cbind(1:dim(parameters$ln_srv_sel_pars)[1], 1,2,1)
  na_map = fix_pars(par_list = parameters, pars_to_exclude = c("ln_srv_sel_pars"),
                    array_elements_to_exclude = list(ln_srv_sel_pars = vals_to_exclude))

  ## check NA's are put in correct location
  expect_true(all(is.na(na_map$ln_srv_sel_pars[5:6])))
  expect_true(all(!is.na(na_map$ln_srv_sel_pars[-c(5:6)])))

  vals_to_exclude = cbind(1:dim(parameters$ln_srv_sel_pars)[1], 2,2, 1)
  na_map = fix_pars(par_list = parameters, pars_to_exclude = c("ln_srv_sel_pars"),
                    array_elements_to_exclude = list(ln_srv_sel_pars = vals_to_exclude))
  expect_true(all(is.na(na_map$ln_srv_sel_pars[7:8])))
  expect_true(all(!is.na(na_map$ln_srv_sel_pars[-c(7:8)])))

  ## the first two slope params
  vals_to_exclude = cbind(1:dim(parameters$ln_srv_sel_pars)[1], 1,1, 1)
  na_map = fix_pars(par_list = parameters, pars_to_exclude = c("ln_srv_sel_pars"),
                    array_elements_to_exclude = list(ln_srv_sel_pars = vals_to_exclude))
  expect_true(all(is.na(na_map$ln_srv_sel_pars[1:2])))
  expect_true(all(!is.na(na_map$ln_srv_sel_pars[-c(1:2)])))
  ## the first two slope params
  vals_to_exclude = cbind(1:dim(parameters$ln_srv_sel_pars)[1], 2,1, 1)
  na_map = fix_pars(par_list = parameters, pars_to_exclude = c("ln_srv_sel_pars"),
                    array_elements_to_exclude = list(ln_srv_sel_pars = vals_to_exclude))
  expect_true(all(is.na(na_map$ln_srv_sel_pars[3:4])))
  expect_true(all(!is.na(na_map$ln_srv_sel_pars[-c(3:4)])))

  ## test fix-pars on vectors
  data$n_init_rec_devs = 5
  parameters$ln_init_rec_dev = rep(log(2), data$n_init_rec_devs)
  na_map = fix_pars(par_list = parameters, pars_to_exclude = c("ln_init_rec_dev"),
                    vec_elements_to_exclude = list(ln_init_rec_dev = c(1,2)))
  expect_true(all(is.na(na_map$ln_init_rec_dev[1:2])))
  expect_true(all(!is.na(na_map$ln_init_rec_dev[-c(1:2)])))

})

#' test-shared-survey-selectivity-pars
#' @description tests that set_up_parameters works as expected
#'
test_that("test-shared-survey-selectivity-pars", {
  load(system.file("testdata", "MockSablefishModel.RData",package="SpatialSablefishAssessment"))
  data$F_method = 1
  ## fix survey sel slope parameters
  na_map = set_up_parameters(data= data, parameters = parameters, srv_sel_first_param_shared_by_sex = T)

  start_values = exp(parameters$ln_srv_sel_pars)

  test_model <- TMB::MakeADFun(data = data,
                               parameters = parameters, map = na_map,
                               DLL = "SpatialSablefishAssessment_TMBExports", silent  = T)


  ## there should be 3 srv selectivity pars as we have fixed one
  srv_sel_pars = test_model$par[names(test_model$par) %in% "ln_srv_sel_pars"]
  expect_true(length(srv_sel_pars) == 3)
  ## check
  test_report = test_model$report()
  expect_true(test_report$srv_sel_pars[1,1,1,1] == test_report$srv_sel_pars[1,1,2,1])

  # alternative params
  new_survey_slope_par = 3
  test_pars = test_model$par
  test_pars[which(names(test_model$par) %in% "ln_srv_sel_pars")[1]] = log(new_survey_slope_par)
  new_test_rep = test_model$report(test_pars)
  expect_equal(new_test_rep$srv_sel_pars[1,1,1,1], new_survey_slope_par, tolerance = 0.00001)
  expect_equal(new_test_rep$srv_sel_pars[1,1,1,1], new_test_rep$srv_sel_pars[1,1,2, 1], tolerance = 0.00001)


  ## fix survey sel slope and delta parameters
  na_map = set_up_parameters(data= data, parameters = parameters, srv_sel_first_param_shared_by_sex = T, srv_sel_second_param_shared_by_sex = T)
  test_model <- TMB::MakeADFun(data = data,
                               parameters = parameters, map = na_map,
                               DLL = "SpatialSablefishAssessment_TMBExports", silent  = T)

  test_report = test_model$report()
  ## there should be 2 srv selectivity pars as we have fixed male and female
  srv_sel_pars = test_model$par[names(test_model$par) %in% "ln_srv_sel_pars"]
  expect_true(length(srv_sel_pars) == 2)

  expect_equal(test_report$srv_sel_pars[1,1,1,1], test_report$srv_sel_pars[1,1,2,1], tolerance = 0.00001)
  expect_equal(test_report$srv_sel_pars[1,2,1,1], test_report$srv_sel_pars[1,2,2,1], tolerance = 0.00001)
  ## the selectivities should have identical values
  expect_true(all(test_report$sel_srv_f ==  test_report$sel_srv_m))


  ## add an additional time-block
  data$srv_sel_type = matrix(0, nrow = 2, ncol = data$n_surveys)
  data$srv_sel_by_year_indicator[5:length(data$srv_sel_by_year_indicator)] = 1
  parameters$ln_srv_sel_pars = abind::abind(parameters$ln_srv_sel_pars, parameters$ln_srv_sel_pars, along = 1)


  ## fix survey sel slope parameters for both time-blocks
  na_map = set_up_parameters(data = data, parameters = parameters, srv_sel_first_param_shared_by_sex = T)
  test_model <- TMB::MakeADFun(data = data,
                               parameters = parameters, map = na_map,
                               DLL = "SpatialSablefishAssessment_TMBExports", silent  = T)
  test_report = test_model$report()
  expect_true(length(test_model$par[names(test_model$par) %in% "ln_srv_sel_pars"]) == 6)

  expect_equal(test_report$srv_sel_pars[1,1,1,1], test_report$srv_sel_pars[1,1,2,1], tolerance = 0.00001)
  expect_equal(test_report$srv_sel_pars[2,1,1,1], test_report$srv_sel_pars[2,1,2,1], tolerance = 0.00001)
  ## FIX slope and delta
  na_map = set_up_parameters(data= data, parameters = parameters, srv_sel_first_param_shared_by_sex = T, srv_sel_second_param_shared_by_sex = T)

  test_model <- TMB::MakeADFun(data = data,
                               parameters = parameters, map = na_map,
                               DLL = "SpatialSablefishAssessment_TMBExports", silent  = T)


  ## there should be 3 srv selectivity pars as we have fixed one
  srv_sel_pars = test_model$par[names(test_model$par) %in% "ln_srv_sel_pars"]
  expect_true(length(srv_sel_pars) == 4)

  pars = test_model$par
  pars[names(test_model$par) %in% "ln_srv_sel_pars"] = log(c(1.3,1.5,0.5,0.6))
  test_report = test_model$report(pars)


  expect_equal(test_report$srv__sel_pars[1,1,1,1], test_report$srv__sel_pars[1,1,2,1], tolerance = 0.00001)
  expect_equal(test_report$srv__sel_pars[1,2,1,1], test_report$srv__sel_pars[1,2,2,1], tolerance = 0.00001)
  expect_equal(test_report$srv__sel_pars[2,1,1,1], test_report$srv__sel_pars[2,1,2,1], tolerance = 0.00001)
  expect_equal(test_report$srv__sel_pars[2,2,1,1], test_report$srv__sel_pars[2,2,2,1], tolerance = 0.00001)

  #######
  ## repeat for fixed gear fishery
  #######
  load(system.file("testdata", "MockSablefishModel.RData",package="SpatialSablefishAssessment"))
  data$F_method = 1
  ## fix survey sel slope parameters
  na_map = set_up_parameters(data= data, parameters = parameters, fixed_sel_first_shared_by_sex = T, fixed_sel_second_shared_by_sex = T)
  test_model <- TMB::MakeADFun(data = data,
                               parameters = parameters, map = na_map,
                               DLL = "SpatialSablefishAssessment_TMBExports", silent  = T)

  ## there should be 3 srv selectivity pars as we have fixed one
  srv_sel_pars = test_model$par[names(test_model$par) %in% "ln_fixed_sel_pars"]
  expect_true(length(srv_sel_pars) == 2)
  ## check
  test_report = test_model$report()
  expect_true(all(test_report$sel_fixed_f == test_report$sel_fixed_m))

  #######
  ## repeat for trawl gear fishery
  #######
  load(system.file("testdata", "MockSablefishModel.RData",package="SpatialSablefishAssessment"))
  data$F_method = 1
  ## fix survey sel slope parameters
  na_map = set_up_parameters(data= data, parameters = parameters, trwl_sel_first_shared_by_sex = T, trwl_sel_second_shared_by_sex = T)
  test_model <- TMB::MakeADFun(data = data,
                               parameters = parameters, map = na_map,
                               DLL = "SpatialSablefishAssessment_TMBExports", silent  = T)

  ## there should be 3 srv selectivity pars as we have fixed one
  srv_sel_pars = test_model$par[names(test_model$par) %in% "ln_trwl_sel_pars"]
  expect_true(length(srv_sel_pars) == 2)
  ## check
  test_report = test_model$report()
  expect_true(all(test_report$sel_trwl_m == test_report$sel_trwl_f))

  #######
  ## tag reporting rates
  #######
  load(system.file("testdata", "MockSablefishModel.RData",package="SpatialSablefishAssessment"))
  data$F_method = 1
  ## fix survey sel slope parameters
  na_map = set_up_parameters(data= data, parameters = parameters, tag_reporting_rate = "constant")
  test_model <- TMB::MakeADFun(data = data,
                               parameters = parameters, map = na_map,
                               DLL = "SpatialSablefishAssessment_TMBExports", silent  = T)

  ## there should be 3 srv selectivity pars as we have fixed one
  logis_tag_pars = test_model$par[names(test_model$par) %in% "logistic_tag_reporting_rate"]
  expect_true(length(logis_tag_pars) == 1)
  test_report = test_model$report()
  pars =  test_model$par
  pars[names(pars) %in% "logistic_tag_reporting_rate"] = qlogis(0.5)
  test_rep = test_model$report(pars)
  expect_true(all(test_rep$tag_reporting_rate == 0.5))

  ## annual
  na_map = set_up_parameters(data= data, parameters = parameters, tag_reporting_rate = "time")
  test_model <- TMB::MakeADFun(data = data,
                               parameters = parameters, map = na_map,
                               DLL = "SpatialSablefishAssessment_TMBExports", silent  = T)

  ## there should be 3 srv selectivity pars as we have fixed one
  logis_tag_pars = test_model$par[names(test_model$par) %in% "logistic_tag_reporting_rate"]
  expect_true(length(logis_tag_pars) == sum(data$tag_recovery_indicator_by_year))

  pars =  test_model$par
  report_vals = seq(from = 0.1, to = 0.9, length = sum(names(pars) %in% "logistic_tag_reporting_rate"))
  pars[names(pars) %in% "logistic_tag_reporting_rate"] = qlogis(report_vals)

  test_rep = test_model$report(pars)
  ## check all regions are the same
  for(i in 1:length(report_vals)) {
    expect_equal(test_rep$tag_reporting_rate[,i], rep(report_vals[i], nrow(test_rep$tag_reporting_rate)), tolerance = 0.001)
  }

  ## check if there are no tag recovery observations that tag related parameters are turned off.
  load(system.file("testdata", "MockSablefishModel.RData",package="SpatialSablefishAssessment"))
  data$tag_recovery_indicator = rep(0, length(data$tag_recovery_indicator_by_year))
  na_map = set_up_parameters(data= data, parameters = parameters, tag_reporting_rate = "constant")
  expect_true(is.na(na_map$ln_tag_phi))
  expect_true(all(is.na(na_map$logistic_tag_reporting_rate)))

  ## check if F_method == 1 then don't estimate F parameters
  load(system.file("testdata", "MockSablefishModel.RData",package="SpatialSablefishAssessment"))
  data$F_method = 1
  ## fix survey sel slope parameters
  na_map = set_up_parameters(data= data, parameters = parameters)
  expect_true(all(is.na(na_map$ln_fixed_F_devs)))
  expect_true(all(is.na(na_map$ln_trwl_F_devs)))
  expect_true(all(is.na(na_map$ln_fixed_F_avg)))
  expect_true(all(is.na(na_map$ln_trwl_F_avg)))



  ## tag-reporting-space
  load(system.file("testdata", "MockSablefishModel.RData",package="SpatialSablefishAssessment"))
  data$F_method = 1
  ## fix survey sel slope parameters
  na_map = set_up_parameters(data= data, parameters = parameters, srv_sel_first_param_shared_by_sex = T, tag_reporting_rate = "space")

  test_model <- TMB::MakeADFun(data = data,
                               parameters = parameters, map = na_map,
                               DLL = "SpatialSablefishAssessment_TMBExports", silent  = T)


  ##
  tag_report_pars = test_model$par[names(test_model$par) %in% "logistic_tag_reporting_rate"]
  expect_true(length(tag_report_pars) == data$n_regions)
  ## check
  test_report = test_model$report()
  expect_true(test_report$srv_sel_pars[1,1,1, 1] == test_report$srv_sel_pars[1,1,2, 1])

  # alternative params
  new_tag_report_pars = (c(0.1,0.2,0.3,0.4,0.5))
  test_pars = test_model$par
  test_pars[which(names(test_model$par) %in% "logistic_tag_reporting_rate")] = logit(new_tag_report_pars)
  new_test_rep = test_model$report(test_pars)

  for(r in 1:data$n_regions)
    expect_equal(new_test_rep$tag_reporting_rate[r,], rep(new_tag_report_pars[r],ncol(new_test_rep$tag_reporting_rate)), tolerance = 0.0001)


  ## proportion male recruits
  load(system.file("testdata", "MockSablefishModel.RData",package="SpatialSablefishAssessment"))
  data$F_method = 1
  ## fix survey sel slope parameters
  na_map = set_up_parameters(data= data, parameters = parameters, srv_sel_first_param_shared_by_sex = T,
                             tag_reporting_rate = "space",
                             est_prop_male_recruit = "constant")

  test_model <- TMB::MakeADFun(data = data,
                               parameters = parameters, map = na_map,
                               DLL = "SpatialSablefishAssessment_TMBExports", silent  = T)


  ##
  prop_male_pars = test_model$par[names(test_model$par) %in% "logistic_prop_recruit_male"]
  expect_true(length(prop_male_pars) == 1)
  ## check
  test_report = test_model$report()
  expect_true(all(test_report$prop_recruit_male == test_report$prop_recruit_female))

  # alternative params
  new_prop_recruit_pars = 0.3
  test_pars = test_model$par
  test_pars[which(names(test_model$par) %in% "logistic_prop_recruit_male")] = logit(new_prop_recruit_pars)
  new_test_rep = test_model$report(test_pars)

  expect_equal(new_test_rep$prop_recruit_male, rep(new_prop_recruit_pars, length(data$years)), tolerance = 0.0001)

  ## proportion male recruits time-blocks
  prop_time_blocks = c(2011, 2014, 2018)
  load(system.file("testdata", "MockSablefishModel.RData",package="SpatialSablefishAssessment"))
  data$F_method = 1
  ## fix survey sel slope parameters
  na_map = set_up_parameters(data= data, parameters = parameters, srv_sel_first_param_shared_by_sex = T,
                             tag_reporting_rate = "space",
                             est_prop_male_recruit = prop_time_blocks)
  test_model <- TMB::MakeADFun(data = data,
                               parameters = parameters, map = na_map,
                               DLL = "SpatialSablefishAssessment_TMBExports", silent  = T)


  ##
  prop_male_pars = test_model$par[names(test_model$par) %in% "logistic_prop_recruit_male"]
  expect_true(length(prop_male_pars) == length(prop_time_blocks))

  # alternative params
  new_prop_recruit_pars = c(0.3, 0.6,0.8)
  test_pars = test_model$par
  test_pars[which(names(test_model$par) %in% "logistic_prop_recruit_male")] = logit(new_prop_recruit_pars)
  new_test_rep = test_model$report(test_pars)
  expect_equal(new_test_rep$prop_recruit_male, c(0.5, rep(0.3, 3), rep(0.6, 4), rep(0.8, 3)), tolerance = 0.0001)

  ## Constant over space - Q
  load(system.file("testdata", "MockSablefishModel.RData",package="SpatialSablefishAssessment"))
  data$F_method = 1
  data$srv_q_by_year_indicator[5:11] = 1
  parameters$trans_srv_q = array(0, dim = c(data$n_regions,length(unique(data$srv_q_by_year_indicator)), data$n_surveys))
  ## fix survey sel slope parameters
  na_map = set_up_parameters(data= data, parameters = parameters, srv_sel_first_param_shared_by_sex = T,
                             srv_q_spatial = F)
  test_model <- TMB::MakeADFun(data = data,
                               parameters = parameters, map = na_map,
                               DLL = "SpatialSablefishAssessment_TMBExports", silent  = T)

  ##
  trans_q = test_model$par[names(test_model$par) %in% "trans_srv_q"]
  expect_true(length(trans_q) == 2)

  # alternative params
  new_q = c(0.3, 0.6)
  test_pars = test_model$par
  test_pars[which(names(test_model$par) %in% "trans_srv_q")] = logit(new_q)
  new_test_rep = test_model$report(test_pars)
  for(i in 1:length(new_q))
    expect_equal(new_test_rep$srv_q[,i,1], rep(new_q[i], data$n_regions), tolerance = 0.0001)
})



#' test-fix_pars_for_multiple_surveys
#' @description test the setup function with multiple surveys
#'
test_that("test-shared-survey-selectivity-pars", {
  load(system.file("testdata", "MockSablefishModel.RData",package="SpatialSablefishAssessment"))
  ## 2 surveys 1 time-block
  data$n_surveys = 2
  data$srv_bio_indicator = array(1, dim = c(data$n_regions, length(data$years), data$n_surveys))
  data$srv_catchatage_indicator = array(1, dim = c(data$n_regions, length(data$years), data$n_surveys))
  data$obs_srv_bio =  array(50, dim = c(data$n_regions, length(data$years), data$n_surveys))
  data$obs_srv_se =  array(0.3, dim = c(data$n_regions, length(data$years), data$n_surveys))
  data$obs_srv_catchatage = array(5, dim = c(2*length(data$ages),data$n_regions, length(data$years), data$n_surveys))
  data$srv_sel_by_year_indicator = array(0, dim = c(length(data$years), data$n_surveys))
  data$srv_sel_type = array(0, dim = c(1, data$n_surveys))

  data$srv_q_by_year_indicator = array(0, dim = c(length(data$years), data$n_surveys))
  data$q_is_nuisance = rep(0, data$n_surveys)
  data$srv_obs_is_abundance = rep(1, data$n_surveys)

  data$srv_q_transformation = rep(1, data$n_surveys)
  data$srv_catchatage_comp_likelihood = rep(0, data$n_surveys)
  # change parameters
  parameters$trans_srv_q = abind::abind(parameters$trans_srv_q,  parameters$trans_srv_q, along = 3)
  parameters$ln_srv_sel_pars = abind::abind(parameters$ln_srv_sel_pars, parameters$ln_srv_sel_pars, along = 4)

  #parameters$srv_sel_by_year_indicator = array(0, dim = c(length(data$years), data$n_surveys))
  expect_true(validate_input_data_and_parameters(data, parameters))
  test_model <- TMB::MakeADFun(data = data,
                               parameters = parameters,
                               DLL = "SpatialSablefishAssessment_TMBExports", silent  = T)

  test_report = test_model$report()

  expect_equal(test_report$pred_srv_bio[,,1], test_report$pred_srv_bio[,,2], tolerance = 0.0001)

  ## different Q's among surveys same selectiviies
  data$srv_q_transformation = rep(0, data$n_surveys)
  ## survey 1 Q = 1, survey 2 Q = 0.5
  parameters$trans_srv_q = abind::abind(array(log(1), dim = c(data$n_regions, 1,1)),  array(log(0.5), dim = c(data$n_regions, 1,1)), along = 3)
  test_model <- TMB::MakeADFun(data = data,
                               parameters = parameters,
                               DLL = "SpatialSablefishAssessment_TMBExports", silent  = T)

  test_report = test_model$report()
  ## Survey 2 should be half the size of different Q's
  expect_equal(test_report$pred_srv_bio[,,1] * 0.5, test_report$pred_srv_bio[,,2], tolerance = 0.0001)

  ## 2 surveys, Survey 1 has two time-blocks for Q and selectivity and Survey 2 has one time-block for Q and selectivity
  ## Survey 1 has years 6-11 with different selectivities and Q's
  data$srv_sel_by_year_indicator[6:11, 1] = 1
  data$srv_q_by_year_indicator[6:11, 1] = 1

  ## should be false due to the following message
  expect_equal(validate_input_data_and_parameters(data, parameters), "srv_sel_type: different number of elements for dimension 1. Array had 1 we expected 2")
  ##
  data$srv_sel_type = array(0, dim = c(2, data$n_surveys))
  parameters$trans_srv_q = array(log(0.2), dim = c(data$n_regions, 2, data$n_surveys))
  parameters$ln_srv_sel_pars = abind::abind(array(c(2.111, -0.711,1.576, -0.7118), dim = c(1,2, 2, 1)), array(c(3.111, -0.511,1.976, -0.9118), dim = c(1,2, 2, 1)), along = 1)
  parameters$ln_srv_sel_pars = abind::abind(parameters$ln_srv_sel_pars, parameters$ln_srv_sel_pars, along = 4)

  ## run setup code and see if the Q's for the second Survey get turned off
  na_map = set_up_parameters(data = data, parameters = parameters)
  expect_true(length(levels(na_map$trans_srv_q)) == 3) ## should be 3 levels
  # the last 5 values should all be NA's
  expect_true(all(is.na(na_map$trans_srv_q[(length(na_map$trans_srv_q) - 4):length(na_map$trans_srv_q)]))) ## should be 3 levels
  # a q for each region and survey
  na_map = set_up_parameters(data = data, parameters = parameters, srv_q_spatial = T)
  expect_true(length(levels(na_map$trans_srv_q)) == data$n_regions * 3) ## should be 3 levels
  # again the last 5 values should all be NA's
  expect_true(all(is.na(na_map$trans_srv_q[(length(na_map$trans_srv_q) - 4):length(na_map$trans_srv_q)]))) ## should be 3 levels
  ## check the survey selectivity
  na_map = set_up_parameters(data = data, parameters = parameters, srv_sel_first_param_shared_by_sex = T)
  test_model <- TMB::MakeADFun(data = data, map = na_map,
                               parameters = parameters,
                               DLL = "SpatialSablefishAssessment_TMBExports", silent  = T)

  test_report = test_model$report()

  expect_equal(test_report$srv_sel_pars[1,1,1,1], test_report$srv_sel_pars[1,1,2,1], tolerance = 0.001)
  expect_equal(test_report$srv_sel_pars[2,1,1,1], test_report$srv_sel_pars[2,1,2,1], tolerance = 0.001)
  expect_equal(test_report$srv_sel_pars[1,1,1,2], test_report$srv_sel_pars[1,1,2,2], tolerance = 0.001)
  ## this survey has no second selectivity time-block so needs to be turned off.
  expect_true(test_report$srv_sel_pars[2,1,1,2] != test_report$srv_sel_pars[2,1,2,2])

})

#' test-get_negloglike
#' @description tests get_negloglike  function
#'
test_that("test-get_negloglike", {
  load(system.file("testdata", "MockSablefishModel.RData",package="SpatialSablefishAssessment"))
  test_model <- TMB::MakeADFun(data = data,
                               parameters = parameters,
                               DLL = "SpatialSablefishAssessment_TMBExports", silent  = T)

  mle_report = test_model$report()
  ## when we change the neg log likelihood check this utility function doesn't break
  nll = get_negloglike(mle_report)
  expect_true(nrow(nll) == 11) ## 11 slots
})


#' test-get_tmb_parameter_element
#' @description tests get_tmb_parameter_element  function
#'
test_that("test-get_tmb_parameter_element", {

  parameters = list(
    scalar_param = 2,
    vector_param = 1:3,
    array_2d = matrix(1:6, ncol = 3),
    array_3d = array(1:18, dim = c(2,3,6))
  )

  ## scalar
  expect_true(get_tmb_parameter_element(parameters = parameters, parameter_label = "scalar_param", element = 1) == 1)
  expect_error(get_tmb_parameter_element(parameters = parameters, parameter_label = "scalar_param", element = 2))
  ## vector
  expect_true(get_tmb_parameter_element(parameters = parameters, parameter_label = "vector_param", element = 1) == 1)
  expect_true(get_tmb_parameter_element(parameters = parameters, parameter_label = "vector_param", element = 2) == 2)
  expect_true(get_tmb_parameter_element(parameters = parameters, parameter_label = "vector_param", element = 3) == 3)
  expect_error(get_tmb_parameter_element(parameters = parameters, parameter_label = "vector_param", element = 4))
  ## 2D array
  expect_true(get_tmb_parameter_element(parameters = parameters, parameter_label = "array_2d", element = matrix(c(1,1), nrow = 1)) == 1)
  expect_true(get_tmb_parameter_element(parameters = parameters, parameter_label = "array_2d", element = matrix(c(2,1), nrow = 1)) == 2)
  expect_true(get_tmb_parameter_element(parameters = parameters, parameter_label = "array_2d", element = matrix(c(1,2), nrow = 1)) == 3)
  expect_true(get_tmb_parameter_element(parameters = parameters, parameter_label = "array_2d", element = matrix(c(2,2), nrow = 1)) == 4)
  expect_true(get_tmb_parameter_element(parameters = parameters, parameter_label = "array_2d", element = matrix(c(1,3), nrow = 1)) == 5)
  expect_true(get_tmb_parameter_element(parameters = parameters, parameter_label = "array_2d", element = matrix(c(2,3), nrow = 1)) == 6)
  ## 3D array
  expect_true(get_tmb_parameter_element(parameters = parameters, parameter_label = "array_3d", element = matrix(c(1,1, 1), nrow = 1)) == 1)
  expect_true(get_tmb_parameter_element(parameters = parameters, parameter_label = "array_3d", element = matrix(c(2,1, 1), nrow = 1)) == 2)
  expect_true(get_tmb_parameter_element(parameters = parameters, parameter_label = "array_3d", element = matrix(c(1,2, 1), nrow = 1)) == 3)
  expect_true(get_tmb_parameter_element(parameters = parameters, parameter_label = "array_3d", element = matrix(c(2,2, 1), nrow = 1)) == 4)
  expect_true(get_tmb_parameter_element(parameters = parameters, parameter_label = "array_3d", element = matrix(c(1,3, 1), nrow = 1)) == 5)
  expect_true(get_tmb_parameter_element(parameters = parameters, parameter_label = "array_3d", element = matrix(c(2,3, 1), nrow = 1)) == 6)

  expect_true(get_tmb_parameter_element(parameters = parameters, parameter_label = "array_3d", element = matrix(c(1,1, 2), nrow = 1)) == 7)
  expect_true(get_tmb_parameter_element(parameters = parameters, parameter_label = "array_3d", element = matrix(c(2,1, 2), nrow = 1)) == 8)
  expect_true(get_tmb_parameter_element(parameters = parameters, parameter_label = "array_3d", element = matrix(c(1,2, 2), nrow = 1)) == 9)
  expect_true(get_tmb_parameter_element(parameters = parameters, parameter_label = "array_3d", element = matrix(c(2,2, 2), nrow = 1)) == 10)
  expect_true(get_tmb_parameter_element(parameters = parameters, parameter_label = "array_3d", element = matrix(c(1,3, 2), nrow = 1)) ==11)
  expect_true(get_tmb_parameter_element(parameters = parameters, parameter_label = "array_3d", element = matrix(c(2,3, 2), nrow = 1)) == 12)

  expect_true(get_tmb_parameter_element(parameters = parameters, parameter_label = "array_3d", element = matrix(c(1,1, 3), nrow = 1)) == 13)
  expect_true(get_tmb_parameter_element(parameters = parameters, parameter_label = "array_3d", element = matrix(c(2,1, 3), nrow = 1)) == 14)
  expect_true(get_tmb_parameter_element(parameters = parameters, parameter_label = "array_3d", element = matrix(c(1,2, 3), nrow = 1)) == 15)
  expect_true(get_tmb_parameter_element(parameters = parameters, parameter_label = "array_3d", element = matrix(c(2,2, 3), nrow = 1)) == 16)
  expect_true(get_tmb_parameter_element(parameters = parameters, parameter_label = "array_3d", element = matrix(c(1,3, 3), nrow = 1)) == 17)
  expect_true(get_tmb_parameter_element(parameters = parameters, parameter_label = "array_3d", element = matrix(c(2,3, 3), nrow = 1)) == 18)

  expect_error(get_tmb_parameter_element(parameters = parameters, parameter_label = "scalar_param", element = matrix(c(2,3, 3), nrow = 1)))
  expect_error(get_tmb_parameter_element(parameters = parameters, parameter_label = "array_3d", element = 1))
})



#' test-profile_param
#' @description tests profile_param  function
#'
test_that("test-profile_param", {
  load(system.file("testdata", "MockSablefishModel.RData",package="SpatialSablefishAssessment"))

  ## start with fixing a recruitment dev
  data$F_method = 1
  na_map = set_up_parameters(data= data, parameters = parameters, trwl_sel_first_shared_by_sex = F, trwl_sel_second_shared_by_sex = F)

  test_model <- TMB::MakeADFun(data = data,
                               parameters = parameters, map = na_map,
                               DLL = "SpatialSablefishAssessment_TMBExports", silent  = T)

  test_rep = test_model$report()
  rec_vals = c(-1,0,1)
  rec_ndx = 11
  test_rec = profile_param(parameters, test_model, na_map, profile_param_label = "trans_rec_dev", element = matrix(c(1,rec_ndx),nrow = 1), profile_values = rec_vals,  same_element = NULL, no_estimation = T, verbose= F)
  for(i in 1:length(rec_vals)) {
    expect_equal(test_rec$profile_mle[[i]]$recruitment_multipliers[,rec_ndx],  rep(exp(rec_vals[i] - 0.5 * test_rep$sigma_R^2), nrow(test_rec$profile_mle[[i]]$recruitment_multipliers)), tolerance = 0.0001)
  }

  ## start by profiling a mean recruitment element
  mean_rec_profiles = c(13, 14, 15, 16)
  mean_rec_ndx = 2
  test_rec = profile_param(parameters, test_model, na_map, profile_param_label = "ln_mean_rec", element = mean_rec_ndx, profile_values = mean_rec_profiles,  same_element = NULL, no_estimation = T, verbose= F)
  for(i in 1:length(mean_rec_profiles)) {
    ## check the profiled value was set to the correct value
    expect_equal(test_rec$profile_mle[[i]]$mean_rec[mean_rec_ndx],  exp(mean_rec_profiles[i]), tolerance = 0.0001)
    ## check the others didn't change
    expect_equal(test_rec$profile_mle[[i]]$mean_rec[1],  exp(parameters$ln_mean_rec[1]), tolerance = 0.0001)

  }

  ## check the same element works
  ## profile catchability
  profile_q_vals = (c(0.1, 0.2, 0.3,0.4))
  q_ndx = 1
  test_profile = profile_param(parameters, test_model, na_map, profile_param_label = "trans_srv_q", element = matrix(c(q_ndx,1, 1),nrow = 1), profile_values = qlogis(profile_q_vals),  same_element = matrix(c(2,1,1,3,1,1, 4, 1,1, 5, 1,1), byrow = T, ncol = 3), no_estimation = T, verbose= F)
  for(i in 1:length(profile_q_vals)) {
    ## check the profiled value was set to the correct value
    expect_equal(as.numeric(test_profile$profile_mle[[i]]$srv_q),  rep(profile_q_vals[i], nrow(test_profile$profile_mle[[i]]$srv_q)), tolerance = 0.0001)
  }
})
