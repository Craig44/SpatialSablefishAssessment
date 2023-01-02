#' test-fix-pars
#' @description tests fix pars function
#'
test_that("test-fix-pars", {
  ## add an additional time-block
  ## test fix-pars on 3D arrays
  data$srv_dom_ll_sel_type = rep(0,2)
  data$srv_dom_ll_sel_by_year_indicator[5:length(data$srv_dom_ll_sel_by_year_indicator)] = 1
  parameters$ln_srv_dom_ll_sel_pars = aperm(array(parameters$ln_srv_dom_ll_sel_pars, dim = c(2,2,2)))
  vals_to_exclude = cbind(1:dim(parameters$ln_srv_dom_ll_sel_pars)[1], 1,2)
  na_map = fix_pars(par_list = parameters, pars_to_exclude = c("ln_srv_dom_ll_sel_pars"),
                    array_elements_to_exclude = list(ln_srv_dom_ll_sel_pars = vals_to_exclude))

  ## check NA's are put in correct location
  expect_true(all(is.na(na_map$ln_srv_dom_ll_sel_pars[5:6])))
  expect_true(all(!is.na(na_map$ln_srv_dom_ll_sel_pars[-c(5:6)])))

  vals_to_exclude = cbind(1:dim(parameters$ln_srv_dom_ll_sel_pars)[1], 2,2)
  na_map = fix_pars(par_list = parameters, pars_to_exclude = c("ln_srv_dom_ll_sel_pars"),
                    array_elements_to_exclude = list(ln_srv_dom_ll_sel_pars = vals_to_exclude))
  expect_true(all(is.na(na_map$ln_srv_dom_ll_sel_pars[7:8])))
  expect_true(all(!is.na(na_map$ln_srv_dom_ll_sel_pars[-c(7:8)])))

  ## the first two slope params
  vals_to_exclude = cbind(1:dim(parameters$ln_srv_dom_ll_sel_pars)[1], 1,1)
  na_map = fix_pars(par_list = parameters, pars_to_exclude = c("ln_srv_dom_ll_sel_pars"),
                    array_elements_to_exclude = list(ln_srv_dom_ll_sel_pars = vals_to_exclude))
  expect_true(all(is.na(na_map$ln_srv_dom_ll_sel_pars[1:2])))
  expect_true(all(!is.na(na_map$ln_srv_dom_ll_sel_pars[-c(1:2)])))
  ## the first two slope params
  vals_to_exclude = cbind(1:dim(parameters$ln_srv_dom_ll_sel_pars)[1], 2,1)
  na_map = fix_pars(par_list = parameters, pars_to_exclude = c("ln_srv_dom_ll_sel_pars"),
                    array_elements_to_exclude = list(ln_srv_dom_ll_sel_pars = vals_to_exclude))
  expect_true(all(is.na(na_map$ln_srv_dom_ll_sel_pars[3:4])))
  expect_true(all(!is.na(na_map$ln_srv_dom_ll_sel_pars[-c(3:4)])))

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

  start_values = exp(parameters$ln_srv_dom_ll_sel_pars)

  test_model <- TMB::MakeADFun(data = data,
                               parameters = parameters, map = na_map,
                               DLL = "SpatialSablefishAssessment_TMBExports", silent  = T)


  ## there should be 3 srv selectivity pars as we have fixed one
  srv_sel_pars = test_model$par[names(test_model$par) %in% "ln_srv_dom_ll_sel_pars"]
  expect_true(length(srv_sel_pars) == 3)
  ## check
  test_report = test_model$report()
  expect_true(test_report$srv_dom_ll_sel_pars[1,1,1] == test_report$srv_dom_ll_sel_pars[1,1,2])

  # alternative params
  new_survey_slope_par = 3
  test_pars = test_model$par
  test_pars[which(names(test_model$par) %in% "ln_srv_dom_ll_sel_pars")[1]] = log(new_survey_slope_par)
  new_test_rep = test_model$report(test_pars)
  expect_true(new_test_rep$srv_dom_ll_sel_pars[1,1,1] == new_survey_slope_par)
  expect_true(new_test_rep$srv_dom_ll_sel_pars[1,1,1] == new_test_rep$srv_dom_ll_sel_pars[1,1,2])


  ## fix survey sel slope and delta parameters
  na_map = set_up_parameters(data= data, parameters = parameters, srv_sel_first_param_shared_by_sex = T, srv_sel_second_param_shared_by_sex = T)
  test_model <- TMB::MakeADFun(data = data,
                               parameters = parameters, map = na_map,
                               DLL = "SpatialSablefishAssessment_TMBExports", silent  = T)

  test_report = test_model$report()
  ## there should be 2 srv selectivity pars as we have fixed male and female
  srv_sel_pars = test_model$par[names(test_model$par) %in% "ln_srv_dom_ll_sel_pars"]
  expect_true(length(srv_sel_pars) == 2)

  expect_true(test_report$srv_dom_ll_sel_pars[1,1,1] == test_report$srv_dom_ll_sel_pars[1,1,2])
  expect_true(test_report$srv_dom_ll_sel_pars[1,2,1] == test_report$srv_dom_ll_sel_pars[1,2,2])
  ## the selectivities should have identical values
  expect_true(all(test_report$sel_srv_dom_ll_f == test_report$sel_srv_dom_ll_m))


  ## add an additional time-block
  data$srv_dom_ll_sel_type = rep(0,2)
  data$srv_dom_ll_sel_by_year_indicator[5:length(data$srv_dom_ll_sel_by_year_indicator)] = 1
  parameters$ln_srv_dom_ll_sel_pars = aperm(array(parameters$ln_srv_dom_ll_sel_pars, dim = c(2,2,2)))
  ## fix survey sel slope parameters for both time-blocks
  na_map = set_up_parameters(data= data, parameters = parameters, srv_sel_first_param_shared_by_sex = T)
  test_model <- TMB::MakeADFun(data = data,
                               parameters = parameters, map = na_map,
                               DLL = "SpatialSablefishAssessment_TMBExports", silent  = T)
  test_report = test_model$report()
  expect_true(length(test_model$par[names(test_model$par) %in% "ln_srv_dom_ll_sel_pars"]) == 6)

  expect_true(test_report$srv_dom_ll_sel_pars[1,1,1] == test_report$srv_dom_ll_sel_pars[1,1,2])
  expect_true(test_report$srv_dom_ll_sel_pars[2,1,1] == test_report$srv_dom_ll_sel_pars[2,1,2])

  ## FIX slope and delta
  na_map = set_up_parameters(data= data, parameters = parameters, srv_sel_first_param_shared_by_sex = T, srv_sel_second_param_shared_by_sex = T)

  test_model <- TMB::MakeADFun(data = data,
                               parameters = parameters, map = na_map,
                               DLL = "SpatialSablefishAssessment_TMBExports", silent  = T)


  ## there should be 3 srv selectivity pars as we have fixed one
  srv_sel_pars = test_model$par[names(test_model$par) %in% "ln_srv_dom_ll_sel_pars"]
  expect_true(length(srv_sel_pars) == 4)

  pars = test_model$par
  pars[names(test_model$par) %in% "ln_srv_dom_ll_sel_pars"] = log(c(1.3,1.5,0.5,0.6))
  test_report = test_model$report(pars)


  expect_true(test_report$srv_dom_ll_sel_pars[1,1,1] == test_report$srv_dom_ll_sel_pars[1,1,2])
  expect_true(test_report$srv_dom_ll_sel_pars[1,2,1] == test_report$srv_dom_ll_sel_pars[1,2,2])
  expect_true(test_report$srv_dom_ll_sel_pars[2,1,1] == test_report$srv_dom_ll_sel_pars[2,1,2])
  expect_true(test_report$srv_dom_ll_sel_pars[2,2,1] == test_report$srv_dom_ll_sel_pars[2,2,2])

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
  expect_true(length(logis_tag_pars) == sum(data$tag_recovery_indicator))

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
  data$tag_recovery_indicator = rep(0, length(data$tag_recovery_indicator))
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


})

