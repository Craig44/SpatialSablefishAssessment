

#' test-Selectivities
#' @description this is to test selectivities
#'
test_that("test-Selectivities", {
  ## Read in mock data
  load(system.file("testdata", "MockSablefishModel.RData",package="SpatialSablefishAssessment"))
  data$model = "TagIntegratedValidate"
  ## no Z or movement
  data$apply_Z_on_tagged_fish = 0
  ## two time-blocks
  data$trwl_sel_by_year_indicator[5:length(data$trwl_sel_by_year_indicator)] = 1
  data$trwl_sel_type = c(1,1)
  data$apply_fixed_movement = 1 ## no movement initial age-structure should be only a function of ageing and M
  parameters$ln_trwl_sel_pars = array(0, dim = c(2, 2, 2))
  parameters$ln_trwl_sel_pars[1,1,1] = 2.311
  parameters$ln_trwl_sel_pars[1,2,1] = 2.2150
  parameters$ln_trwl_sel_pars[1,1,2] = 2.011
  parameters$ln_trwl_sel_pars[1,2,2] = 2.1150
  parameters$ln_trwl_sel_pars[2,1,1] = 1.311
  parameters$ln_trwl_sel_pars[2,2,1] = 2.4150
  parameters$ln_trwl_sel_pars[2,1,2] = 1.011
  parameters$ln_trwl_sel_pars[2,2,2] = 2.21150

  ## srv with 3 time-block
  data$srv_dom_ll_sel_type = c(0,0,0)

  parameters$ln_srv_dom_ll_sel_pars  = array(0, dim = c(3, 2, 2))
  parameters$ln_srv_dom_ll_sel_pars[1,1,1] = 2.111
  parameters$ln_srv_dom_ll_sel_pars[1,2,1] = -0.711
  parameters$ln_srv_dom_ll_sel_pars[1,1,2] = 1.576
  parameters$ln_srv_dom_ll_sel_pars[1,2,2] = -0.7118
  parameters$ln_srv_dom_ll_sel_pars[2,1,1] = 2.111
  parameters$ln_srv_dom_ll_sel_pars[2,2,1] = -0.711
  parameters$ln_srv_dom_ll_sel_pars[2,1,2] = 1.576
  parameters$ln_srv_dom_ll_sel_pars[2,2,2] = -0.7118
  parameters$ln_srv_dom_ll_sel_pars[3,1,1] = 2.111
  parameters$ln_srv_dom_ll_sel_pars[3,2,1] = -0.711
  parameters$ln_srv_dom_ll_sel_pars[3,1,2] = 1.576
  parameters$ln_srv_dom_ll_sel_pars[3,2,2] = -0.7118


  test_model <- TMB::MakeADFun(data = data,
                               parameters = parameters,
                               DLL = "SpatialSablefishAssessment_TMBExports", silent  = T)

  test_report = test_model$report()

  ## double normal
  # Time-block 1
  # Female
  expect_equal(test_report$sel_trwl_f[,1], double_normal(x = data$ages, x50 = exp(parameters$ln_trwl_sel_pars[1,1,2]), delta = exp(parameters$ln_trwl_sel_pars[1,2,2])))
  # Male
  expect_equal(test_report$sel_trwl_m[,1], double_normal(x = data$ages, x50 = exp(parameters$ln_trwl_sel_pars[1,1,1]), delta = exp(parameters$ln_trwl_sel_pars[1,2,1])))
  # Time-block 2
  # Female
  expect_equal(test_report$sel_trwl_f[,2], double_normal(x = data$ages, x50 = exp(parameters$ln_trwl_sel_pars[2,1,2]), delta = exp(parameters$ln_trwl_sel_pars[2,2,2])))
  # Male
  expect_equal(test_report$sel_trwl_m[,2], double_normal(x = data$ages, x50 = exp(parameters$ln_trwl_sel_pars[2,1,1]), delta = exp(parameters$ln_trwl_sel_pars[2,2,1])))

  ## Logistic
  # Time-block 1
  # Female
  expect_equal(test_report$sel_fixed_f[,1], logis_alt(x = data$ages, x50 = exp(parameters$ln_fixed_sel_pars[1,1,2]), delta = exp(parameters$ln_fixed_sel_pars[1,2,2])))
  # Male
  expect_equal(test_report$sel_fixed_m[,1], logis_alt(x = data$ages, x50 = exp(parameters$ln_fixed_sel_pars[1,1,1]), delta = exp(parameters$ln_fixed_sel_pars[1,2,1])))

  ## test log and exp change
  expect_true(all(test_report$fixed_sel_pars == exp(parameters$ln_fixed_sel_pars)))
  expect_true(all(!(test_report$trwl_sel_pars - exp(parameters$ln_trwl_sel_pars)) > 0.000001))
  expect_true(all(!(test_report$srv_dom_ll_sel_pars - exp(parameters$ln_srv_dom_ll_sel_pars)) > 0.000001))

})


#' test-Selectivities
#' @description test new double normal selectivity
#'
test_that("test-Selectivities-DoubleNormalthreeParam", {
  ## Read in mock data
  load(system.file("testdata", "MockSablefishModel.RData",package="SpatialSablefishAssessment"))
  data$model = "TagIntegratedValidate"
  ## no Z or movement
  data$apply_Z_on_tagged_fish = 0
  ## two time-blocks
  data$trwl_sel_type = c(5)
  data$apply_fixed_movement = 1 ## no movement initial age-structure should be only a function of ageing and M
  parameters$ln_trwl_sel_pars = array(0, dim = c(1, 3, 2))
  parameters$ln_trwl_sel_pars[1,1,1] = log(10)
  parameters$ln_trwl_sel_pars[1,2,1] = log(1.4)
  parameters$ln_trwl_sel_pars[1,3,1] = log(4.3)
  parameters$ln_trwl_sel_pars[1,1,2] = log(12)
  parameters$ln_trwl_sel_pars[1,2,2] = log(2)
  parameters$ln_trwl_sel_pars[1,3,2] = log(6.3)



  test_model <- TMB::MakeADFun(data = data,
                               parameters = parameters,
                               DLL = "SpatialSablefishAssessment_TMBExports", silent  = T)

  test_report = test_model$report()

  ## double normal
  # Time-block 1
  # Female
  expect_equal(test_report$sel_trwl_f[,1], d_norm_alt_sel(bins = data$ages, mu = exp(parameters$ln_trwl_sel_pars[1,1,2]), sig_r = exp(parameters$ln_trwl_sel_pars[1,2,2]), sig_l = exp(parameters$ln_trwl_sel_pars[1,3,2])))
  # Male
  expect_equal(test_report$sel_trwl_m[,1], d_norm_alt_sel(bins = data$ages, mu = exp(parameters$ln_trwl_sel_pars[1,1,1]), sig_r = exp(parameters$ln_trwl_sel_pars[1,2,1]), sig_l = exp(parameters$ln_trwl_sel_pars[1,3,1])))
})
