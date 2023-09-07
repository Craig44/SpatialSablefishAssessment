

#' test-TagIntegratedInitialisation-no_movement
#' @description this will check the initialisation works as expected
#'
test_that("test-TagIntegratedInitialisation-no_movement", {
  ## Read in mock data
  load(system.file("testdata", "MockSablefishModel.RData",package="SpatialSablefishAssessment"))
  data$model = "TagIntegratedValidate"
  ## no Z or movement
  data$apply_Z_on_tagged_fish = 0
  data$apply_fixed_movement = 1 ## no movement initial age-structure should be only a function of ageing and M


  test_model <- TMB::MakeADFun(data = data,
                               parameters = parameters,
                               DLL = "SpatialSablefishAssessment_TMBExports", silent  = T)

  test_report = test_model$report()
  M = 0.104884
  expect_init_age = NULL
  for(r in 1:data$n_regions) {
    this_region_init = c(test_report$mean_rec[r], test_report$mean_rec[r] * exp(-M * 1:(length(data$ages) - 1)))
    ## plus group
    this_region_init[length(data$ages)] = test_report$mean_rec[r] * exp(-(length(data$ages) - 1)*M) / (1 - exp(-M))
    ## save it
    expect_init_age = cbind(expect_init_age, this_region_init / 2) # divide by 2 split between male and female
  }
  ## run test
  expect_equal(expect_init_age, test_report$init_natage_f, tolerance = 0.001)
})


#' test-TagIntegratedInitialisation-movement
#' @description this will check the initialisation works as expected
#'
test_that("test-TagIntegratedInitialisation-movement", {
  ## Read in mock data
  load(system.file("testdata", "MockSablefishModel.RData",package="SpatialSablefishAssessment"))
  data$model = "TagIntegratedValidate"
  ## no Z or movement
  data$apply_Z_on_tagged_fish = 0
  data$apply_fixed_movement = 0 ## no movement initial age-structure should be only a function of ageing and M


  test_model <- TMB::MakeADFun(data = data,
                               parameters = parameters,
                               DLL = "SpatialSablefishAssessment_TMBExports", silent  = T)

  test_report = test_model$report()

  M = 0.104884
  plus_group_age = length(data$ages)
  N_age = matrix(0, nrow = data$n_regions, ncol = plus_group_age)
  update_N_age = N_age

  for(i in 1:(plus_group_age)) {
    ## interpolate SSB
    current_ssb = rowSums(sweep(N_age/2 * exp(-M)^data$spawning_time_proportion[1], MARGIN = 2, test_report$weight_maturity_prod_f[,1], FUN = "*"))

    # recruitment
    update_N_age[,1] = test_report$mean_rec
    # ageing and mortality
    update_N_age[,2:plus_group_age] = N_age[,1:(plus_group_age - 1)] * exp(-M)
    # plus group
    update_N_age[,plus_group_age] = update_N_age[,plus_group_age] + N_age[,plus_group_age] * exp(-M)
    # movement
    N_age = t(test_report$movement_matrix[,,1]) %*% update_N_age
  }
  ## calculate one more annual cycle
  # recruitment
  update_N_age[,1] = test_report$mean_rec
  # ageing and mortality
  update_N_age[,2:plus_group_age] = N_age[,1:(plus_group_age - 1)] * exp(-M)
  # plus group
  update_N_age[,plus_group_age] = update_N_age[,plus_group_age] + N_age[,plus_group_age] * exp(-M)

  # movement
  update_N_age = t(test_report$movement_matrix[,,1]) %*% update_N_age
  ## approximate!
  c = update_N_age[,plus_group_age] / N_age[,plus_group_age] - 1
  update_N_age[,plus_group_age] = N_age[,plus_group_age] * 1 / (1 - c)

  ## run test on partitiion
  expect_equal(t(update_N_age/2), test_report$init_natage_f, tolerance = 0.001)

  ## test Bzero calculation
  Bzero = rowSums(sweep(update_N_age/2 * exp(-M)^data$spawning_time_proportion[1], MARGIN = 2, test_report$weight_maturity_prod_f[,1], FUN = "*"))

  expect_equal(test_report$Bzero, Bzero, tolerance = 0.001)

})

#' test-TagIntegratedInitialisation-no_movement-init_rec_devs
#' @description this will check the initialisation work with init-devs
#'
test_that("test-TagIntegratedInitialisation-no_movement-init_rec_devs", {
  ## Read in mock data
  load(system.file("testdata", "MockSablefishModel.RData",package="SpatialSablefishAssessment"))
  data$model = "TagIntegratedValidate"
  ## no Z or movement
  data$apply_Z_on_tagged_fish = 0
  data$apply_fixed_movement = 1 ## no movement initial age-structure should be only a function of ageing and M
  ## start with one recruitment dev this should be applied to all age elements 2:(plus_group - 1)
  data$n_init_rec_devs = 1
  parameters$ln_init_rec_dev = log(1.2)

  test_model <- TMB::MakeADFun(data = data,
                               parameters = parameters,
                               DLL = "SpatialSablefishAssessment_TMBExports", silent  = T)

  test_report = test_model$report()
  M = 0.104884
  equlibrium_init_age = NULL
  for(r in 1:data$n_regions) {
    this_region_init = c(test_report$mean_rec[r], test_report$mean_rec[r] * exp(-M * 1:(length(data$ages) - 1)))
    ## plus group
    this_region_init[length(data$ages)] = test_report$mean_rec[r] * exp(-(length(data$ages) - 1)*M) / (1 - exp(-M))
    ## save it
    equlibrium_init_age = cbind(equlibrium_init_age, this_region_init / 2) # divide by 2 split between male and female
  }
  ## account for multiplier
  expect_init_age = equlibrium_init_age
  expect_init_age[2:(length(data$ages) - 1),] = expect_init_age[2:(length(data$ages) - 1),] * 1.2
  ## run test
  expect_equal(expect_init_age, test_report$init_natage_f, tolerance = 0.001)

  ## 10 devs
  devs = rnorm(n = 10, mean = 0, sd = 1)
  data$n_init_rec_devs = 10
  parameters$ln_init_rec_dev = (devs)

  test_model <- TMB::MakeADFun(data = data,
                               parameters = parameters,
                               DLL = "SpatialSablefishAssessment_TMBExports", silent  = T)

  test_report = test_model$report()
  ## account for multiplier
  expect_init_age = equlibrium_init_age
  expect_init_age[2:10,] = expect_init_age[2:10,] * exp( devs[1:9])
  ## the last one applies to a range of older ages
  expect_init_age[11:(length(data$ages) - 1),] = expect_init_age[11:(length(data$ages) - 1),] * exp(devs[10])

  ## run test
  expect_equal(expect_init_age, test_report$init_natage_f, tolerance = 0.001)

})
