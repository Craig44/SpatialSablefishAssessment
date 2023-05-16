#'
#' test stock recruit functions
#'


#' average-with-devs
#' @description test the status quo used in sablefish assessment
#'
test_that("stock-recruit-average-with-devs", {
  ## Read in mock data
  load(system.file("testdata", "MockSablefishModel.RData",package="SpatialSablefishAssessment"))
  data$SrType = 3

  test_model <- TMB::MakeADFun(data = data,
                         parameters = parameters,
                         DLL = "SpatialSablefishAssessment_TMBExports", silent  = T)

  test_report = test_model$report()
  expect_results = matrix(test_report$mean_rec * test_report$recruitment_multipliers, ncol = ncol(test_report$recruitment_yr), byrow = T)
  expect_equal(test_report$recruitment_yr, expect_results, tolerance = 0.001)
})

#' BH-with-steepness
#' @description test the Beverton Holt stock recruit relationshisp
#'
test_that("stock-recruit-BH-with-steepness", {
  ## Read in mock data
  load(system.file("testdata", "MockSablefishModel.RData",package="SpatialSablefishAssessment"))
  data$SrType = 2
  parameters$trans_SR_pars = rep(qlogis(0.8)) ## steepness
  test_model <- TMB::MakeADFun(data = data,
                               parameters = parameters,
                               DLL = "SpatialSablefishAssessment_TMBExports", silent  = T)

  test_report = test_model$report()
  expect_results = matrix(0, ncol = ncol(test_report$recruitment_yr), nrow = nrow(test_report$recruitment_yr), byrow = T)

  for(y in 1:length(data$years)) {
    for(i in 1:data$n_regions) {
      SR = 0
      if(y <= min(data$ages)) {
        SR = BH(test_report$Binit[i], test_report$Bzero[i], test_report$SR_pars)
      } else {
        SR = BH(test_report$SSB_yr[y - min(data$ages), i], test_report$Bzero[i], test_report$SR_pars)
      }
      expect_results[y,i] = test_report$mean_rec[i] * SR * test_report$recruitment_multipliers[i,y]
    }
  }

  expect_equal(test_report$recruitment_yr, expect_results, tolerance = 0.001)
})
