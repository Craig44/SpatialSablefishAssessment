

#' test-TagIntegratedAnnualCycle-fixed-movement
#' @description this is just to tell me when something changes the annual cycle.
#'
test_that("test-TagIntegratedAnnualCycle-fixed-movement", {
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
  ## freeze the last years partition. I haven't actually independently checked this.
  ## this will only tell you when things have changed, not whether anything is correct

  ## paste(round(test_report$natage_f[,,12],3), collapse = ", ")
  expected = c(0, 3.362, 2.88, 2.434, 2.035, 1.687, 1.388, 1.138, 0.931, 0.761, 0.623, 0.51, 0.893, 0.773, 0.68, 0.606, 0.547, 0.498, 0.456, 0.42, 0.388, 0.359, 0.331, 0.306, 0.283, 0.261, 0.24, 0.22, 0.202, 1.967, 0, 3.303, 2.829, 2.392, 2, 1.657, 1.364, 1.118, 0.914, 0.748, 0.612, 0.501, 0.877, 0.759, 0.668, 0.595, 0.537, 0.489, 0.448, 0.413, 0.381, 0.352, 0.326, 0.301, 0.278, 0.256, 0.236, 0.216, 0.199, 1.932, 0, 2.345, 2.008, 1.698, 1.42, 1.176, 0.968, 0.794, 0.649, 0.531, 0.434, 0.356, 0.623, 0.539, 0.474, 0.423, 0.381, 0.347, 0.318, 0.293, 0.271, 0.25, 0.231, 0.214, 0.197, 0.182, 0.167, 0.154, 0.141, 1.372, 0, 4.056, 3.473, 2.937, 2.455, 2.034, 1.675, 1.373, 1.123, 0.918, 0.751, 0.616, 1.077, 0.932, 0.82, 0.731, 0.66, 0.601, 0.551, 0.507, 0.468, 0.432, 0.4, 0.369, 0.341, 0.314, 0.289, 0.266, 0.244, 2.372, 0, 2.976, 2.548, 2.155, 1.801, 1.493, 1.229, 1.007, 0.824, 0.674, 0.551, 0.452, 0.79, 0.684, 0.602, 0.536, 0.484, 0.441, 0.404, 0.372, 0.343, 0.317, 0.293, 0.271, 0.25, 0.231, 0.212, 0.195, 0.179, 1.74)
  expected_mat = matrix(expected, ncol = data$n_regions, byrow = F)

  ## run test
  expect_equal(expected_mat, test_report$natage_f[,,12], tolerance = 0.001)
})


#' test-TagIntegratedAnnualCycle-movement
#' @description this will check the initialisation works as expected
#'
test_that("test-TagIntegratedAnnualCycle-movement", {
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
  ## paste(round(test_report$natage_f[,,12],3), collapse = ", ")
  expected = c(0, 3.354, 2.865, 2.416, 2.016, 1.666, 1.369, 1.12, 0.914, 0.746, 0.609, 0.498, 0.869, 0.751, 0.66, 0.587, 0.529, 0.481, 0.441, 0.405, 0.373, 0.345, 0.318, 0.294, 0.271, 0.25, 0.229, 0.211, 0.193, 1.867, 0, 3.289, 2.804, 2.361, 1.966, 1.623, 1.331, 1.087, 0.886, 0.722, 0.589, 0.482, 0.837, 0.723, 0.634, 0.564, 0.508, 0.461, 0.422, 0.388, 0.357, 0.33, 0.304, 0.28, 0.258, 0.238, 0.219, 0.201, 0.184, 1.768, 0, 2.409, 2.114, 1.829, 1.562, 1.32, 1.107, 0.923, 0.767, 0.637, 0.529, 0.439, 0.787, 0.689, 0.612, 0.551, 0.502, 0.461, 0.426, 0.395, 0.367, 0.342, 0.318, 0.296, 0.274, 0.254, 0.235, 0.217, 0.2, 2.012, 0, 4.006, 3.39, 2.833, 2.343, 1.921, 1.566, 1.271, 1.031, 0.835, 0.678, 0.551, 0.95, 0.817, 0.714, 0.633, 0.567, 0.514, 0.468, 0.429, 0.394, 0.363, 0.334, 0.307, 0.282, 0.259, 0.238, 0.218, 0.199, 1.895, 0, 2.986, 2.565, 2.175, 1.824, 1.516, 1.251, 1.028, 0.843, 0.69, 0.566, 0.465, 0.816, 0.708, 0.624, 0.557, 0.503, 0.459, 0.421, 0.388, 0.359, 0.332, 0.307, 0.284, 0.262, 0.242, 0.223, 0.205, 0.188, 1.8422)
  expected_mat = matrix(expected, ncol = data$n_regions, byrow = F)

  ## run test
  expect_equal(expected_mat, test_report$natage_f[,,12], tolerance = 0.001)

})



#' test-TagIntegratedAnnualCycle-spatial_recruitment-movement
#' @description this is just to tell me when something changes the annual cycle.
#'
test_that("test-TagIntegratedAnnualCycle-spatial_recruitment-movement", {
  ## Read in mock data
  load(system.file("testdata", "MockSablefishModel.RData",package="SpatialSablefishAssessment"))
  data$model = "TagIntegratedValidate"
  ## no Z or movement
  data$apply_Z_on_tagged_fish = 0
  data$apply_fixed_movement = 1 ## no movement initial age-structure should be only a function of ageing and M
  data$global_rec_devs = 0
  parameters$trans_rec_dev = matrix(0, ncol = length(data$years), nrow = data$n_regions)

  test_model <- TMB::MakeADFun(data = data,
                               parameters = parameters,
                               DLL = "SpatialSablefishAssessment_TMBExports", silent  = T)

  test_report = test_model$report()
  ## freeze the last years partition. I haven't actually independently checked this.
  ## this will only tell you when things have changed, not whether anything is correct

  ## paste(round(test_report$natage_f[,,12],3), collapse = ", ")
  expected = c(0, 3.362, 2.88, 2.434, 2.035, 1.687, 1.388, 1.138, 0.931, 0.761, 0.623, 0.51, 0.893, 0.773, 0.68, 0.606, 0.547, 0.498, 0.456, 0.42, 0.388, 0.359, 0.331, 0.306, 0.283, 0.261, 0.24, 0.22, 0.202, 1.967, 0, 3.303, 2.829, 2.392, 2, 1.657, 1.364, 1.118, 0.914, 0.748, 0.612, 0.501, 0.877, 0.759, 0.668, 0.595, 0.537, 0.489, 0.448, 0.413, 0.381, 0.352, 0.326, 0.301, 0.278, 0.256, 0.236, 0.216, 0.199, 1.932, 0, 2.345, 2.008, 1.698, 1.42, 1.176, 0.968, 0.794, 0.649, 0.531, 0.434, 0.356, 0.623, 0.539, 0.474, 0.423, 0.381, 0.347, 0.318, 0.293, 0.271, 0.25, 0.231, 0.214, 0.197, 0.182, 0.167, 0.154, 0.141, 1.372, 0, 4.056, 3.473, 2.937, 2.455, 2.034, 1.675, 1.373, 1.123, 0.918, 0.751, 0.616, 1.077, 0.932, 0.82, 0.731, 0.66, 0.601, 0.551, 0.507, 0.468, 0.432, 0.4, 0.369, 0.341, 0.314, 0.289, 0.266, 0.244, 2.372, 0, 2.976, 2.548, 2.155, 1.801, 1.493, 1.229, 1.007, 0.824, 0.674, 0.551, 0.452, 0.79, 0.684, 0.602, 0.536, 0.484, 0.441, 0.404, 0.372, 0.343, 0.317, 0.293, 0.271, 0.25, 0.231, 0.212, 0.195, 0.179, 1.74)
  expected_mat = matrix(expected, ncol = data$n_regions, byrow = F)

  ## run test
  expect_equal(expected_mat, test_report$natage_f[,,12], tolerance = 0.001)
})


#' test-TagIntegratedAnnualCycle-spatial_recruitment-recruits-movement
#' @description this is just to tell me when something changes the annual cycle.
#'
test_that("test-TagIntegratedAnnualCycle-spatial_recruitment-recruits-movement", {
  ## Read in mock data
  load(system.file("testdata", "MockSablefishModel.RData",package="SpatialSablefishAssessment"))
  data$model = "TagIntegratedValidate"
  ## no Z or movement
  data$apply_Z_on_tagged_fish = 0
  data$apply_fixed_movement = 0
  data$global_rec_devs = 0
  parameters$trans_rec_dev = matrix(0, ncol = length(data$years), nrow = data$n_regions)
  data$do_recruits_move = 1
  test_model <- TMB::MakeADFun(data = data,
                               parameters = parameters,
                               DLL = "SpatialSablefishAssessment_TMBExports", silent  = T)

  test_report = test_model$report()
  ## freeze the last years partition. I haven't actually independently checked this.
  ## this will only tell you when things have changed, not whether anything is correct

  ## paste(round(test_report$natage_f[,,12],3), collapse = ", ")
  expected = c(0, 3.354, 2.865, 2.416, 2.016, 1.666, 1.369, 1.12, 0.914, 0.746, 0.609, 0.498, 0.869, 0.751, 0.66, 0.587, 0.529, 0.481, 0.441, 0.405, 0.373, 0.345, 0.318, 0.294, 0.271, 0.25, 0.229, 0.211, 0.193, 1.867, 0, 3.289, 2.804, 2.361, 1.966, 1.623, 1.331, 1.087, 0.886, 0.722, 0.589, 0.482, 0.837, 0.723, 0.634, 0.564, 0.508, 0.461, 0.422, 0.388, 0.357, 0.33, 0.304, 0.28, 0.258, 0.238, 0.219, 0.201, 0.184, 1.768, 0, 2.409, 2.114, 1.829, 1.562, 1.32, 1.107, 0.923, 0.767, 0.637, 0.529, 0.439, 0.787, 0.689, 0.612, 0.551, 0.502, 0.461, 0.426, 0.395, 0.367, 0.342, 0.318, 0.296, 0.274, 0.254, 0.235, 0.217, 0.2, 2.012, 0, 4.006, 3.39, 2.833, 2.343, 1.921, 1.566, 1.271, 1.031, 0.835, 0.678, 0.551, 0.95, 0.817, 0.714, 0.633, 0.567, 0.514, 0.468, 0.429, 0.394, 0.363, 0.334, 0.307, 0.282, 0.259, 0.238, 0.218, 0.199, 1.895, 0, 2.986, 2.565, 2.175, 1.824, 1.516, 1.251, 1.028, 0.843, 0.69, 0.566, 0.465, 0.816, 0.708, 0.624, 0.557, 0.503, 0.459, 0.421, 0.388, 0.359, 0.332, 0.307, 0.284, 0.262, 0.242, 0.223, 0.205, 0.188, 1.842)
  expected_mat = matrix(expected, ncol = data$n_regions, byrow = F)

  ## run test
  expect_equal(expected_mat, test_report$natage_f[,,12], tolerance = 0.001)
})



#' test-TagIntegratedAnnualCycle-recruitment-sumtozero
#' @description test sum to zero constraint for recruit devs
#'
test_that("test-TagIntegratedAnnualCycle-recruitment-sumtozero", {
  ## Read in mock data
  load(system.file("testdata", "MockSablefishModel.RData",package="SpatialSablefishAssessment"))
  data$model = "TagIntegrated"
  data$rec_devs_sum_to_zero = 1
  ## no Z or movement
  data$apply_Z_on_tagged_fish = 0
  data$apply_fixed_movement = 1 ## no movement initial age-structure should be only a function of ageing and M
  data$global_rec_devs = 0
  parameters$trans_rec_dev = matrix(rnorm(data$n_regions * (length(data$years) - 1)), ncol = length(data$years) - 1, nrow = data$n_regions)

  test_model <- TMB::MakeADFun(data = data,
                               parameters = parameters,
                               DLL = "SpatialSablefishAssessment_TMBExports", silent  = T)

  test_report = test_model$report()

  for(r in 1:data$n_regions) {
    expect_equal(sum(test_report$recruitment_devs[r, ]), 0, tolerance = 0.000001)
    # compare with R function
    expect_equal(test_report$recruitment_devs[r, ], sum_to_zero_QR(parameters$trans_rec_dev[r,]), tolerance = 0.000001)

  }
  ## not global rec devs
  data$global_rec_devs = 1
  parameters$trans_rec_dev = matrix(rnorm( (length(data$years) - 1)), ncol = length(data$years) - 1, nrow = 1)

  test_model <- TMB::MakeADFun(data = data,
                               parameters = parameters,
                               DLL = "SpatialSablefishAssessment_TMBExports", silent  = T)

  test_report = test_model$report()
  for(r in 1:data$n_regions) {
    expect_equal(sum(test_report$recruitment_devs[r, ]), 0, tolerance = 0.000001)
  }
  # compare with R function
  expect_equal(test_report$recruitment_devs[1, ], sum_to_zero_QR(parameters$trans_rec_dev), tolerance = 0.000001)

})
