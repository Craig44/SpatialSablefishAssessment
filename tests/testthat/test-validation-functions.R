

#' check_dims
#' @description this will validate some of the checking helper functions
#'
test_that("check_dims", {
  n_regions = 10
  ## 2 dimensional case
  move_matrix = matrix(0, nrow = n_regions, ncol = n_regions)
  ## check that its true when correct
  expect_true(check_dim(move_matrix, c(n_regions, n_regions))$result)
  ## check that if fails when dims are wrong
  expect_false(check_dim(move_matrix, c(n_regions, n_regions - 1))$result)

  ## 3 dimensional case
  n_ages = 20
  n_regions = 5
  n_years = 30

  n_at_age = array(20, dim = c(n_ages, n_regions, n_years))
  expect_true(check_dim(n_at_age, c(n_ages, n_regions, n_years))$result)
  ## if wrong
  expect_false(check_dim(n_at_age, c(n_regions, n_ages, n_years))$result)

  ## check vector function
  vec_vals = rep(0, n_years)
  expect_true(check_length(vec_vals, n_years)$result)
  expect_false(check_length(0, n_years)$result)
})



#' validate-inputs
#' @description this will validate some of the checking helper functions
#'
test_that("validate-inputs", {
  load(system.file("testdata", "MockSablefishModel.RData",package="SpatialSablefishAssessment"))
  ## check the models
  expect_true(validate_input_data_and_parameters(data, parameters))
})
