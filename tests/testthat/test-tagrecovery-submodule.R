#'
#' This unit-test will run some basic unit-tests on the tag-recovery submodule
#'

#' single-release-test-recovery-reporting-rate
#' @description this will check the predicted tag-recovery observations are working as expected
#'
test_that("single-release-test-recovery-reporting-rate", {
  ## Read in mock data
  load(system.file("testdata", "MockSablefishModel.RData",package="SpatialSablefishAssessment"))
  data$model = "TagIntegratedValidate"
  ## no Z or movement
  data$apply_Z_on_tagged_fish = 0
  data$apply_fixed_movement = 1
  ## this assumes no movement
  fixed_movement_matrix = matrix(0, nrow = data$n_regions, ncol = data$n_regions);
  diag(fixed_movement_matrix) = 1
  data$fixed_movement_matrix = array(0, dim = c(data$n_regions,data$n_regions,1))
  data$fixed_movement_matrix[,,1] = fixed_movement_matrix
  data$apply_fishery_tag_reporting = 0 ## all tagged fish will be recovered not a function of F
  ## turn off tag shedding and initial mortality
  data$initial_tag_induced_mortality = rep(0.0, sum(data$tag_release_event_this_year))
  data$annual_tag_shedding_rate = 0.0
  ##
  parameters$logistic_tag_reporting_rate = matrix(logit(0.2), nrow = data$n_regions, ncol = sum(data$tag_recovery_indicator))

  ## because there is no movement we need to turn off tag-recoveries in non release years
  tag_recovery_years = 2011:2020
  tag_release_years = data$years[which(data$tag_release_event_this_year == 1)]
  # drop any recovery years before release years
  #tag_recovery_years = tag_recovery_years[which(tag_recovery_years %in% (tag_release_years + 1))] ## the plus one is because we don't allow a recovery unless after a year at release
  #
  data$tag_recovery_indicator_by_year = rep(0, length(data$years)) ## no tag releases
  data$tag_recovery_indicator_by_year[data$years %in% tag_recovery_years] = 1
  data$tag_recovery_indicator = array(0, dim = c(data$n_regions * (data$n_years_to_retain_tagged_cohorts_for + 1), data$n_regions, length(tag_recovery_years)))
  data$obs_tag_recovery = array(10, dim = c(data$n_regions * (data$n_years_to_retain_tagged_cohorts_for + 1), data$n_regions, length(tag_recovery_years)))
  data$apply_tag_reporting_rate = 1;

  ## track each tag cohort for 6 years opst release
  for(y_ndx in 1:length(tag_release_years)) {
    for(r_ndx in 1:data$n_regions) {
      for(recovery_yr_ndx in 1:length(tag_recovery_years)) {
        diff_ = tag_recovery_years[recovery_yr_ndx] - tag_release_years
        diff_ = min(c(diff_ + 1, data$n_years_to_retain_tagged_cohorts_for + 1))
        release_event_ndx = get_tag_release_ndx(r_ndx, diff_, data$n_regions)
        data$tag_recovery_indicator[release_event_ndx, r_ndx, recovery_yr_ndx] = 1
      }
    }
  }
  #sum(data$tag_recovery_indicator)

  test_model <- TMB::MakeADFun(data = data,
                               parameters = parameters,
                               DLL = "SpatialSablefishAssessment_TMBExports", silent  = T)

  test_report = test_model$report()
  release_event_label = "2010-BS"

  ## since no movement check all recoveries are in release regions
  ## no mortality or movement. So tag-recoveries n should equal releases plus ageing with reporting rate
  plt = plot_tag_recovery_obs(MLE_report = test_report, region_key = region_key, sex = "both", release_ndx_to_plot = 1:10)
  recovery_df = plt$data %>% dplyr::filter(recovery_region == release_region, release_event == release_event_label)
  ## with no movement the release and recovery periods should be the same
  final_year = max(data$years)
  ## get the tag-release_event ndx in the final year and check the tagged partition is as expected
  tag_release_years_in_final_year = final_year - data$years[which(data$tag_release_event_this_year == 1)] + 1
  ## we actually age fish at the end of the final year so incerment these
  tag_release_years_in_final_year = tag_release_years_in_final_year + 1
  ## compare numbers at age for one of the releases that were released in region 1
  init_n_age = data$male_tagged_cohorts_by_age[,1,1]
  ## age this partition
  release_year = data$years[which(data$tag_release_event_this_year == 1)][1]
  ##
  times_to_age = final_year - release_year + 1 ## plus one for the last year
  n_ages = length(data$ages)
  for(i in 1:times_to_age) {
    temp_tag_partition = init_n_age
    init_n_age = vector(length = length(init_n_age)) ## clear it
    init_n_age[2:length(init_n_age)] = temp_tag_partition[1:(length(init_n_age) - 1)]
    init_n_age[length(init_n_age)] = temp_tag_partition[length(init_n_age) - 1] +  temp_tag_partition[length(init_n_age)]
    expected_result = init_n_age * 0.2 ## annual shedding rate
    this_recovery_year = release_year + i
    test_data = recovery_df %>% dplyr::filter(recovery_year == this_recovery_year) %>% dplyr::select(predicted)
    if(nrow(test_data) > 0) {
      expect_equal(test_data$predicted, sum(expected_result) * 2, tolerance = 0.001)
    }
  }
  tmp = plt$data  %>% group_by(release_event, recovery_region, recovery_year) %>% summarise(obs = sum(observed), pred = sum(predicted))
  ## check that the Poisson likelihood evaluation is as expected
  expect_equal(test_report$nll[8], -1 * sum(dpois(tmp$obs, tmp$pred, log = T)), tolerance = 0.0001)


  ## with negtive binomial likelihood
  data$tag_likelihood = 1
  dispersion_param = 2
  parameters$ln_tag_phi = log(dispersion_param)

  test_model <- TMB::MakeADFun(data = data,
                               parameters = parameters,
                               DLL = "SpatialSablefishAssessment_TMBExports", silent  = T)

  test_report = test_model$report()
  expect_equal(test_report$nll[8], -1 * sum(dnbinom(tmp$obs, size = dispersion_param, mu = tmp$pred, log = T)), tolerance = 0.0001)

  dispersion_param = 0.001
  parameters$ln_tag_phi = log(dispersion_param)

  test_model <- TMB::MakeADFun(data = data,
                               parameters = parameters,
                               DLL = "SpatialSablefishAssessment_TMBExports", silent  = T)

  test_report = test_model$report()
  expect_equal(test_report$nll[8], -1 * sum(dnbinom(tmp$obs, size = dispersion_param, mu = tmp$pred, log = T)), tolerance = 0.0001)

})


#' single-release-F-reporting
#' @description this will check the predicted tag-recovery observations are working when returns are from the fishery
#'
test_that("single-release-F-reporting", {
  ## Read in mock data
  load(system.file("testdata", "MockSablefishModel.RData",package="SpatialSablefishAssessment"))
  data$model = "TagIntegratedValidate"
  ## no Z or movement
  data$apply_Z_on_tagged_fish = 0
  data$apply_fixed_movement = 1
  ## this assumes no movement
  fixed_movement_matrix = matrix(0, nrow = data$n_regions, ncol = data$n_regions);
  diag(fixed_movement_matrix) = 1
  data$fixed_movement_matrix = array(0, dim = c(data$n_regions,data$n_regions,1))
  data$fixed_movement_matrix[,,1] = fixed_movement_matrix
  data$apply_fishery_tag_reporting = 1 ## all tagged fish will be recovered by fixed gear fishery
  ## turn off tag shedding and initial mortality
  data$initial_tag_induced_mortality = rep(0.0, sum(data$tag_release_event_this_year))
  data$annual_tag_shedding_rate = 0.0
  data$apply_tag_reporting_rate = 1;
  ##
  report_rate = 0.2
  parameters$logistic_tag_reporting_rate = matrix(logit(report_rate), nrow = data$n_regions, ncol = sum(data$tag_recovery_indicator))


  test_model <- TMB::MakeADFun(data = data,
                               parameters = parameters,
                               DLL = "SpatialSablefishAssessment_TMBExports", silent  = T)

  test_report = test_model$report()
  release_event_label = "2010-BS"

  ## since no movement check all recoveries are in release regions
  ## no mortality or movement. So tag-recoveries n should equal releases plus ageing with reporting rate
  plt = plot_tag_recovery_obs(MLE_report = test_report, region_key = region_key, sex = "both", release_ndx_to_plot = 1:5)
  recovery_df = plt$data %>% dplyr::filter(recovery_region == release_region, release_event == release_event_label)
  ## with no movement the release and recovery periods should be the same
  final_year = max(data$years)
  ## get the tag-release_event ndx in the final year and check the tagged partition is as expected
  tag_release_years_in_final_year = final_year - data$years[which(data$tag_release_event_this_year == 1)] + 1
  ## we actually age fish at the end of the final year so incerment these
  tag_release_years_in_final_year = tag_release_years_in_final_year + 1
  ## compare numbers at age for one of the releases that were released in region 1
  init_n_age = data$male_tagged_cohorts_by_age[,1,1]
  ## age this partition
  release_year = data$years[which(data$tag_release_event_this_year == 1)][1]
  ##
  times_to_age = final_year - release_year + 1 ## plus one for the last year
  n_ages = length(data$ages)
  for(i in 1:times_to_age) {
    temp_tag_partition = init_n_age
    init_n_age = vector(length = length(init_n_age)) ## clear it
    init_n_age[2:length(init_n_age)] = temp_tag_partition[1:(length(init_n_age) - 1)]
    init_n_age[length(init_n_age)] = temp_tag_partition[length(init_n_age) - 1] +  temp_tag_partition[length(init_n_age)]
    expected_result = (init_n_age * test_report$F_fixed_m[,1,i] / test_report$Z_m[,1,i] * (1 - test_report$S_m[,1,i]) + init_n_age * test_report$F_fixed_f[,1,i] / test_report$Z_f[,1,i] * (1 - test_report$S_f[,1,i])) * report_rate## annual shedding rate
    this_recovery_year = release_year + i
    test_data = recovery_df %>% dplyr::filter(recovery_year == this_recovery_year) %>% dplyr::select(predicted)
    if(nrow(test_data) > 0) {
      expect_equal(test_data$predicted, sum(expected_result), tolerance = 0.001)
    }
  }
})


#' single-release-test-multinomial-likelihood-no-movement
#' @description this will check the predicted and likelihood calculation for when tag_likelihood == 3
#'
test_that("single-release-test-multinomial-likelihood-no-movement", {
  ## Read in mock data
  load(system.file("testdata", "MockSablefishModel.RData",package="SpatialSablefishAssessment"))
  data$model = "TagIntegratedValidate"
  ## no Z or movement
  data$apply_Z_on_tagged_fish = 1
  data$apply_fixed_movement = 1
  ## this assumes no movement
  fixed_movement_matrix = matrix(0, nrow = data$n_regions, ncol = data$n_regions);
  diag(fixed_movement_matrix) = 1
  data$fixed_movement_matrix = array(0, dim = c(data$n_regions,data$n_regions,1))
  data$fixed_movement_matrix[,,1] = fixed_movement_matrix
  data$apply_fishery_tag_reporting = 1 ## all tagged fish will be recovered not a function of F
  ## turn off tag shedding and initial mortality
  data$initial_tag_induced_mortality = rep(0.0, sum(data$tag_release_event_this_year))
  data$annual_tag_shedding_rate = 0.0
  ##
  report_rate = 0.9
  parameters$logistic_tag_reporting_rate = matrix(logit(report_rate), nrow = data$n_regions, ncol = sum(data$tag_recovery_indicator))

  ## because there is no movement we need to turn off tag-recoveries in non release years
  tag_recovery_years = 2011:2020
  tag_release_years = data$years[which(data$tag_release_event_this_year == 1)]
  # drop any recovery years before release years
  #tag_recovery_years = tag_recovery_years[which(tag_recovery_years %in% (tag_release_years + 1))] ## the plus one is because we don't allow a recovery unless after a year at release
  data$tag_recovery_indicator_by_year = rep(0, length(data$years)) ## no tag releases
  data$tag_recovery_indicator_by_year[data$years %in% tag_recovery_years] = 1
  data$tag_recovery_indicator = array(0, dim = c(length(data$years), data$n_regions))
  data$tag_recovery_indicator[1,] = 1
  data$obs_tag_recovery = array(10, dim = c(data$n_regions * data$n_years_to_retain_tagged_cohorts_for + 1, data$n_regions, length(data$years)))

  for(r in 1:data$n_regions)
    data$obs_tag_recovery[,r,1] = (data$obs_tag_recovery[,r,1] / sum(data$obs_tag_recovery[,r,1])) * (sum(data$male_tagged_cohorts_by_age[,r,1]) + sum(data$female_tagged_cohorts_by_age[,r,1]))

  data$apply_tag_reporting_rate = 0;
  data$tag_likelihood = 2
  #sum(data$tag_recovery_indicator)

  test_model <- TMB::MakeADFun(data = data,
                               parameters = parameters,
                               DLL = "SpatialSablefishAssessment_TMBExports", silent  = T)

  test_report = test_model$report()
  release_event_label = "2010-BS"

  ## check predictions sum = 1
  for(r in 1:data$n_regions)
    expect_equal(sum(test_report$pred_tag_recovery[,r,1]), 1, 0.001)

  ## since no movement check all recoveries are in release regions
  tag_df = get_tag_recovery_obs_fitted_values(test_report, region_key)

  init_m_age = data$male_tagged_cohorts_by_age[,1,1]
  init_f_age = data$female_tagged_cohorts_by_age[,1,1]
  released_fish = sum(data$male_tagged_cohorts_by_age[,1,1]) + sum(data$female_tagged_cohorts_by_age[,1,1])
  times_to_age = data$n_years_to_retain_tagged_cohorts_for
  n_ages = length(data$ages)
  expected_recaptures = c()
  for(i in 1:times_to_age) {
    temp_m_tag_partition = init_m_age
    temp_f_tag_partition = init_f_age

    expected_result = (init_m_age * test_report$F_fixed_m[,1,i] / test_report$Z_m[,1,i] * (1 - test_report$S_m[,1,i]) + init_f_age * test_report$F_fixed_f[,1,i] / test_report$Z_f[,1,i] * (1 - test_report$S_f[,1,i]))## annual shedding rate
    result = posfun(sum(expected_result), 0.00001)

    expected_recaptures = c(expected_recaptures, c(result, rep(posfun(0, 0.00001), 4)));

    ## age and Z
    init_m_age = vector(length = length(init_m_age)) ## clear it
    init_m_age[2:length(init_m_age)] = temp_m_tag_partition[1:(length(init_m_age) - 1)] * test_report$S_m[1:(length(init_m_age) - 1),1,i]
    init_m_age[length(init_m_age)] = temp_m_tag_partition[length(init_m_age) - 1] * test_report$S_m[length(init_m_age) - 1, 1, i] +  temp_m_tag_partition[length(init_m_age)] * test_report$S_m[length(init_m_age),1,i]

    init_f_age = vector(length = length(init_f_age)) ## clear it
    init_f_age[2:length(init_f_age)] = temp_f_tag_partition[1:(length(init_f_age) - 1)] * test_report$S_f[1:(length(init_m_age) - 1),1,i]
    init_f_age[length(init_f_age)] = temp_f_tag_partition[length(init_f_age) - 1] * test_report$S_f[length(init_m_age) - 1, 1, i] +  temp_f_tag_partition[length(init_f_age)] * test_report$S_m[length(init_m_age),1,i]
  }
  ## divide by released fish
  expected_recaptures = expected_recaptures / released_fish
  ## not-recaptured group
  expected_recaptures = c(expected_recaptures, 1 -  sum(expected_recaptures))
  tmp = tag_df %>% dplyr::filter(release_year == 2010, release_region == "BS")

  ## check that the expected values are the same
  expect_equal(expected_recaptures, tmp$predicted, tolerance = 0.0001)
})



#' single-release-test-multinomial-likelihood-movement
#' @description this will check the predicted and likelihood calculation for when tag_likelihood == 3
#'
test_that("single-release-test-multinomial-likelihood-movement", {
  ## Read in mock data
  load(system.file("testdata", "MockSablefishModel.RData",package="SpatialSablefishAssessment"))
  data$model = "TagIntegratedValidate"
  ## no Z or movement
  data$apply_Z_on_tagged_fish = 1
  data$apply_fixed_movement = 0
  data$apply_fishery_tag_reporting = 1 ## all tagged fish will be recovered not a function of F
  ## turn off tag shedding and initial mortality
  data$initial_tag_induced_mortality = rep(0.0, sum(data$tag_release_event_this_year))
  data$annual_tag_shedding_rate = 0.0
  ##
  report_rate = 0.9
  parameters$logistic_tag_reporting_rate = matrix(logit(report_rate), nrow = data$n_regions, ncol = sum(data$tag_recovery_indicator))

  ## because there is no movement we need to turn off tag-recoveries in non release years
  tag_recovery_years = 2011:2020
  tag_release_years = data$years[which(data$tag_release_event_this_year == 1)]
  # drop any recovery years before release years
  #tag_recovery_years = tag_recovery_years[which(tag_recovery_years %in% (tag_release_years + 1))] ## the plus one is because we don't allow a recovery unless after a year at release
  data$tag_recovery_indicator_by_year = rep(0, length(data$years)) ## no tag releases
  data$tag_recovery_indicator_by_year[data$years %in% tag_recovery_years] = 1
  data$tag_recovery_indicator = array(0, dim = c(length(data$years), data$n_regions))
  data$tag_recovery_indicator[1,] = 1
  data$obs_tag_recovery = array(10, dim = c(data$n_regions * data$n_years_to_retain_tagged_cohorts_for + 1, data$n_regions, length(data$years)))

  for(r in 1:data$n_regions)
    data$obs_tag_recovery[,r,1] = (data$obs_tag_recovery[,r,1] / sum(data$obs_tag_recovery[,r,1])) * (sum(data$male_tagged_cohorts_by_age[,r,1]) + sum(data$female_tagged_cohorts_by_age[,r,1]))

  data$apply_tag_reporting_rate = 0;
  data$tag_likelihood = 2
  #sum(data$tag_recovery_indicator)

  test_model <- TMB::MakeADFun(data = data,
                               parameters = parameters,
                               DLL = "SpatialSablefishAssessment_TMBExports", silent  = T)

  test_report = test_model$report()
  release_event_label = "2010-BS"

  ## check predictions sum = 1
  for(r in 1:data$n_regions)
    expect_equal(sum(test_report$pred_tag_recovery[,r,1]), 1, 0.001)

  ## since no movement check all recoveries are in release regions
  tag_df = get_tag_recovery_obs_fitted_values(test_report, region_key)
  nll = 0;
  for(j in 1:data$n_regions) {
    init_m_age = init_f_age = matrix(0, nrow = data$n_regions, ncol = length(data$ages))
    init_m_age[j,] = data$male_tagged_cohorts_by_age[,j,1]
    init_f_age[j,] = data$female_tagged_cohorts_by_age[,j,1]

    released_fish = sum(data$male_tagged_cohorts_by_age[,j,1]) + sum(data$female_tagged_cohorts_by_age[,j,1])
    times_to_age = data$n_years_to_retain_tagged_cohorts_for
    n_ages = length(data$ages)
    expected_recaptures = c()
    for(i in 1:times_to_age) {
      temp_m_tag_partition = init_m_age
      temp_f_tag_partition = init_f_age
      for(r in 1:data$n_regions) {
        expected_result = (init_m_age[r,] * test_report$F_fixed_m[,r,i] / test_report$Z_m[,r,i] * (1 - test_report$S_m[,r,i]) + init_f_age[r,] * test_report$F_fixed_f[,r,i] / test_report$Z_f[,r,i] * (1 - test_report$S_f[,r,i]))
        result = sum(expected_result)
        result = posfun(result, 0.00001)
        expected_recaptures = c(expected_recaptures, result);
        ## age and Z
        init_m_age[r,] = vector(length = ncol(init_m_age)) ## clear it
        init_m_age[r,2:ncol(init_m_age)] = temp_m_tag_partition[r,1:(ncol(init_m_age) - 1)] * test_report$S_m[1:(ncol(init_m_age) - 1),r,i]
        init_m_age[r,ncol(init_m_age)] = temp_m_tag_partition[r,ncol(init_m_age) - 1] * test_report$S_m[ncol(init_m_age) - 1, r, i] +  temp_m_tag_partition[r,ncol(init_m_age)] * test_report$S_m[ncol(init_m_age),r,i]

        init_f_age[r,] = vector(length = ncol(init_f_age)) ## clear it
        init_f_age[r,2:ncol(init_f_age)] = temp_f_tag_partition[r,1:(ncol(init_f_age) - 1)] * test_report$S_f[1:(ncol(init_m_age) - 1),1,i]
        init_f_age[r,ncol(init_f_age)] = temp_f_tag_partition[r,ncol(init_f_age) - 1] * test_report$S_f[ncol(init_m_age) - 1, 1, i] +  temp_f_tag_partition[r,ncol(init_f_age)] * test_report$S_m[ncol(init_m_age),1,i]
      }
      ## movement
      init_m_age = test_report$movement_matrix[,,1] %*%  init_m_age
      init_f_age = test_report$movement_matrix[,,1] %*%  init_f_age
    }
    ## divide by released fish
    expected_recaptures = expected_recaptures / released_fish
    ## not-recaptured group
    expected_recaptures = c(expected_recaptures, 1 -  sum(expected_recaptures))
    tmp = tag_df %>% dplyr::filter(release_year == 2010, release_region == region_key$area[j])
    expect_equal(expected_recaptures, tmp$predicted, tolerance = 0.01)
  }
  nll = 0
  for(j in 1:data$n_regions)
    nll = nll - dmultinom_upd(test_report$obs_tag_recovery[,j,1], prob = test_report$pred_tag_recovery[,j,1], log = T)

  ## check that the expected values are the same
  expect_equal(nll, test_report$nll[8], 0.0001);
})

