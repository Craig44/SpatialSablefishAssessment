#'
#' These unit-tests will run some basic unit-tests on the tagging submodule, to make sure
#' The C++ code is working how we think it should be. It is checked against manual R-code for validation
#'


#' single-release-no-movement-and-Z
#' @description this will track taged cohorts in a model with no movement or mortality
#' to check that tagged fish move as we expect. Because they have a slightly abstract
#' partition design this is a useful first validation
#'
test_that("single-release-no-movement-and-Z", {
  ## Read in mock data
  load(system.file("testdata", "MockSablefishModel.RData",package="SpatialSablefishAssessment"))

  data$model = "TagIntegratedValidate"
  ## no Z or movement
  data$apply_Z_on_tagged_fish = 0
  data$apply_fixed_movement = 1
  ## this assumes no movement

  data$apply_fishery_tag_reporting = 0 ## all tagged fish will be recovered not a function of F
  ## turn off tag shedding and initial mortality
  data$initial_tag_induced_mortality = rep(0.0, sum(data$tag_release_event_this_year))
  data$annual_tag_shedding_rate = 0.0



  test_model <- TMB::MakeADFun(data = data,
                         parameters = parameters,
                         DLL = "SpatialSablefishAssessment_TMBExports", silent  = T)

  test_report = test_model$report()
  ## test that the same number of tagged fish are in the tagged partition as was at the beginning
  all_releases = sum(data$male_tagged_cohorts_by_age) + sum(data$female_tagged_cohorts_by_age)
  releases_at_the_end = sum(test_report$tagged_natage_f) + sum(test_report$tagged_natage_m)
  ## should be equal because no Z
  expect_true(all_releases == releases_at_the_end, label = "tagged fish at end of model not as expected")

  ## with no movement the release and recovery periods should be the same
  final_year = max(data$years)
  ## get the tag-release_event ndx in the final year and check the tagged partition is as expected
  tag_release_years_in_final_year = final_year - data$years[which(data$tag_release_event_this_year == 1)] + 1
  ## we actually age fish at the end of the final year so incerment these
  tag_release_years_in_final_year = tag_release_years_in_final_year + 1

  ## if values over n_years_to_retain_tagged_cohorts_for we cache them in a pooled group
  tag_release_years_in_final_year[tag_release_years_in_final_year > (data$n_years_to_retain_tagged_cohorts_for + 1)] = data$n_years_to_retain_tagged_cohorts_for + 1
  get_all_tag_release_events = vector()
  for(i in 1:length(tag_release_years_in_final_year)) {
    for(rel_region_ndx in 1:data$n_regions) {
      ndx = get_tag_release_ndx(rel_region_ndx, tag_release_years_in_final_year[i], data$n_regions)
      get_all_tag_release_events = c(get_all_tag_release_events, ndx)
    }
  }
  ## check again by sex
  expect_true(all(apply(test_report$tagged_natage_f, 3, sum)[get_all_tag_release_events] == apply(data$male_tagged_cohorts_by_age, 2, sum)))
  expect_true(all(apply(test_report$tagged_natage_m, 3, sum)[get_all_tag_release_events] == apply(data$female_tagged_cohorts_by_age, 2, sum)))

  ## compare numbers at age for one of the releases that were released in region 1
  init_n_age = data$male_tagged_cohorts_by_age[,1,1]
  ## age this partition
  release_year = data$years[which(data$tag_release_event_this_year == 1)][1]
  ##
  times_to_age = final_year - release_year + 1 ## plus one for the last year
  for(i in 1:times_to_age) {
    temp_tag_partition = init_n_age
    init_n_age = vector(length = length(init_n_age)) ## clear it
    init_n_age[2:length(init_n_age)] = temp_tag_partition[1:(length(init_n_age) - 1)]
    init_n_age[length(init_n_age)] = temp_tag_partition[length(init_n_age) - 1] +  temp_tag_partition[length(init_n_age)]
  }
  ## check the ageing works as expected
  for(age_ndx in 1:length(init_n_age))
    expect_true(test_report$tagged_natage_m[age_ndx,1,get_all_tag_release_events[1]] == init_n_age[age_ndx], info = paste0("age = ", data$ages[age_ndx]))
  ## since no movement check all recoveries are in release regions
  rel_plt = plot_frequency_of_tag_release_and_recoveries(data, region_key = region_key, release_ndx_to_plot = 1:5)
  tmp_rel_def = rel_plt$data %>% dplyr::group_by(release_event) %>% dplyr::summarise(releases = unique(n_releases))
  ## no mortality or movement. So tag-recoveries n should equal releases
  plt = plot_tag_recovery_obs(test_report, region_key = region_key, sex = "both", release_ndx_to_plot = 1:5)
  tmp_rec_def = plt$data %>% dplyr::filter(recovery_region == release_region) %>% dplyr::group_by(release_event, recovery_year) %>% dplyr::summarise(recoveries = sum(predicted))
  for(i in 1:nrow(tmp_rel_def))
    expect_true(all(tmp_rec_def$recoveries[tmp_rec_def$release_event %in% tmp_rel_def$release_event[i]] == tmp_rel_def$releases[i]), info = paste0("release event = ", tmp_rel_def$release_event[i]))

})


#' single-release-no-movement-with-Z
#' @description this will track taged cohorts in a model with no movement but now apply mortality
#' check Z and ageing are being applied correctly
#' We only test movement for one release event of the 5 in the model.
#'
test_that("single-release-no-movement-with-Z", {
  ## Read in mock data
  load(system.file("testdata", "MockSablefishModel.RData",package="SpatialSablefishAssessment"))

  data$model = "TagIntegratedValidate"
  ## we apply Z now
  data$apply_Z_on_tagged_fish = 1
  ## still no movement
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



  test_model <- TMB::MakeADFun(data = data,
                               parameters = parameters,
                               DLL = "SpatialSablefishAssessment_TMBExports", silent  = T)

  test_report = test_model$report()
  ## test that the same number of tagged fish are in the tagged partition as was at the beginning
  all_releases = sum(data$male_tagged_cohorts_by_age) + sum(data$female_tagged_cohorts_by_age)
  releases_at_the_end = sum(test_report$tagged_natage_f) + sum(test_report$tagged_natage_m)
  ## should be equal because no Z
  expect_false(all_releases == releases_at_the_end, label = "tagged fish at end of model should not be the same as input as we have Z now")

  ## with no movement the release and recovery periods should be the same
  final_year = max(data$years)

  ## compare numbers at age for one of the releases that were released in region 1
  init_n_age = data$male_tagged_cohorts_by_age[,1,1]
  ## age this partition
  release_year = data$years[which(data$tag_release_event_this_year == 1)][1]
  ##
  times_to_age = final_year - release_year + 1 ## plus one for the last year
  for(i in 1:times_to_age) {
    temp_tag_partition = init_n_age
    init_n_age = vector(length = length(init_n_age)) ## clear it
    init_n_age[2:length(init_n_age)] = temp_tag_partition[1:(length(init_n_age) - 1)] * test_report$S_m[1:(length(init_n_age) - 1),1,i]
    init_n_age[length(init_n_age)] = temp_tag_partition[length(init_n_age) - 1] * test_report$S_m[length(init_n_age) - 1,1,i] +  temp_tag_partition[length(init_n_age)] * test_report$S_m[length(init_n_age),1,i]
  }
  ## check the ageing and Z are being applied correctly
  for(age_ndx in 1:length(init_n_age))
    expect_true(test_report$tagged_natage_m[age_ndx,1,get_tag_release_ndx(1, data$n_years_to_retain_tagged_cohorts_for + 1, data$n_regions)] == init_n_age[age_ndx], info = paste0("age = ", data$ages[age_ndx]), label = "male")

  ## compare numbers at age for one of the releases that were released in region 1
  init_n_age = data$female_tagged_cohorts_by_age[,1,1]
  ## age this partition
  release_year = data$years[which(data$tag_release_event_this_year == 1)][1]
  ##
  times_to_age = final_year - release_year + 1 ## plus one for the last year
  for(i in 1:times_to_age) {
    temp_tag_partition = init_n_age
    init_n_age = vector(length = length(init_n_age)) ## clear it
    init_n_age[2:length(init_n_age)] = temp_tag_partition[1:(length(init_n_age) - 1)] * test_report$S_f[1:(length(init_n_age) - 1),1,i]
    init_n_age[length(init_n_age)] = temp_tag_partition[length(init_n_age) - 1] * test_report$S_f[length(init_n_age) - 1,1,i] +  temp_tag_partition[length(init_n_age)] * test_report$S_f[length(init_n_age),1,i]
  }
  ## check the ageing and Z are being applied correctly
  for(age_ndx in 1:length(init_n_age))
    expect_true(test_report$tagged_natage_f[age_ndx,1,get_tag_release_ndx(1, data$n_years_to_retain_tagged_cohorts_for + 1, data$n_regions)] == init_n_age[age_ndx], info = paste0("age = ", data$ages[age_ndx]), label = "female")

})


#' single-release-with-movement-no-Z
#' @description this will track taged cohorts in a model with movement and no mortality
#' We only test movement and Z for one release event of the 5 in the model.
#'
test_that("single-release-with-movement-no-Z", {
  ## Read in mock data
  load(system.file("testdata", "MockSablefishModel.RData",package="SpatialSablefishAssessment"))

  data$model = "TagIntegratedValidate"
  ## we apply Z now
  data$apply_Z_on_tagged_fish = 0
  ## still no movement
  data$apply_fixed_movement = 1
  ## this assumes the movment matrix in data$movement_matrix
  data$fixed_movement_matrix  = data$movement_matrix
  data$apply_fishery_tag_reporting = 0 ## all tagged fish will be recovered not a function of F
  ## turn off tag shedding and initial mortality
  data$initial_tag_induced_mortality = rep(0.0, sum(data$tag_release_event_this_year))
  data$annual_tag_shedding_rate = 0.0



  test_model <- TMB::MakeADFun(data = data,
                               parameters = parameters,
                               DLL = "SpatialSablefishAssessment_TMBExports", silent  = T)

  test_report = test_model$report()
  ## test that the same number of tagged fish are in the tagged partition as was at the beginning
  all_releases = sum(data$male_tagged_cohorts_by_age) + sum(data$female_tagged_cohorts_by_age)
  releases_at_the_end = sum(test_report$tagged_natage_f) + sum(test_report$tagged_natage_m)
  ## should be equal because no Z
  ## change to expect equal because of precision difference between C++ and R
  expect_equal(all_releases, releases_at_the_end, tolerance = 0.001, label = "tagged fish at end of model should be the same as input there is no Z")

  ## with no movement the release and recovery periods should be the same
  final_year = max(data$years)

  ## compare numbers at age for one of the releases that were released in region 1
  init_n_age = data$male_tagged_cohorts_by_age[,1,1]
  ## age this partition
  release_year = data$years[which(data$tag_release_event_this_year == 1)][1]
  ##
  times_to_age = final_year - release_year + 1 ## plus one for the last year
  n_ages = length(data$ages)
  numbers_by_age_and_region = matrix(0, nrow = n_ages,  ncol = data$n_regions)
  numbers_by_age_and_region[,1] = init_n_age
  for(i in 1:times_to_age) {
    ## ageing first
    for(reg_ndx in 1:data$n_regions) {
      temp_tag_partition = numbers_by_age_and_region[,reg_ndx]
      numbers_by_age_and_region[,reg_ndx] = rep(0, n_ages)
      numbers_by_age_and_region[2:n_ages, reg_ndx] = temp_tag_partition[1:(n_ages - 1)]
      numbers_by_age_and_region[n_ages, reg_ndx] = temp_tag_partition[n_ages - 1] +  temp_tag_partition[n_ages]

    }
    ## now movement
    numbers_by_age_and_region = numbers_by_age_and_region %*% data$movement_matrix[,,1]
  }
  ## check the ageing and Z are being applied correctly
  for(reg_ndx in 1:data$n_regions) {
    for(age_ndx in 1:length(init_n_age))
      expect_true(test_report$tagged_natage_m[age_ndx,reg_ndx,get_tag_release_ndx(1, data$n_years_to_retain_tagged_cohorts_for + 1, data$n_regions)] == numbers_by_age_and_region[age_ndx, reg_ndx], info = paste0("age = ", data$ages[age_ndx]), label = "male")
  }

  ## repeat for females
  init_n_age = data$female_tagged_cohorts_by_age[,1,1]
  ## age this partition
  release_year = data$years[which(data$tag_release_event_this_year == 1)][1]
  ##
  times_to_age = final_year - release_year + 1 ## plus one for the last year
  n_ages = length(data$ages)
  numbers_by_age_and_region = matrix(0, nrow = n_ages,  ncol = data$n_regions)
  numbers_by_age_and_region[,1] = init_n_age
  for(i in 1:times_to_age) {
    ## ageing first
    for(reg_ndx in 1:data$n_regions) {
      temp_tag_partition = numbers_by_age_and_region[,reg_ndx]
      numbers_by_age_and_region[,reg_ndx] = rep(0, n_ages)
      numbers_by_age_and_region[2:n_ages, reg_ndx] = temp_tag_partition[1:(n_ages - 1)]
      numbers_by_age_and_region[n_ages, reg_ndx] = temp_tag_partition[n_ages - 1] +  temp_tag_partition[n_ages]
    }
    ## now movement
    numbers_by_age_and_region = numbers_by_age_and_region %*% data$movement_matrix[,,1]
  }
  ## check the ageing and Z are being applied correctly
  for(reg_ndx in 1:data$n_regions) {
    for(age_ndx in 1:length(init_n_age))
      expect_true(test_report$tagged_natage_f[age_ndx,reg_ndx,get_tag_release_ndx(1, data$n_years_to_retain_tagged_cohorts_for + 1, data$n_regions)] == numbers_by_age_and_region[age_ndx, reg_ndx], info = paste0("age = ", data$ages[age_ndx]), label = "female")
  }
})



#' single-release-with-movement-and-Z
#' @description this will track taged cohorts in a model with movement and no mortaltity
#'
test_that("single-release-with-movement-and-Z", {
  ## Read in mock data
  load(system.file("testdata", "MockSablefishModel.RData",package="SpatialSablefishAssessment"))

  data$model = "TagIntegratedValidate"
  ## we apply Z now
  data$apply_Z_on_tagged_fish = 1
  ## still no movement
  data$apply_fixed_movement = 1
  ## this assumes the movment matrix in data$movement_matrix
  data$fixed_movement_matrix  = data$movement_matrix
  data$apply_fishery_tag_reporting = 0 ## all tagged fish will be recovered not a function of F
  ## turn off tag shedding and initial mortality
  data$initial_tag_induced_mortality = rep(0.0, sum(data$tag_release_event_this_year))
  data$annual_tag_shedding_rate = 0.0



  test_model <- TMB::MakeADFun(data = data,
                               parameters = parameters,
                               DLL = "SpatialSablefishAssessment_TMBExports", silent  = T)

  test_report = test_model$report()
  ## with no movement the release and recovery periods should be the same
  final_year = max(data$years)

  ## compare numbers at age for one of the releases that were released in region 1
  init_n_age = data$male_tagged_cohorts_by_age[,1,1]
  ## age this partition
  release_year = data$years[which(data$tag_release_event_this_year == 1)][1]
  ##
  times_to_age = final_year - release_year + 1 ## plus one for the last year
  n_ages = length(data$ages)
  numbers_by_age_and_region = matrix(0, nrow = n_ages,  ncol = data$n_regions)
  numbers_by_age_and_region[,1] = init_n_age
  for(i in 1:times_to_age) {
    ## ageing first
    for(reg_ndx in 1:data$n_regions) {
      temp_tag_partition = numbers_by_age_and_region[,reg_ndx]
      numbers_by_age_and_region[,reg_ndx] = rep(0, n_ages)
      numbers_by_age_and_region[2:n_ages, reg_ndx] = temp_tag_partition[1:(n_ages - 1)] * test_report$S_m[1:(n_ages - 1),reg_ndx,i]
      numbers_by_age_and_region[n_ages, reg_ndx] = temp_tag_partition[n_ages - 1] * test_report$S_m[n_ages - 1,reg_ndx,i] +  temp_tag_partition[n_ages] *  test_report$S_m[n_ages,reg_ndx,i]

    }
    ## now movement
    numbers_by_age_and_region = numbers_by_age_and_region %*% data$movement_matrix[,,1]
  }
  ## check the ageing and Z are being applied correctly
  for(reg_ndx in 1:data$n_regions) {
    for(age_ndx in 1:length(init_n_age))
      expect_true(test_report$tagged_natage_m[age_ndx,reg_ndx,get_tag_release_ndx(1, data$n_years_to_retain_tagged_cohorts_for + 1, data$n_regions)] == numbers_by_age_and_region[age_ndx, reg_ndx], info = paste0("age = ", data$ages[age_ndx]), label = "male")
  }

  ## repeat for females
  init_n_age = data$female_tagged_cohorts_by_age[,1,1]
  ## age this partition
  release_year = data$years[which(data$tag_release_event_this_year == 1)][1]
  ##
  times_to_age = final_year - release_year + 1 ## plus one for the last year
  n_ages = length(data$ages)
  numbers_by_age_and_region = matrix(0, nrow = n_ages,  ncol = data$n_regions)
  numbers_by_age_and_region[,1] = init_n_age
  for(i in 1:times_to_age) {
    ## ageing first
    for(reg_ndx in 1:data$n_regions) {
      temp_tag_partition = numbers_by_age_and_region[,reg_ndx]
      numbers_by_age_and_region[,reg_ndx] = rep(0, n_ages)
      numbers_by_age_and_region[2:n_ages, reg_ndx] = temp_tag_partition[1:(n_ages - 1)] * test_report$S_f[1:(n_ages - 1),reg_ndx,i]
      numbers_by_age_and_region[n_ages, reg_ndx] = temp_tag_partition[n_ages - 1] * test_report$S_f[n_ages - 1,reg_ndx,i] +  temp_tag_partition[n_ages] *  test_report$S_f[n_ages,reg_ndx,i]
    }
    ## now movement
    numbers_by_age_and_region = numbers_by_age_and_region %*% data$movement_matrix[,,1]
  }
  ## check the ageing and Z are being applied correctly
  for(reg_ndx in 1:data$n_regions) {
    for(age_ndx in 1:length(init_n_age))
      expect_true(test_report$tagged_natage_f[age_ndx,reg_ndx,get_tag_release_ndx(1, data$n_years_to_retain_tagged_cohorts_for + 1, data$n_regions)] == numbers_by_age_and_region[age_ndx, reg_ndx], info = paste0("age = ", data$ages[age_ndx]), label = "female")
  }
})



#' single-release-initial-mortality
#' @description validate the initial mortality is applied correctly
#'
test_that("single-release-initial-mortality", {
  ## Read in mock data
  load(system.file("testdata", "MockSablefishModel.RData",package="SpatialSablefishAssessment"))

  data$model = "TagIntegratedValidate"
  ## no Z or movement
  data$apply_Z_on_tagged_fish = 0
  data$apply_fixed_movement = 1
  ## this assumes no movement
  data$apply_fishery_tag_reporting = 0 ## all tagged fish will be recovered not a function of F
  ## turn off tag shedding and initial mortality
  data$initial_tag_induced_mortality = rep(0.1, sum(data$tag_release_event_this_year))
  data$annual_tag_shedding_rate = 0.0

  test_model <- TMB::MakeADFun(data = data,
                               parameters = parameters,
                               DLL = "SpatialSablefishAssessment_TMBExports", silent  = T)

  test_report = test_model$report()
  ## test that the same number of tagged fish are in the tagged partition as was at the beginning
  all_releases = sum(data$male_tagged_cohorts_by_age) + sum(data$female_tagged_cohorts_by_age)
  releases_at_the_end = sum(test_report$tagged_natage_f) + sum(test_report$tagged_natage_m)
  ## should be equal because no Z
  expect_equal((all_releases * exp(-unique(data$initial_tag_induced_mortality))), releases_at_the_end, tolerance = 0.001, label = "tagged fish at end of model not as expected")

})


#' single-release-tag-shedding
#' @description validate the annual tag-shedding rate is applied correctly
#'
test_that("single-release-tag-shedding", {
  ## Read in mock data
  load(system.file("testdata", "MockSablefishModel.RData",package="SpatialSablefishAssessment"))

  data$model = "TagIntegratedValidate"
  ## no Z or movement
  data$apply_Z_on_tagged_fish = 0
  data$apply_fixed_movement = 1
  ## this assumes no movement
  data$apply_fishery_tag_reporting = 0 ## all tagged fish will be recovered not a function of F
  ## turn off tag shedding and initial mortality
  data$initial_tag_induced_mortality = rep(0.0, sum(data$tag_release_event_this_year))
  data$annual_tag_shedding_rate = 0.05 ## approx 5%

  test_model <- TMB::MakeADFun(data = data,
                               parameters = parameters,
                               DLL = "SpatialSablefishAssessment_TMBExports", silent  = T)

  test_report = test_model$report()

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
  for(i in 1:times_to_age) {
    temp_tag_partition = init_n_age
    init_n_age = vector(length = length(init_n_age)) ## clear it
    init_n_age[2:length(init_n_age)] = temp_tag_partition[1:(length(init_n_age) - 1)]
    init_n_age[length(init_n_age)] = temp_tag_partition[length(init_n_age) - 1] +  temp_tag_partition[length(init_n_age)]
    init_n_age = init_n_age * exp(-0.05) ## annual shedding rate
  }
  ## check the ageing works as expected
  for(age_ndx in 1:length(init_n_age))
    expect_true(test_report$tagged_natage_m[age_ndx,1,get_tag_release_ndx(1, data$n_years_to_retain_tagged_cohorts_for + 1, data$n_regions)] == init_n_age[age_ndx], info = paste0("age = ", data$ages[age_ndx]))
  ## repeat for the females
  ## compare numbers at age for one of the releases that were released in region 1
  init_n_age = data$female_tagged_cohorts_by_age[,1,1]
  ## age this partition
  release_year = data$years[which(data$tag_release_event_this_year == 1)][1]
  ##
  times_to_age = final_year - release_year + 1 ## plus one for the last year
  for(i in 1:times_to_age) {
    temp_tag_partition = init_n_age
    init_n_age = vector(length = length(init_n_age)) ## clear it
    init_n_age[2:length(init_n_age)] = temp_tag_partition[1:(length(init_n_age) - 1)]
    init_n_age[length(init_n_age)] = temp_tag_partition[length(init_n_age) - 1] +  temp_tag_partition[length(init_n_age)]
    init_n_age = init_n_age * exp(-0.05) ## annual shedding rate
  }
  for(age_ndx in 1:length(init_n_age))
    expect_true(test_report$tagged_natage_f[age_ndx,1,get_tag_release_ndx(1, data$n_years_to_retain_tagged_cohorts_for + 1, data$n_regions)] == init_n_age[age_ndx], info = paste0("age = ", data$ages[age_ndx]))


})




#' test-get_tag_release_ndx-method
#' @description A fundamental utility function is get_tag_release_ndx
#'
test_that("test-get_tag_release_ndx-method", {
  n_regions = 5
  n_years_to_retain_tagged_cohorts_for = 10

  result = 1:50
  counter = 1;
  for(tag_partition_ndx in 1:n_years_to_retain_tagged_cohorts_for) {
    for(reg_ndx in 1:n_regions) {
      ndx = get_tag_release_ndx(reg_ndx, tag_partition_ndx, n_regions)
      expect_equal(ndx, result[counter])
      counter = counter + 1
    }
  }
})

