#
# A set of functions that are used to check data and parameters are consistent. The TMB model does no sanity checks for efficiency
# reasons. The R functions in this script are responsible for checking the TMB model wont run out of memory and crash.
#

#'
#' get_max_sel_pars
#' @param sel_type_for_blocks vector of selectivity types
#' @return integer specifying the number of parameters for this selectivity parameter
#' @export

get_max_sel_pars = function(sel_type_for_blocks) {
  # sel-types see C++ function inst/include/SetupSelectivities.hpp for types
  # 0 = logistic expected 2 parameters
  # 1 = Double normal expected 2 parameters
  # 2 = power functions 1 parameters
  # 3 = alternative logistic
  # 4 = exponentail
  # 5 = Three parameter double normal
  if(any(sel_type_for_blocks %in% c(0,1, 3))) {
    return(2)
  } else if (any(sel_type_for_blocks %in% c(5))) {
    ## three parameter double normal
    return(3)
  }
  return(1)
}

#'
#' get_number_of_sel_pars
#' @param sel_type integer of selecvitiy types
#' @return integer specifying the number of parameters for this selectivity parameter
#' @export

get_number_of_sel_pars = function(sel_type) {
  # sel-types see C++ function inst/include/SetupSelectivities.hpp for types
  # 0 = logistic expected 2 parameters
  # 1 = Double normal expected 2 parameters
  # 2 = power functions 1 parameters
  # 3 = alternative logistic
  # 4 = exponentail
  # 5 = Three parameter double normal
  if(sel_type == 0) {
    return(2)
  } else if(sel_type == 1) {
    return(2)
  } else if (sel_type == 2) {
    return(1)
  } else if (sel_type == 3) {
    return(2)
  } else if (sel_type == 4) {
    return(1)
  } else if (sel_type == 5) {
    return(3)
  }
  return(NA)
}

#' check_dim checks the dimension of an array are consistent with expected values
#' @param array a matrix or array that we are checking
#' @param expected_dimensions a vector with the expected number of dimensions that the array should have
#' @return a list
#' @export
#'
check_dim = function(array, expected_dimensions) {
  array_dim = dim(array)
  if(length(array_dim) != length(expected_dimensions))
    return(list(result = FALSE, message = paste0("different number of dimensions. Array had ", length(array_dim), " we expected ", length(expected_dimensions))))

  for(i in 1:length(array_dim)) {
    if(array_dim[i] != expected_dimensions[i])
      return(list(result = FALSE, message = paste0("different number of elements for dimension ", i, ". Array had ", array_dim[i], " we expected ", expected_dimensions[i])))
  }

  return(list(result = TRUE, message = "array has expected dimension"))
}
#' check_length checks the length of vector with expected values
#' @param vector_to_check a matrix or array that we are checking
#' @param n_elements the expected length of vectors
#' @return a list
#' @export
#'
check_length = function(vector_to_check, n_elements) {
  if(length(vector_to_check) != n_elements)
    return(list(result = FALSE, message = paste0("different number of elements. Vector had ", length(vector_to_check), " we expected ", n_elements)))

  return(list(result = TRUE, message = "vector has expected number of elements"))
}

#' validate_input_data_and_parameters given a list of data and parameters will the TMB model crash?
#' @param data a list of data inputs for the model
#' @param parameters a list of parameter values for the models
#' @return string telling you if there are problems or if you are good to go
#' @export
validate_input_data_and_parameters = function(data, parameters) {
  if(!data$model %in% c("TagIntegrated", "TagIntegratedValidate","TagIntegratedAgeBasedMovement", "Assessment"))
    return("This validate function is currently only written for models: TagIntegrated, TagIntegratedValidate, TagIntegratedAgeBasedMovement,Assessment")

  n_ages = length(data$ages)
  n_length_bins = length(data$length_bins)
  n_projyears = length(data$years) +  data$n_projections_years
  n_years = length(data$years)
  n_regions = data$n_regions
  n_surveys = data$n_surveys
  projyears = min(data$years):(max(data$years) + data$n_projections_years)

  if(!data$do_projection %in% c(0,1))
    stop("data$do_projection must be either 0 or 1")
  if(min(data$ages) < 1)
    stop("data$ages must be greater than or equal to 1")

  ## check years and ages are consecutive
  for(age_ndx in 1:(length(data$ages) - 1)) {
    if((data$ages[age_ndx] + 1) != data$ages[age_ndx + 1])
      stop("data$ages needs to be strictly increasing and consecutive numbers")
  }
  for(yr_ndx in 1:(length(data$years) - 1)) {
    if((data$years[yr_ndx] + 1) != data$years[yr_ndx + 1])
      stop("data$years needs to be strictly increasing and consecutive numbers")
  }

  ## check model inputs dimensions
  check = check_dim(data$M, c(n_ages, n_projyears))
  if(!check$result)
    return(paste0("M: ", check$message))
  ## maturity
  check = check_dim(data$maturity, c(n_ages, n_projyears))
  if(!check$result)
    return(paste0("maturity: ", check$message))
  ## mean weight male
  check = check_dim(data$male_mean_weight_by_age, c(n_ages, n_projyears))
  if(!check$result)
    return(paste0("male_mean_weight_by_age: ", check$message))
  ## mean weight female
  check = check_dim(data$female_mean_weight_by_age, c(n_ages, n_projyears))
  if(!check$result)
    return(paste0("female_mean_weight_by_age: ", check$message))
  ## male_age_length_transition
  check = check_dim(data$male_age_length_transition, c(n_ages, n_length_bins, n_projyears))
  if(!check$result)
    return(paste0("male_age_length_transition: ", check$message))
  ## female_age_length_transition
  check = check_dim(data$female_age_length_transition, c(n_ages, n_length_bins, n_projyears))
  if(!check$result)
    return(paste0("female_age_length_transition: ", check$message))
  ## spawning_time_proportion
  check = check_length(data$spawning_time_proportion, n_projyears)
  if(!check$result)
    return(paste0("spawning_time_proportion: ", check$message))
  ## ageing_error_matrix
  check = check_dim(data$ageing_error_matrix, c(n_ages, n_ages))
  if(!check$result)
    return(paste0("ageing_error_matrix: ", check$message))


  ## Assessment model specialties
  if(data$model == "Assessment") {
    ##
    ## parameters Containers
    ##
    ## ln_rec_dev
    check = check_length(parameters$ln_rec_dev, n_years)
    if(!check$result)
      return(paste0("parameters$ln_rec_dev: ", check$message))
    ## ln_ll_F_devs
    check = check_length(parameters$ln_ll_F_devs, n_years)
    if(!check$result)
      return(paste0("parameters$ln_ll_F_devs: ", check$message))
    ## ln_trwl_F_devs
    check = check_length(parameters$ln_trwl_F_devs, n_years)
    if(!check$result)
      return(paste0("parameters$ln_trwl_F_devs: ", check$message))
    ##
    ## Data containers
    ##
    ## Historical Japanese fishery age-length transition matrix
    check = check_dim(data$historical_jap_fishery_ll_age_length_transition, c(n_ages, n_length_bins))
    if(!check$result)
      return(paste0("historical_jap_fishery_ll_age_length_transition: ", check$message))
    ## proportion_male
    check = check_length(data$proportion_male, n_projyears)
    if(!check$result)
      return(paste0("proportion_male: ", check$message))
    ## proportion_male2
    check = check_length(data$proportion_male2, n_projyears)
    if(!check$result)
      return(paste0("proportion_male2: ", check$message))
    ## ll_fishery_catch
    check = check_length(data$ll_fishery_catch, n_projyears)
    if(!check$result)
      return(paste0("ll_fishery_catch: ", check$message))
    ## trwl_fishery_catch
    check = check_length(data$trwl_fishery_catch, n_projyears)
    if(!check$result)
      return(paste0("trwl_fishery_catch: ", check$message))
    ## ll_sel_by_year_indicator
    check = check_length(data$ll_sel_by_year_indicator, n_projyears)
    if(!check$result)
      return(paste0("ll_sel_by_year_indicator: ", check$message))
    ## trwl_sel_by_year_indicator
    check = check_length(data$trwl_sel_by_year_indicator, n_projyears)
    if(!check$result)
      return(paste0("trwl_sel_by_year_indicator: ", check$message))
    ## srv_dom_ll_sel_by_year_indicator
    check = check_length(data$srv_dom_ll_sel_by_year_indicator, n_projyears)
    if(!check$result)
      return(paste0("srv_dom_ll_sel_by_year_indicator: ", check$message))
    ## srv_dom_ll_q_by_year_indicator
    check = check_length(data$srv_dom_ll_q_by_year_indicator, n_projyears)
    if(!check$result)
      return(paste0("srv_dom_ll_q_by_year_indicator: ", check$message))
    ## srv_jap_ll_sel_by_year_indicator
    check = check_length(data$srv_jap_ll_sel_by_year_indicator, n_projyears)
    if(!check$result)
      return(paste0("srv_jap_ll_sel_by_year_indicator: ", check$message))
    ## srv_jap_ll_q_by_year_indicator
    check = check_length(data$srv_jap_ll_q_by_year_indicator, n_projyears)
    if(!check$result)
      return(paste0("srv_jap_ll_q_by_year_indicator: ", check$message))
    ## srv_nmfs_trwl_sel_by_year_indicator
    check = check_length(data$srv_nmfs_trwl_sel_by_year_indicator, n_projyears)
    if(!check$result)
      return(paste0("srv_nmfs_trwl_sel_by_year_indicator: ", check$message))
    ## srv_nmfs_trwl_q_by_year_indicator
    check = check_length(data$srv_nmfs_trwl_q_by_year_indicator, n_projyears)
    if(!check$result)
      return(paste0("srv_nmfs_trwl_q_by_year_indicator: ", check$message))
    ## srv_jap_fishery_ll_sel_by_year_indicator
    check = check_length(data$srv_jap_fishery_ll_sel_by_year_indicator, n_projyears)
    if(!check$result)
      return(paste0("srv_jap_fishery_ll_sel_by_year_indicator: ", check$message))
    ## srv_jap_fishery_ll_q_by_year_indicator
    check = check_length(data$srv_jap_fishery_ll_q_by_year_indicator, n_projyears)
    if(!check$result)
      return(paste0("srv_jap_fishery_ll_q_by_year_indicator: ", check$message))
    ## ll_cpue_q_by_year_indicator
    check = check_length(data$ll_cpue_q_by_year_indicator, n_projyears)
    if(!check$result)
      return(paste0("ll_cpue_q_by_year_indicator: ", check$message))
    ## ll_catchatage_indicator
    check = check_length(data$ll_catchatage_indicator, n_projyears)
    if(!check$result)
      return(paste0("ll_catchatage_indicator: ", check$message))
    ## obs_ll_catchatage
    check = check_dim(data$obs_ll_catchatage, c(n_ages, sum(data$ll_catchatage_indicator)))
    if(!check$result)
      return(paste0("obs_ll_catchatage: ", check$message))
    ## ll_catchatlgth_indicator
    check = check_length(data$ll_catchatlgth_indicator, n_projyears)
    if(!check$result)
      return(paste0("ll_catchatlgth_indicator: ", check$message))
    ## obs_ll_catchatlgth_m
    if(sum(data$ll_catchatlgth_indicator) != 0) { # only check if observations are supplied
      check = check_dim(data$obs_ll_catchatlgth_m, c(n_length_bins, sum(data$ll_catchatlgth_indicator)))
      if(!check$result)
        return(paste0("obs_ll_catchatlgth_m: ", check$message))
      ## obs_ll_catchatlgth_f
      check = check_dim(data$obs_ll_catchatlgth_f, c(n_length_bins, sum(data$ll_catchatlgth_indicator)))
      if(!check$result)
        return(paste0("obs_ll_catchatlgth_f: ", check$message))
    }
    ## trwl_catchatlgth_indicator
    check = check_length(data$trwl_catchatlgth_indicator, n_projyears)
    if(!check$result)
      return(paste0("trwl_catchatlgth_indicator: ", check$message))
    ## obs_trwl_catchatlgth_m
    if(sum(data$trwl_catchatlgth_indicator) != 0) { # only check if observations are supplied
      check = check_dim(data$obs_trwl_catchatlgth_m, c(n_length_bins, sum(data$trwl_catchatlgth_indicator)))
      if(!check$result)
        return(paste0("obs_trwl_catchatlgth_m: ", check$message))
      ## obs_trwl_catchatlgth_f
      check = check_dim(data$obs_trwl_catchatlgth_f, c(n_length_bins, sum(data$trwl_catchatlgth_indicator)))
      if(!check$result)
        return(paste0("obs_trwl_catchatlgth_f: ", check$message))
    }
    ## srv_dom_ll_bio_indicator
    check = check_length(data$srv_dom_ll_bio_indicator, n_projyears)
    if(!check$result)
      return(paste0("srv_dom_ll_bio_indicator: ", check$message))
    ## obs_dom_ll_bio
    if(sum(data$srv_dom_ll_bio_indicator) != 0) { # only check if observations are supplied
      check = check_length(data$obs_dom_ll_bio, sum(data$srv_dom_ll_bio_indicator))
      if(!check$result)
        return(paste0("obs_dom_ll_bio: ", check$message))
      ## se_dom_ll_bio
      check = check_length(data$se_dom_ll_bio, sum(data$srv_dom_ll_bio_indicator))
      if(!check$result)
        return(paste0("se_dom_ll_bio: ", check$message))
    }
    ## srv_jap_ll_bio_indicator
    check = check_length(data$srv_jap_ll_bio_indicator, n_projyears)
    if(!check$result)
      return(paste0("srv_jap_ll_bio_indicator: ", check$message))
    ## obs_jap_ll_bio
    if(sum(data$srv_jap_ll_bio_indicator) != 0) { # only check if observations are supplied
      check = check_length(data$obs_jap_ll_bio, sum(data$srv_jap_ll_bio_indicator))
      if(!check$result)
        return(paste0("obs_jap_ll_bio: ", check$message))
      ## se_jap_ll_bio
      check = check_length(data$se_jap_ll_bio, sum(data$srv_jap_ll_bio_indicator))
      if(!check$result)
        return(paste0("se_jap_ll_bio: ", check$message))
    }
    ## srv_nmfs_trwl_bio_indicator
    check = check_length(data$srv_nmfs_trwl_bio_indicator, n_projyears)
    if(!check$result)
      return(paste0("srv_nmfs_trwl_bio_indicator: ", check$message))
    ## obs_nmfs_trwl_bio
    if(sum(data$srv_nmfs_trwl_bio_indicator) != 0) { # only check if observations are supplied
      check = check_length(data$obs_nmfs_trwl_bio, sum(data$srv_nmfs_trwl_bio_indicator))
      if(!check$result)
        return(paste0("obs_nmfs_trwl_bio: ", check$message))
      ## se_nmfs_trwl_bio
      check = check_length(data$se_nmfs_trwl_bio, sum(data$srv_nmfs_trwl_bio_indicator))
      if(!check$result)
        return(paste0("se_nmfs_trwl_bio: ", check$message))
    }
    ## ll_cpue_indicator
    check = check_length(data$ll_cpue_indicator, n_projyears)
    if(!check$result)
      return(paste0("ll_cpue_indicator: ", check$message))
    ## obs_ll_cpue
    if(sum(data$ll_cpue_indicator) != 0) { # only check if observations are supplied
      check = check_length(data$obs_ll_cpue, sum(data$ll_cpue_indicator))
      if(!check$result)
        return(paste0("obs_ll_cpue: ", check$message))
      ## se_ll_cpue
      check = check_length(data$se_ll_cpue, sum(data$ll_cpue_indicator))
      if(!check$result)
        return(paste0("se_ll_cpue: ", check$message))
    }
    ## srv_jap_fishery_ll_bio_indicator
    check = check_length(data$srv_jap_fishery_ll_bio_indicator, n_projyears)
    if(!check$result)
      return(paste0("srv_jap_fishery_ll_bio_indicator: ", check$message))
    ## obs_jap_fishery_ll_bio
    if(sum(data$srv_jap_fishery_ll_bio_indicator) != 0) { # only check if observations are supplied
      check = check_length(data$obs_jap_fishery_ll_bio, sum(data$srv_jap_fishery_ll_bio_indicator))
      if(!check$result)
        return(paste0("obs_jap_fishery_ll_bio: ", check$message))
      ## se_jap_fishery_ll_bio
      check = check_length(data$se_jap_fishery_ll_bio, sum(data$srv_jap_fishery_ll_bio_indicator))
      if(!check$result)
        return(paste0("se_jap_fishery_ll_bio: ", check$message))
    }
    ## srv_dom_ll_age_indicator
    check = check_length(data$srv_dom_ll_age_indicator, n_projyears)
    if(!check$result)
      return(paste0("srv_dom_ll_age_indicator: ", check$message))
    ## obs_srv_dom_ll_age
    if(sum(data$srv_dom_ll_age_indicator) != 0) { # only check if observations are supplied
      check = check_dim(data$obs_srv_dom_ll_age, c(n_ages, sum(data$srv_dom_ll_age_indicator)))
      if(!check$result)
        return(paste0("obs_srv_dom_ll_age: ", check$message))
    }
    ## srv_dom_ll_lgth_indicator
    check = check_length(data$srv_dom_ll_lgth_indicator, n_projyears)
    if(!check$result)
      return(paste0("srv_dom_ll_lgth_indicator: ", check$message))
    ## obs_srv_dom_ll_lgth_m
    if(sum(data$srv_dom_ll_lgth_indicator) != 0) { # only check if observations are supplied
      check = check_dim(data$obs_srv_dom_ll_lgth_m, c(n_length_bins, sum(data$srv_dom_ll_lgth_indicator)))
      if(!check$result)
        return(paste0("obs_srv_dom_ll_lgth_m: ", check$message))
      ## obs_srv_dom_ll_lgth_f
      check = check_dim(data$obs_srv_dom_ll_lgth_f, c(n_length_bins, sum(data$srv_dom_ll_lgth_indicator)))
      if(!check$result)
        return(paste0("obs_srv_dom_ll_lgth_f: ", check$message))
    }
    ## srv_jap_ll_age_indicator
    check = check_length(data$srv_jap_ll_age_indicator, n_projyears)
    if(!check$result)
      return(paste0("srv_jap_ll_age_indicator: ", check$message))
    ## obs_srv_jap_ll_age
    if(sum(data$srv_jap_ll_age_indicator) != 0) { # only check if observations are supplied
      check = check_dim(data$obs_srv_jap_ll_age, c(n_ages, sum(data$srv_jap_ll_age_indicator)))
      if(!check$result)
        return(paste0("obs_srv_jap_ll_age: ", check$message))
    }
    ## srv_jap_ll_lgth_indicator
    check = check_length(data$srv_jap_ll_lgth_indicator, n_projyears)
    if(!check$result)
      return(paste0("srv_jap_ll_lgth_indicator: ", check$message))
    ## obs_srv_jap_ll_lgth_m
    if(sum(data$srv_jap_ll_lgth_indicator) != 0) { # only check if observations are supplied
      check = check_dim(data$obs_srv_jap_ll_lgth_m, c(n_length_bins, sum(data$srv_jap_ll_lgth_indicator)))
      if(!check$result)
        return(paste0("obs_srv_jap_ll_lgth_m: ", check$message))
      ## obs_srv_jap_ll_lgth_f
      check = check_dim(data$obs_srv_jap_ll_lgth_f, c(n_length_bins, sum(data$srv_jap_ll_lgth_indicator)))
      if(!check$result)
        return(paste0("obs_srv_jap_ll_lgth_f: ", check$message))
    }
    ## srv_jap_fishery_ll_lgth_indicator
    check = check_length(data$srv_jap_fishery_ll_lgth_indicator, n_projyears)
    if(!check$result)
      return(paste0("srv_jap_fishery_ll_lgth_indicator: ", check$message))
    ## obs_srv_jap_fishery_ll_lgth
    if(sum(data$srv_jap_fishery_ll_lgth_indicator) != 0) { # only check if observations are supplied
      check = check_dim(data$obs_srv_jap_fishery_ll_lgth, c(n_length_bins, sum(data$srv_jap_fishery_ll_lgth_indicator)))
      if(!check$result)
        return(paste0("obs_srv_jap_fishery_ll_lgth: ", check$message))
    }
    ## srv_nmfs_trwl_age_indicator
    check = check_length(data$srv_nmfs_trwl_age_indicator, n_projyears)
    if(!check$result)
      return(paste0("srv_nmfs_trwl_age_indicator: ", check$message))
    ## obs_srv_dom_ll_age
    if(sum(data$srv_nmfs_trwl_age_indicator) != 0) { # only check if observations are supplied
      check = check_dim(data$obs_srv_nmfs_trwl_age, c(n_ages, sum(data$srv_nmfs_trwl_age_indicator)))
      if(!check$result)
        return(paste0("obs_srv_nmfs_trwl_age: ", check$message))
    }
    ## srv_nmfs_trwl_lgth_indicator
    check = check_length(data$srv_nmfs_trwl_lgth_indicator, n_projyears)
    if(!check$result)
      return(paste0("srv_nmfs_trwl_lgth_indicator: ", check$message))
    ## obs_srv_nmfs_trwl_lgth_m
    if(sum(data$srv_nmfs_trwl_lgth_indicator) != 0) { # only check if observations are supplied
      check = check_dim(data$obs_srv_nmfs_trwl_lgth_m, c(n_length_bins, sum(data$srv_nmfs_trwl_lgth_indicator)))
      if(!check$result)
        return(paste0("obs_srv_nmfs_trwl_lgth_m: ", check$message))
      ## obs_srv_nmfs_trwl_lgth_f
      check = check_dim(data$obs_srv_nmfs_trwl_lgth_f, c(n_length_bins, sum(data$srv_nmfs_trwl_lgth_indicator)))
      if(!check$result)
        return(paste0("obs_srv_nmfs_trwl_lgth_f: ", check$message))
    }
    ## check prior inputs
    check = check_length(data$mu_cpue_q, length(parameters$ln_ll_cpue_q))
    if(!check$result)
      return(paste0("mu_cpue_q: ", check$message, ". this is the length of 'ln_ll_cpue_q'"))
    check = check_length(data$sd_cpue_q, length(parameters$ln_ll_cpue_q))
    if(!check$result)
      return(paste0("sd_cpue_q: ", check$message, ". this is the length of 'ln_ll_cpue_q'"))
    check = check_length(data$mu_srv_jap_fishery_ll_q, length(parameters$ln_srv_jap_fishery_ll_q))
    if(!check$result)
      return(paste0("mu_srv_jap_fishery_ll_q: ", check$message, ". this is the length of 'ln_srv_jap_fishery_ll_q'"))
    check = check_length(data$sd_srv_jap_fishery_ll_q, length(parameters$ln_srv_jap_fishery_ll_q))
    if(!check$result)
      return(paste0("sd_srv_jap_fishery_ll_q: ", check$message, ". this is the length of 'ln_srv_jap_fishery_ll_q'"))
    check = check_length(data$mu_srv_nmfs_trwl_q, length(parameters$ln_srv_nmfs_trwl_q))
    if(!check$result)
      return(paste0("mu_srv_nmfs_trwl_q: ", check$message, ". this is the length of 'ln_srv_nmfs_trwl_q'"))
    check = check_length(data$sd_srv_nmfs_trwl_q, length(parameters$ln_srv_nmfs_trwl_q))
    if(!check$result)
      return(paste0("sd_srv_nmfs_trwl_q: ", check$message, ". this is the length of 'ln_srv_nmfs_trwl_q'"))
    check = check_length(data$mu_srv_jap_ll_q, length(parameters$ln_srv_jap_ll_q))
    if(!check$result)
      return(paste0("mu_srv_jap_ll_q: ", check$message, ". this is the length of 'ln_srv_jap_ll_q'"))
    check = check_length(data$sd_srv_jap_ll_q, length(parameters$ln_srv_jap_ll_q))
    if(!check$result)
      return(paste0("sd_srv_jap_ll_q: ", check$message, ". this is the length of 'ln_srv_jap_ll_q'"))
    check = check_length(data$mu_srv_dom_ll_q, length(parameters$ln_srv_dom_ll_q))
    if(!check$result)
      return(paste0("mu_srv_dom_ll_q: ", check$message, ". this is the length of 'ln_srv_dom_ll_q'"))
    check = check_length(data$sd_srv_dom_ll_q, length(parameters$ln_srv_dom_ll_q))
    if(!check$result)
      return(paste0("sd_srv_dom_ll_q: ", check$message, ". this is the length of 'ln_srv_dom_ll_q'"))

  } else {
    if(data$do_projection == 1) {
      ## check projection variables
      if(data$n_projections_years <= 0)
        stop("data$n_projections_years must be greater than max(data$years)")

      ## future fishing containers
      check = check_dim(data$future_fishing_inputs_fixed, c(n_regions, data$n_projections_years))
      if(!check$result)
        return(paste0("future_fishing_inputs_fixed: ", check$message))
      check = check_dim(data$future_fishing_inputs_trwl, c(n_regions, data$n_projections_years))
      if(!check$result)
        return(paste0("future_fishing_inputs_trwl: ", check$message))

    }

    if(data$model == "TagIntegratedAgeBasedMovement") {
      ## fixed_movement_matrix_old
      check = check_dim(data$fixed_movement_matrix_old, c(n_regions, n_regions))
      if(!check$result)
        return(paste0("fixed_movement_matrix_old: ", check$message))
      ## fixed_movement_matrix_young
      check = check_dim(data$fixed_movement_matrix_young, c(n_regions, n_regions))
      if(!check$result)
        return(paste0("fixed_movement_matrix_young: ", check$message))

    } else {
      ## fixed_movement_matrix
      check = check_dim(data$fixed_movement_matrix, c(n_regions, n_regions))
      if(!check$result)
        return(paste0("fixed_movement_matrix: ", check$message))
    }

    ## fixed_fishery_catch
    check = check_dim(data$fixed_fishery_catch, c(n_regions, n_years))
    if(!check$result)
      return(paste0("fixed_fishery_catch: ", check$message))
    ## trwl_fishery_catch
    check = check_dim(data$trwl_fishery_catch, c(n_regions, n_years))
    if(!check$result)
      return(paste0("trwl_fishery_catch: ", check$message))
    ## fixed_sel_by_year_indicator
    check = check_length(data$fixed_sel_by_year_indicator, n_projyears)
    if(!check$result)
      return(paste0("fixed_sel_by_year_indicator: ", check$message))
    ## trwl_sel_by_year_indicator
    check = check_length(data$trwl_sel_by_year_indicator, n_projyears)
    if(!check$result)
      return(paste0("trwl_sel_by_year_indicator: ", check$message))

    ## srv_sel_by_year_indicator
    check = check_dim(data$srv_sel_by_year_indicator, c(n_projyears, n_surveys))
    if(!check$result)
      return(paste0("srv_sel_by_year_indicator: ", check$message))

    n_survey_sel_time_blocks = max(apply(data$srv_sel_by_year_indicator, MARGIN = 2, function(x){length(unique(x))}))
    ## srv_sel_type
    check = check_dim(data$srv_sel_type, c(n_survey_sel_time_blocks, n_surveys))
    if(!check$result)
      return(paste0("srv_sel_type: ", check$message))
    if(length(unique(as.vector(data$srv_sel_type))) > 1)
      warning("'srv_sel_type', has different selectivity types, you will need to be careful how you turn off values (or map) values for parameters$ln_srv_sel_pars")

    ## srv_obs_is_abundance
    check = check_length(data$srv_obs_is_abundance, n_surveys)
    if(!check$result)
      return(paste0("srv_obs_is_abundance: ", check$message))
    ## q_is_nuisance
    check = check_length(data$q_is_nuisance, n_surveys)
    if(!check$result)
      return(paste0("q_is_nuisance: ", check$message))

    ## tag release
    ## tag_release_event_this_year
    check = check_length(data$tag_release_event_this_year, n_years)
    if(!check$result)
      return(paste0("tag_release_event_this_year: ", check$message))

    n_tag_releases = sum(data$tag_release_event_this_year)
    if(n_tag_releases > 0) {
      ## male_tagged_cohorts_by_age
      check = check_dim(data$male_tagged_cohorts_by_age, c(n_ages, n_regions, n_tag_releases))
      if(!check$result)
        return(paste0("male_tagged_cohorts_by_age: ", check$message))
      ## female_tagged_cohorts_by_age
      check = check_dim(data$female_tagged_cohorts_by_age, c(n_ages, n_regions, n_tag_releases))
      if(!check$result)
        return(paste0("female_tagged_cohorts_by_age: ", check$message))
      ## initial_tag_induced_mortality
      check = check_length(data$initial_tag_induced_mortality, n_tag_releases)
      if(!check$result)
        return(paste0("initial_tag_induced_mortality: ", check$message))
    }

    ## fixed_catchatage_indicator
    check = check_dim(data$fixed_catchatage_indicator, c(n_regions, n_years))
    if(!check$result)
      return(paste0("fixed_catchatage_indicator: ", check$message))
    ## obs_fixed_catchatage
    check = check_dim(data$obs_fixed_catchatage, c(2*n_ages, n_regions, n_years))
    if(!check$result)
      return(paste0("obs_fixed_catchatage: ", check$message))
    if(any(is.na(data$obs_fixed_catchatage)))
      return(paste0("NA's found in obs_fixed_catchatage"))

    ## trwl_catchatlgth_indicator
    check = check_dim(data$trwl_catchatlgth_indicator, c(n_regions, n_years))
    if(!check$result)
      return(paste0("trwl_catchatlgth_indicator: ", check$message))
    ## obs_trwl_catchatlgth
    check = check_dim(data$obs_trwl_catchatlgth, c(2*n_length_bins, n_regions, n_years))
    if(!check$result)
      return(paste0("obs_trwl_catchatlgth: ", check$message))
    if(any(is.na(data$obs_trwl_catchatlgth)))
      return(paste0("NA's found in obs_trwl_catchatlgth"))

    ## fixed_catchatlgth_indicator
    check = check_dim(data$fixed_catchatlgth_indicator, c(n_regions, n_years))
    if(!check$result)
      return(paste0("fixed_catchatlgth_indicator: ", check$message))
    ## obs_fixed_catchatlgth
    check = check_dim(data$obs_fixed_catchatlgth, c(2*n_length_bins, n_regions, n_years))
    if(!check$result)
      return(paste0("obs_fixed_catchatlgth: ", check$message))
    if(any(is.na(data$obs_fixed_catchatlgth)))
      return(paste0("NA's found in obs_fixed_catchatlgth"))

    ## srv_catchatage_indicator
    check = check_dim(data$srv_catchatage_indicator, c(n_regions, n_years, n_surveys))
    if(!check$result)
      return(paste0("srv_catchatage_indicator: ", check$message))
    ## obs_srv_catchatage
    check = check_dim(data$obs_srv_catchatage, c(2*n_ages, n_regions, n_years, n_surveys))
    if(!check$result)
      return(paste0("obs_srv_catchatage: ", check$message))
    if(any(is.na(data$obs_srv_catchatage)))
      return(paste0("NA's found in obs_catchatage"))

    ## srv_bio_indicator
    check = check_dim(data$srv_bio_indicator, c(n_regions, n_years, n_surveys))
    if(!check$result)
      return(paste0("srv_bio_indicator: ", check$message))
    ## obs_srv_bio
    check = check_dim(data$obs_srv_bio, c(n_regions, n_years, n_surveys))
    if(!check$result)
      return(paste0("obs_srv_bio: ", check$message))

    if(any(is.na(data$obs_srv_bio)))
      return(paste0("NA's found in obs_srv_bio"))
    ## obs_srv_se
    check = check_dim(data$obs_srv_se, c(n_regions, n_years, n_surveys))
    if(!check$result)
      return(paste0("obs_srv_se: ", check$message))

    if(any(is.na(data$obs_srv_se)))
      return(paste0("NA's found in obs_srv_se"))
    ## tag recovery observations
    ## tag_recovery_indicator
    check = check_length(data$tag_recovery_indicator_by_year, n_years)
    if(!check$result)
      return(paste0("tag_recovery_indicator_by_year: ", check$message))

    n_tag_recoveries = sum(data$tag_recovery_indicator_by_year)
    if(n_tag_recoveries > 0) {
      if(data$model == "TagIntegratedAgeBasedMovement") {
        if(data$age_based_movement == 1) {
          if (data$tag_likelihood != 2) {
            return("error in 'tag_likelihood'. When age_based_movement == 1, you should only use the multinomial likelihood. See https://craig44.github.io/SableFishResearch/agebasedmovement.html for a simulation on why")
          }
        }
        ## obs_tag_recovery
        check = check_dim(data$obs_tag_recovery, c(n_ages, n_regions * (data$n_years_to_retain_tagged_cohorts_for + 1), n_regions, n_tag_recoveries))
        if(!check$result)
          return(paste0("obs_tag_recovery: ", check$message))

      } else {
        if(data$tag_likelihood %in% c(0,1)) {
          ## tag_recovery_indicator
          check = check_dim(data$tag_recovery_indicator, c(n_regions * (data$n_years_to_retain_tagged_cohorts_for + 1), n_regions, n_tag_recoveries))
          if(!check$result)
            return(paste0("tag_recovery_indicator: ", check$message))
          ## obs_tag_recovery
          check = check_dim(data$obs_tag_recovery, c(n_regions * (data$n_years_to_retain_tagged_cohorts_for + 1), n_regions, n_tag_recoveries))
          if(!check$result)
            return(paste0("obs_tag_recovery: ", check$message))
        } else if(data$tag_likelihood %in% c(2)) {
          ## tag_recovery_indicator
          check = check_dim(data$tag_recovery_indicator, c(n_years, n_regions))
          if(!check$result)
            return(paste0("tag_recovery_indicator: ", check$message))
          ## obs_tag_recovery
          check = check_dim(data$obs_tag_recovery, c(n_regions * data$n_years_to_retain_tagged_cohorts_for + 1, n_regions, n_years))
          if(!check$result)
            return(paste0("obs_tag_recovery: ", check$message))
        }
      }
      if(any(is.na(data$obs_tag_recovery)))
        return("Found NA's in obs_tag_recovery")

    }
    ## check projection inputs. TODO: should only really check these if do_projection == 1
    if(!data$future_recruitment_type %in% c(0,1,2))
      return("Unknown value for future_recruitment_type, expected either 0, 1, or 2")
    check = check_length(data$year_ndx_for_empirical_resampling, 2)
    if(!check$result)
      return(paste0("year_ndx_for_empirical_resampling: ", check$message))

    ## parameters
    # ln_mean_rec
    check = check_length(parameters$ln_mean_rec, n_regions)
    if(!check$result)
      return(paste0("ln_mean_rec: ", check$message))
    # trans_rec_dev
    if(data$rec_devs_sum_to_zero == 0) {
      if(data$global_rec_devs == 1){
        check = check_dim(parameters$trans_rec_dev, c(1,n_years))
        if(!check$result)
          return(paste0("trans_rec_dev: ", check$message))

      } else if(data$global_rec_devs == 0) {
        check = check_dim(parameters$trans_rec_dev, c(n_regions,n_years))
        if(!check$result)
          return(paste0("trans_rec_dev: ", check$message))
      } else {
        return("Unknown input value for data$global_rec_devs")
      }
    } else {
      if(data$global_rec_devs == 1){
        check = check_dim(parameters$trans_rec_dev, c(1,n_years - 1))
        if(!check$result)
          return(paste0("trans_rec_dev: ", check$message))

      } else if(data$global_rec_devs == 0) {
        check = check_dim(parameters$trans_rec_dev, c(n_regions,n_years - 1))
        if(!check$result)
          return(paste0("trans_rec_dev: ", check$message))
      } else {
        return("Unknown input value for data$global_rec_devs")
      }
    }

    # fixed sel pars
    n_fixed_sel_time_blocks = length(unique(data$fixed_sel_by_year_indicator))
    if(!any(data$fixed_sel_by_year_indicator == 0))
      return("Could not find a 0 index in fixed_sel_by_year_indicator, this is likely an error")
    max_sel_pars = get_max_sel_pars(data$fixed_sel_type)
    check = check_dim(parameters$ln_fixed_sel_pars, c(n_fixed_sel_time_blocks, max_sel_pars, 2))
    if(!check$result)
      return(paste0("ln_fixed_sel_pars: ", check$message))

    # trawl sel pars
    n_fixed_sel_time_blocks = length(unique(data$trwl_sel_by_year_indicator))
    if(!any(data$trwl_sel_by_year_indicator == 0))
      return("Could not find a 0 index in trwl_sel_by_year_indicator, this is likely an error")
    max_sel_pars = get_max_sel_pars(data$trwl_sel_type)
    check = check_dim(parameters$ln_trwl_sel_pars, c(n_fixed_sel_time_blocks, max_sel_pars, 2))
    if(!check$result)
      return(paste0("ln_trwl_sel_pars: ", check$message))

    # Survey sel pars
    if(!any(data$srv_sel_by_year_indicator == 0))
      return("Could not find a 0 index in srv_sel_by_year_indicator, this is likely an error")
    max_sel_pars = max(apply(data$srv_sel_type, MARGIN = 2, function(x){get_max_sel_pars(x)}))
    check = check_dim(parameters$ln_srv_sel_pars, c(n_survey_sel_time_blocks, max_sel_pars, 2, n_surveys))
    if(!check$result)
      return(paste0("ln_srv_sel_pars: ", check$message))


    # srv_q_by_year_indicator
    check = check_dim(data$srv_q_by_year_indicator, c(n_years, n_surveys) )
    if(!check$result)
      return(paste0("srv_q_by_year_indicator: ", check$message))

    for(srv_ndx in 1:n_surveys) {
      if(any(data$q_is_nuisance[srv_ndx] == 0)) {
        n_fixed_survey_q_time_blocks = length(unique(data$srv_q_by_year_indicator[,srv_ndx]))
        if(!any(data$srv_q_by_year_indicator[srv_ndx] == 0))
          return("Could not find a 0 index in srv_q_by_year_indicator, this is likely an error")
      } else if(data$q_is_nuisance[srv_ndx] == 1) {
        if(!all(data$srv_q_by_year_indicator[,srv_ndx] == 0))
          return("When data$q_is_nuisance = 1, then srv_q_by_year_indicator must be all zeros")
      }
    }

    check = check_dim(parameters$trans_srv_q, c(n_regions, max(apply(data$srv_q_by_year_indicator, MARGIN = 2, function(x){length(unique(x))})), n_surveys))
    if(!check$result)
      return(paste0("trans_srv_q: ", check$message))

    # movement parameters transformed_movement_pars
    if(data$model == "TagIntegratedAgeBasedMovement") {
      if(n_regions > 1) {
        check = check_dim(parameters$transformed_movement_pars_old, c(n_regions - 1, n_regions))
        if(!check$result)
          return(paste0("transformed_movement_pars_old: ", check$message))
        check = check_dim(parameters$transformed_movement_pars_young, c(n_regions - 1, n_regions))
        if(!check$result)
          return(paste0("transformed_movement_pars_young: ", check$message))
      } else {
        check = check_dim(parameters$transformed_movement_pars_old, c(1, 1))
        if(!check$result)
          return(paste0("transformed_movement_pars_old: ", check$message))
        check = check_dim(parameters$transformed_movement_pars_young, c(1, 1))
        if(!check$result)
          return(paste0("transformed_movement_pars_young: ", check$message))
      }
    } else {
      if(n_regions > 1) {
        check = check_dim(parameters$transformed_movement_pars, c(n_regions - 1, n_regions))
        if(!check$result)
          return(paste0("transformed_movement_pars: ", check$message))
      } else {
        if(n_regions > 1) {
          check = check_dim(parameters$transformed_movement_pars, c(1, 1))
          if(!check$result)
            return(paste0("transformed_movement_pars: ", check$message))
        }
      }
    }
    # ln_fixed_F_devs
    check = check_dim(parameters$ln_fixed_F_devs, c(n_regions, n_years))
    if(!check$result)
      return(paste0("ln_fixed_F_devs: ", check$message))

    # ln_trwl_F_devs
    check = check_dim(parameters$ln_trwl_F_devs, c(n_regions, n_years))
    if(!check$result)
      return(paste0("ln_trwl_F_devs: ", check$message))

    # logistic_tag_reporting_rate
    check = check_dim(parameters$logistic_tag_reporting_rate, c(n_regions, max(n_tag_recoveries, 1)))
    if(!check$result)
      return(paste0("logistic_tag_reporting_rate: ", check$message))
    # ln_init_rec_dev
    if(data$n_init_rec_devs > 0) {
      check = check_length(parameters$ln_init_rec_dev, data$n_init_rec_devs)
      if(!check$result)
        return(paste0("ln_init_rec_dev: ", check$message))
    } else if(data$n_init_rec_devs == 0) {
      check = check_length(parameters$ln_init_rec_dev, 1)
      if(!check$result)
        return(paste0("ln_init_rec_dev: ", check$message))
    }
    ## logistic_prop_recruit_male
    check = check_length(parameters$logistic_prop_recruit_male, n_years)
    if(!check$result)
      return(paste0("logistic_prop_recruit_male: ", check$message))

    ## trans_fixed_catchatage_error
    if(data$fixed_catchatage_comp_likelihood == 1) {
      ## dirichlet multinomial
      check = check_length(parameters$trans_fixed_catchatage_error, 1)
      if(!check$result)
        return(paste0("trans_fixed_catchatage_error: ", check$message))
    }
    ## trans_fixed_catchatage_error
    if(data$fixed_catchatage_comp_likelihood == 1) {
      ## dirichlet multinomial
      check = check_length(parameters$trans_fixed_catchatage_error, 1)
      if(!check$result)
        return(paste0("trans_fixed_catchatage_error: ", check$message))
    }
    ## trans_fixed_catchatage_error
    if(data$fixed_catchatage_comp_likelihood == 1) {
      ## dirichlet multinomial
      check = check_length(parameters$trans_fixed_catchatage_error, 1)
      if(!check$result)
        return(paste0("trans_fixed_catchatage_error: ", check$message))
    }
    ## trans_fixed_catchatlgth_error
    if(data$fixed_catchatlgth_comp_likelihood == 1) {
      ## dirichlet multinomial
      check = check_length(parameters$trans_fixed_catchatlgth_error, 1)
      if(!check$result)
        return(paste0("trans_fixed_catchatlgth_error: ", check$message))
    }
    ## trans_trwl_catchatlgth_error
    if(data$trwl_catchatlgth_comp_likelihood == 1) {
      ## dirichlet multinomial
      check = check_length(parameters$trans_trwl_catchatlgth_error, 1)
      if(!check$result)
        return(paste0("trans_trwl_catchatlgth_error: ", check$message))
    }
    ## trans_trwl_catchatlgth_error
    for(srv_ndx in 1:n_surveys) {
      if(data$srv_catchatage_comp_likelihood[srv_ndx] == 1) {
        ## dirichlet multinomial
        check = check_length(parameters$trans_srv_catchatage_error[srv_ndx], 1)
        if(!check$result)
          return(paste0("trans_srv_catchatage_error: ", check$message))
      }
    }
  }

  print("Success!! Hopefully the model wont crash (no promises though)")
  return(TRUE)
}

