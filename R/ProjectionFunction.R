
#' setup_proj_data
#' @param mle_obj an optimised obj that has been created by `TMB::MakeADFun`
#' @param n_proj_years number of years you want to project out to
#' @param future_recruitment the type of future recruitment during projection phase
#' @details this function will extend all the biological objects and selectivity objects so that they have the correct dimensions.
#' This is done by repeating the value at n_years or dimension for all projection periods
#' this function will return a data list that you can further modify if you wish.
#' @export
#' @return data list that can be used by `TMB::MakeADFun` to make a projection model.
setup_proj_data <- function(mle_obj, n_proj_years = 100, future_recruitment = 2) {
  tmb_data = mle_obj$env$data
  n_regions = tmb_data$n_regions
  n_ages = length(tmb_data$ages)
  n_length_bins = length(tmb_data$length_bins) # the last length bin value is the minimum for a length plus group
  tmb_data$n_projections_years = n_proj_years
  tmb_data$do_projection = 1
  n_projyears = length(tmb_data$years) +  tmb_data$n_projections_years
  n_years = length(tmb_data$years)
  projyears = min(tmb_data$years):(max(tmb_data$years) + tmb_data$n_projections_years)

  n_yrs_to_append =  n_projyears - n_years
  ## now extend all necessary containers
  ## biological values
  tmb_data$M = extend_2darray(tmb_data$M, n_yrs_to_append)
  tmb_data$maturity = extend_2darray(tmb_data$maturity, n_yrs_to_append)
  tmb_data$female_mean_weight_by_age = extend_2darray(tmb_data$female_mean_weight_by_age, n_yrs_to_append)
  tmb_data$male_mean_weight_by_age = extend_2darray(tmb_data$male_mean_weight_by_age, n_yrs_to_append)
  tmb_data$spawning_time_proportion = extend_vec_last_val(tmb_data$spawning_time_proportion, n_yrs_to_append)


  ## age length matricies
  tmb_data$female_age_length_transition = extend_3darray_last_dim(tmb_data$female_age_length_transition, n_yrs_to_append)
  tmb_data$male_age_length_transition = extend_3darray_last_dim(tmb_data$male_age_length_transition, n_yrs_to_append)

  ## selectivity indicators
  tmb_data$fixed_sel_by_year_indicator = extend_vec_last_val(tmb_data$fixed_sel_by_year_indicator, n_yrs_to_append)
  tmb_data$trwl_sel_by_year_indicator = extend_vec_last_val(tmb_data$trwl_sel_by_year_indicator, n_yrs_to_append)
  tmb_data$srv_dom_ll_sel_by_year_indicator = extend_vec_last_val(tmb_data$srv_dom_ll_sel_by_year_indicator, n_yrs_to_append)

  ## Change future inputs
  tmb_data$future_recruitment_type = future_recruitment = 2
  tmb_data$future_fishing_inputs_trwl = matrix(0.1, nrow = tmb_data$n_regions, ncol = tmb_data$n_projections_years)
  tmb_data$future_fishing_inputs_fixed = matrix(0.1, nrow = tmb_data$n_regions, ncol = tmb_data$n_projections_years)
  tmb_data$future_fishing_type = 0
  return(tmb_data)
}




#' run_this_F
#' This function is run by the `find_apportioned_fref` method
#' @param trial_F F value to F that is being trialed
#' @param input_pars vector of input parameters that are passed to proj_obj$report()
#' @param target_percent_Bzero global percent B0 that we are trying to solve for
#' @param proj_obj an obj that has been created by `TMB::MakeADFun`
#' @param fixed_F_allocation_by_region vector of proportions to allocated to each region for the fixed gear fishery
#' @param trwl_F_allocation_by_region vector of proportions to allocated to each region for the trawl gear fishery
#' @details the `trial_F` is allocated among both `fixed_F_allocation_by_region` and `trwl_F_allocation_by_region` which must sum = 0
#' @return sum of squares difference between the the trial F's Bzero and the value we are trying to solve for target_percent_Bzero
#' @export
run_this_F = function(trial_F, input_pars, target_percent_Bzero, proj_obj, fixed_F_allocation_by_region, trwl_F_allocation_by_region) {
  ## set future F's
  proj_obj$env$data$future_fishing_inputs_fixed = matrix(trial_F * fixed_F_allocation_by_region, nrow = proj_obj$env$data$n_regions, ncol = proj_obj$env$data$n_projections_years, byrow = F)
  proj_obj$env$data$future_fishing_inputs_trwl = matrix(trial_F * trwl_F_allocation_by_region, nrow = proj_obj$env$data$n_regions, ncol = proj_obj$env$data$n_projections_years, byrow = F)
  ## run model
  proj_rep = proj_obj$report(input_pars)
  ## get SSBS
  this_ssb = proj_rep$SSB_yr
  ## get global depletion
  this_depletion = rowSums(this_ssb) / sum(proj_rep$Bzero)
  ## return sum of squares
  return((target_percent_Bzero - (this_depletion[length(this_depletion)] * 100))^2)
}


#' find_apportioned_fref
#' deterministically search for an F that achieves some future Bzero value
#' @param proj_obj an obj that has been created by `TMB::MakeADFun`. It is assumed this obj object is setup for projection e.g., `do_projection = 1` etc
#' @param proj_pars vector of input parameters that are passed to proj_obj$report()
#' @param percent_Bzero This is the target biomass that will trigger the Fref value reported.
#' @details This algorithm will search over a range of F values to find one that closely results in percent B0 specified by percent_Bzero in the terminal projection year.
#' Consider using `setup_proj_data` when taking an estimation model and setting up a `proj_obj`
#' @return list with an F that has solved for a global F to achieve some percent_Bzero
#' @export

find_apportioned_fref <- function(proj_obj, proj_pars, percent_Bzero = 40) {
  if(percent_Bzero < 0 | percent_Bzero > 100)
    stop("percent_Bzero is a percentage and needs to between 0 and 100")
  tmp_data = proj_obj$env$data
  proj_rep = proj_obj$report(proj_pars)
  if(tmp_data$do_projection == 0) {
    ## set the do_projection = 1
    proj_obj$env$data$do_projection = 1
  }
  if(tmp_data$n_projections_years <= 0)
    stop("proj_obj: is not set up for projections. Found data$n_projections_years <= 0")

  ## set future fishing to F mode
  proj_obj$env$data$future_fishing_type = 0
  ## set future recruitment to mean recruitment no stochasticity
  proj_obj$env$data$future_recruitment_type = 2

  ## find F ratio between fisheries and regions
  terminal_F_by_fishery_and_region = c(proj_rep$annual_F_fixed[,length(tmp_data$years)], proj_rep$annual_F_trwl[,length(tmp_data$years)])

  fixed_F_allocation_by_region = proj_rep$annual_F_fixed[,length(tmp_data$years)] / sum(terminal_F_by_fishery_and_region)
  trwl_F_allocation_by_region = proj_rep$annual_F_trwl[,length(tmp_data$years)] / sum(terminal_F_by_fishery_and_region)
  if(abs(sum(c(fixed_F_allocation_by_region, trwl_F_allocation_by_region)) - 1.0) > 0.0001)
    stop("code error allocation should be close to 1.0. This was not the case")
  # Optimize for F that reaches percent Bzero
  optim_F = nlminb(start = 0.1, objective = run_this_F, lower = 0.000001,input_pars = proj_pars, target_percent_Bzero = percent_Bzero, proj_obj = proj_obj, fixed_F_allocation_by_region = fixed_F_allocation_by_region, trwl_F_allocation_by_region = trwl_F_allocation_by_region)
  # Run with reference F
  ## set future F's
  proj_obj$env$data$future_fishing_inputs_fixed = matrix(optim_F$par * fixed_F_allocation_by_region, nrow = proj_obj$env$data$n_regions, ncol = proj_obj$env$data$n_projections_years, byrow = F)
  proj_obj$env$data$future_fishing_inputs_trwl = matrix(optim_F$par * trwl_F_allocation_by_region, nrow = proj_obj$env$data$n_regions, ncol = proj_obj$env$data$n_projections_years, byrow = F)
  ## run model
  proj_rep = proj_obj$report(proj_pars)
  ## get SSBS
  this_ssb = proj_rep$SSB_yr
  ## get global depletion
  this_depletion = rowSums(this_ssb) / sum(proj_rep$Bzero)
  ##
  return(list(proj_data = proj_obj$env$data, F_ref = optim_F$par, global_depletion = this_depletion, proj_rep = proj_rep))
}


#' run_multi_F
#' This function is run by the `find_multi_fref` method
#' @param trial_Fs a vector of F values to be trialed should have length n_regions * 2. The first n_regions are for the fixed gear fishery and the second set are for the trawl
#' @param input_pars vector of input parameters that are passed to proj_obj$report()
#' @param target_percent_Bzero regional percent B0 that we are trying to solve for in all regions
#' @param proj_obj an obj that has been created by `TMB::MakeADFun`
#' @return sum of squares difference between the Bzero resulting from the trial F's and the value we are trying to solve for target_percent_Bzero
#' @export
run_multi_F = function(trial_Fs, input_pars, target_percent_Bzero, proj_obj) {
  ## set future F's
  n_regions =  proj_obj$env$data$n_regions
  proj_obj$env$data$future_fishing_inputs_fixed = matrix(trial_Fs[1:n_regions], nrow = n_regions, ncol = proj_obj$env$data$n_projections_years, byrow = F)
  proj_obj$env$data$future_fishing_inputs_trwl = matrix(trial_Fs[(n_regions + 1):(2*n_regions)], nrow = n_regions, ncol = proj_obj$env$data$n_projections_years, byrow = F)
  ## run model
  proj_rep = proj_obj$report(input_pars)
  ## get SSBS
  this_ssb = proj_rep$SSB_yr
  ## get global depletion
  this_depletion = this_ssb[nrow(this_ssb)] /proj_rep$Bzero
  ## return sum of squares
  return(sum((rep(target_percent_Bzero, n_regions) - this_depletion * 100)^2))
}
#' find_regional_fref
#' deterministically search for a vector of Fs that achieves some regional percent b0 in all regions
#' @param proj_obj an obj that has been created by `TMB::MakeADFun`. It is assumed this obj object is setup for projection e.g., `do_projection = 1` etc
#' @param proj_pars vector of input parameters that are passed to proj_obj$report()
#' @param percent_Bzero This is the target biomass that will trigger the Fref value reported.
#' @details This algorithm will search over a vector of  F values to find values that result in percent B0 specified by percent_Bzero in the terminal projection year.
#' Consider using `setup_proj_data` when taking an estimation model and setting up a `proj_obj`
#' @return list with estimated Fs that have been solved to achieve some percent_Bzero in all regions
#' @export

find_regional_fref <- function(proj_obj, proj_pars, percent_Bzero = 40) {
  if(percent_Bzero < 0 | percent_Bzero > 100)
    stop("percent_Bzero is a percentage and needs to between 0 and 100")
  tmp_data = proj_obj$env$data
  n_regions = tmp_data$n_regions
  proj_rep = proj_obj$report(proj_pars)
  if(tmp_data$do_projection == 0) {
    ## set the do_projection = 1
    proj_obj$env$data$do_projection = 1
  }
  if(tmp_data$n_projections_years <= 0)
    stop("proj_obj: is not set up for projections. Found data$n_projections_years <= 0")

  ## set future fishing to F mode
  proj_obj$env$data$future_fishing_type = 0
  ## set future recruitment to mean recruitment no stochasticity
  proj_obj$env$data$future_recruitment_type = 2

  start_Fs = c(proj_rep$annual_F_fixed[,length(tmp_data$years)], proj_rep$annual_F_trwl[,length(tmp_data$years)])
  # Optimize for F that reaches percent Bzero
  optim_F = nlminb(start = start_Fs, objective = run_multi_F,input_pars = proj_pars, target_percent_Bzero = percent_Bzero, proj_obj = proj_obj, lower = rep(0.0000001, length(start_Fs)))
  # Run with reference F
  ## set future F's
  proj_obj$env$data$future_fishing_inputs_fixed = matrix(optim_F$par[1:n_regions], nrow = n_regions, ncol = proj_obj$env$data$n_projections_years, byrow = F)
  proj_obj$env$data$future_fishing_inputs_trwl = matrix(optim_F$par[(n_regions + 1):length(optim_F$par)], nrow = n_regions, ncol = proj_obj$env$data$n_projections_years, byrow = F)
  ## run model
  proj_rep = proj_obj$report(proj_pars)
  return(list(proj_data = proj_obj$env$data, F_ref = optim_F$par, proj_rep = proj_rep))
}
