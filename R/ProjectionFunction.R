
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
  warning("this function is deprecated. Please use 'find_regional_Fspr'")

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
#' @param trial_Fs a vector of F values to be trialed should have length n_regions.
#' @param fixed_F_proportion for each region (trawl is is 1 - fixed_F_proportion)
#' @param input_pars vector of input parameters that are passed to proj_obj$report()
#' @param target_percent_Bzero regional percent B0 that we are trying to solve for in all regions
#' @param proj_obj an obj that has been created by `TMB::MakeADFun`
#' @param trace boolean print information as you are trialling Fs to help debug algorithm
#' @return sum of squares difference between the Bzero resulting from the trial F's and the value we are trying to solve for target_percent_Bzero
#' @export
run_multi_F = function(trial_Fs, fixed_F_proportion, input_pars, target_percent_Bzero, proj_obj, trace = F) {
  if(trace)
    cat("trialling Fs ", trial_Fs, "\n")
  ## set future F's
  n_regions =  proj_obj$env$data$n_regions
  proj_obj$env$data$future_fishing_inputs_fixed = matrix(trial_Fs * fixed_F_proportion, nrow = n_regions, ncol = proj_obj$env$data$n_projections_years, byrow = F)
  proj_obj$env$data$future_fishing_inputs_trwl = matrix(trial_Fs * (1 - fixed_F_proportion), nrow = n_regions, ncol = proj_obj$env$data$n_projections_years, byrow = F)
  ## run model
  proj_rep = proj_obj$report(input_pars)
  ## get SSBS
  this_ssb = proj_rep$SSB_yr
  ## get global depletion
  this_depletion = this_ssb[nrow(this_ssb),] / proj_rep$Bzero * 100
  ## return sum of squares
  return(sum((rep(target_percent_Bzero, n_regions) - this_depletion)^2))
}
#' find_regional_fref
#' deterministically search for a vector of Fs that achieves some regional percent b0 in all regions
#' @param proj_obj an obj that has been created by `TMB::MakeADFun`. It is assumed this obj object is setup for projection e.g., `do_projection = 1` etc
#' @param proj_pars vector of input parameters that are passed to proj_obj$report()
#' @param n_years_for_fleet_ratio an integer specifying how many terminal years to use to calculate F proportion by gear within each region
#' @param percent_Bzero This is the target biomass that will trigger the Fref value reported.
#' @param trace boolean print information as you are trialling Fs to help debug algorithm
#' @details This algorithm will search over a vector of  F values for each region to find values that result in percent B0 specified by percent_Bzero in the terminal projection year.
#' Consider using `setup_proj_data` when taking an estimation model and setting up a `proj_obj`. The ratio of F among fishing fleets within a region is specified by the input parameter `n_years_for_fleet_ratio`
#' @return list with estimated Fs that have been solved to achieve some percent_Bzero in all regions
#' @export

find_regional_fref <- function(proj_obj, proj_pars, n_years_for_fleet_ratio = 2, percent_Bzero = 40, trace = F) {
  warning("this function is deprecated. Please use 'find_regional_Fspr'")

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
  mean_fixed_F_by_region = mean_trwl_F_by_region = NULL
  if(n_regions == 1) {
    mean_fixed_F_by_region = mean(proj_rep$annual_F_fixed[,(length(tmp_data$years) - n_years_for_fleet_ratio + 1):length(tmp_data$years)])
    mean_trwl_F_by_region = mean(proj_rep$annual_F_trwl[,(length(tmp_data$years) - n_years_for_fleet_ratio + 1):length(tmp_data$years)])
  } else {
    if(n_years_for_fleet_ratio == 1) {
      mean_fixed_F_by_region = proj_rep$annual_F_fixed[,(length(tmp_data$years) - n_years_for_fleet_ratio + 1):length(tmp_data$years)]
      mean_trwl_F_by_region = proj_rep$annual_F_trwl[,(length(tmp_data$years) - n_years_for_fleet_ratio + 1):length(tmp_data$years)]

    } else {
      mean_fixed_F_by_region = apply(proj_rep$annual_F_fixed[,(length(tmp_data$years) - n_years_for_fleet_ratio + 1):length(tmp_data$years)], 1, mean)
      mean_trwl_F_by_region = apply(proj_rep$annual_F_trwl[,(length(tmp_data$years) - n_years_for_fleet_ratio + 1):length(tmp_data$years)], 1, mean)
    }
  }

  fixed_F_proportion = mean_fixed_F_by_region / (mean_fixed_F_by_region + mean_trwl_F_by_region)
  ## set future fishing to F mode
  proj_obj$env$data$future_fishing_type = 0
  ## set future recruitment to mean recruitment no stochasticity
  proj_obj$env$data$future_recruitment_type = 2

  # Optimize for F that reaches percent Bzero
  optim_F = nlminb(start = rep(mean(c(proj_rep$annual_F_fixed[,length(tmp_data$years)])), n_regions), objective = run_multi_F, input_pars = proj_pars, fixed_F_proportion = fixed_F_proportion,target_percent_Bzero = percent_Bzero, proj_obj = proj_obj, trace= trace, lower = rep(0.001, n_regions), control = list(x.tol = 1e-4))
  # Run with reference F
  test = run_multi_F(trial_Fs = optim_F$par, input_pars = proj_pars, fixed_F_proportion = fixed_F_proportion,target_percent_Bzero = percent_Bzero, proj_obj = proj_obj)
  ## set future F's
  proj_obj$env$data$future_fishing_inputs_fixed = matrix(optim_F$par * fixed_F_proportion, nrow = n_regions, ncol = proj_obj$env$data$n_projections_years, byrow = F)
  proj_obj$env$data$future_fishing_inputs_trwl = matrix(optim_F$par * (1 - fixed_F_proportion), nrow = n_regions, ncol = proj_obj$env$data$n_projections_years, byrow = F)
  ## run model
  proj_rep = proj_obj$report(proj_pars)
  return(list(proj_data = proj_obj$env$data, F_ref = optim_F$par, proj_rep = proj_rep))
}


#' trail_Fspr
#' This function is run by the `find_regional_Fspr` method to solve for an Fspr
#' @param log_trial_Fs a vector of log F values to be trialed should have length n_regions.
#' @param fixed_F_proportion for each region (trawl is is 1 - fixed_F_proportion)
#' @param n_future_years an integer, specifying how long to simulate ahead by for calculations. Must be long enough for the model to reach som convergence.
#' @param init_numbers_at_age matrix of initial numbers at age first value is assumed to be equal to 1. dim n_age x n_region
#' @param R0 vector of R0 parameters for each region
#' @param natural_mortality natural mortality vector for ages
#' @param maturity_weight_at_age vector of maturity and weight at age
#' @param trwl_sel vector of selectivities by age for the trawl fishery
#' @param fixed_sel vector of selectivities by age for the fixed gear fishery
#' @param movement_matrix matrix of annual movement
#' @param target_percent_Bzero regional percent B0 that we are trying to solve for in all regions
#' @param prop_Z_in_ssb scalar specifying the proportion of Z to account for when calculating SSB during the year
#' @param do_recruits_move 1 = yes 0 = no
#' @param verbose boolean print information as you are trialling Fs to help debug algorithm
#' @return sum of squares difference between the Bzero resulting from the trial F's and the value we are trying to solve for target_percent_Bzero
#' @export
trail_Fspr = function(log_trial_Fs, fixed_F_proportion, n_future_years = 60, init_numbers_at_age, R0, natural_mortality, maturity_weight_at_age, trwl_sel,fixed_sel, movement_matrix, target_percent_Bzero = 40, prop_Z_in_ssb = 0, verbose = F, do_recruits_move) {

  this_ssb = get_terminal_ssbs(trial_Fs = exp(log_trial_Fs), fixed_F_proportion = fixed_F_proportion, n_future_years = n_future_years, init_numbers_at_age = init_numbers_at_age, R0 = R0, natural_mortality = natural_mortality, maturity_weight_at_age = maturity_weight_at_age, trwl_sel = trwl_sel, fixed_sel = fixed_sel, movement_matrix = movement_matrix, target_percent_Bzero = target_percent_Bzero, prop_Z_in_ssb = prop_Z_in_ssb, do_recruits_move = do_recruits_move)

  SSE_terminal_Catch = sum((rep(target_percent_Bzero, length(log_trial_Fs)) - this_ssb$terminal_ssb/this_ssb$B0s * 100)^2)
  if(verbose)
    cat("trialling Fs ", exp(log_trial_Fs)," SSE = ", SSE_terminal_Catch, "\n")
  ## return sum of squares
  return(SSE_terminal_Catch)
}

#' get_terminal_ssbs
#' This function is run by the `find_regional_Fspr` method to solve for an Fspr. Given a an F for each region will return the percent B0 in each region at the terminal year
#' @param trial_Fs a vector of F values to be trialed should have length n_regions.
#' @param fixed_F_proportion for each region (trawl is is 1 - fixed_F_proportion)
#' @param n_future_years an integer, specifying how long to simulate ahead by for calculations. Must be long enough for the model to reach som convergence.
#' @param init_numbers_at_age matrix of initial numbers at age first value is assumed to be equal to 1. dim n_age x n_region
#' @param R0 vector of R0 parameters for each region
#' @param natural_mortality natural mortality vector for ages
#' @param maturity_weight_at_age vector of maturity and weight at age
#' @param trwl_sel vector of selectivities by age for the trawl fishery
#' @param fixed_sel vector of selectivities by age for the fixed gear fishery
#' @param movement_matrix matrix of annual movement
#' @param target_percent_Bzero regional percent B0 that we are trying to solve for in all regions
#' @param prop_Z_in_ssb scalar specifiying the proprotion of Z to account for when calculating SSB during the year
#' @param do_recruits_move 1 = yes 0 = no
#' @return sum of squares difference between the Bzero resulting from the trial F's and the value we are trying to solve for target_percent_Bzero
#' @export

get_terminal_ssbs <- function(trial_Fs, fixed_F_proportion, n_future_years = 60, init_numbers_at_age, R0, natural_mortality, maturity_weight_at_age, trwl_sel,fixed_sel, movement_matrix, target_percent_Bzero = 40, prop_Z_in_ssb = 0, do_recruits_move) {
  B0s = rowSums(sweep(init_numbers_at_age, MARGIN = 2, maturity_weight_at_age, FUN = "*"))
  n_ages = length(natural_mortality)
  n_regions = length(trial_Fs)
  Z_age = matrix(natural_mortality, nrow = length(trial_Fs), ncol = n_ages, byrow = T)
  fixed_F = matrix(fixed_sel, nrow = length(trial_Fs), ncol = n_ages, byrow = T)
  trwl_F = matrix(trwl_sel, nrow = length(trial_Fs), ncol = n_ages, byrow = T)
  fixed_F = sweep(fixed_F, MARGIN = 1, trial_Fs * fixed_F_proportion, FUN = "*")
  trwl_F = sweep(trwl_F, MARGIN = 1,trial_Fs * (1 - fixed_F_proportion), FUN = "*")
  Z_age = Z_age + trwl_F + fixed_F
  S_for_ssb = exp(-Z_age)^prop_Z_in_ssb
  update_N_age = N_age = init_numbers_at_age
  percent_ssbs = matrix(0,nrow = length(trial_Fs), ncol = n_future_years)
  current_ssb = NULL
  for(i in 1:n_future_years) {
    # interpolate SSB from start year
    current_ssb = rowSums(sweep(N_age * S_for_ssb, MARGIN = 2, maturity_weight_at_age, FUN = "*"))
    percent_ssbs[,i] = current_ssb / B0s * 100

    # recruitment
    if(do_recruits_move == 1) {
      update_N_age[,1] = R0 #* exp(-Z_age[,1])
    } else {
      update_N_age[,1] = N_age[,1] #* exp(-Z_age[,1])
    }
    # ageing and mortality
    update_N_age[,2:n_ages] = N_age[,1:(n_ages - 1)] * exp(-Z_age[,1:(n_ages - 1)] )
    # plus group
    update_N_age[,n_ages] = update_N_age[,n_ages] + N_age[,n_ages] * exp(-Z_age[,n_ages])
    ##
    # movement
    N_age = t(movement_matrix) %*% update_N_age
    ## if recruits don't move add them after movement occurs
    if(do_recruits_move == 0)
      N_age[,1] = R0
  }

  return(list(terminal_ssb = current_ssb, ssb_proj = percent_ssbs, B0s = B0s))
}

#' find_regional_Fspr
#' deterministically search for a vector of Fs (a value for each region) that achieves some regional percent b0
#' @param MLE_report a list that is output from obj$report() usually once an optimsation routine has been done.
#' @param n_years_for_fleet_ratio an integer specifying how many terminal years to average over to calculate F proportion by gear within each region
#' @param percent_Bzero This is the target biomass that will trigger the Fref value reported.
#' @param n_future_years an integer, specifying how long to simulate ahead by for calculations. Must be long enough for the model to reach som convergence.
#' @param verbose boolean print information as you are trialling Fs to help debug algorithm
#' @details This algorithm will search over a vector of  F values for each region to find values that result in percent B0 specified by percent_Bzero in the terminal projection year.
#' Consider using `setup_proj_data` when taking an estimation model and setting up a `proj_obj`. The ratio of F among fishing fleets within a region is specified by the input parameter `n_years_for_fleet_ratio`
#' @return list with estimated Fs that have been solved to achieve some percent_Bzero in all regions
#' @export

find_regional_Fspr <- function(data, MLE_report, n_years_for_fleet_ratio = 5, percent_Bzero = 40, n_future_years = 60, verbose = FALSE) {
  if(percent_Bzero < 0 | percent_Bzero > 100)
    stop("percent_Bzero is a percentage and needs to between 0 and 100")
  n_years = length(MLE_report$years)
  n_regions = data$n_regions
  n_ages = length(data$ages)
  ## get terminal maturity/weight at age for females
  weight_maturity_age = MLE_report$weight_maturity_prod_f[,n_years] ## note this will effect our B0 calculation/interpretation using the terminal growth

  ## get terminal selectivity by gear
  female_trwl_sel = MLE_report$sel_trwl_f[,ncol(MLE_report$sel_trwl_f)]
  female_fixed_sel = MLE_report$sel_fixed_f[,ncol(MLE_report$sel_fixed_f)]
  move_matrix = MLE_report$movement_matrix
  if(data$apply_fixed_movement == 1)
    move_matrix =  MLE_report$fixed_movement_matrix

  ## get terminal M matrix
  natural_M = data$M[,n_years]
  ## calcualte F ratio by gear type for each region
  mean_fixed_F_by_region = mean_trwl_F_by_region = NULL
  if(data$n_regions == 1) {
    mean_fixed_F_by_region = mean(MLE_report$annual_F_fixed[,(n_years - n_years_for_fleet_ratio + 1):n_years])
    mean_trwl_F_by_region = mean(MLE_report$annual_F_trwl[,(n_years - n_years_for_fleet_ratio + 1):n_years])
  } else {
    if(n_years_for_fleet_ratio == 1) {
      mean_fixed_F_by_region = MLE_report$annual_F_fixed[,(n_years - n_years_for_fleet_ratio + 1):n_years]
      mean_trwl_F_by_region = MLE_report$annual_F_trwl[,(n_years - n_years_for_fleet_ratio + 1):n_years]

    } else {
      mean_fixed_F_by_region = apply(MLE_report$annual_F_fixed[,(n_years - n_years_for_fleet_ratio + 1):n_years], 1, mean)
      mean_trwl_F_by_region = apply(MLE_report$annual_F_trwl[,(n_years - n_years_for_fleet_ratio + 1):n_years], 1, mean)
    }
  }
  fixed_F_proportion = mean_fixed_F_by_region / (mean_fixed_F_by_region + mean_trwl_F_by_region)

  ## calculate initial numbers at age
  init_numbers_at_age = calculate_initial_numbers_at_age(data$n_regions, n_ages, R0 = 0.5 * MLE_report$mean_rec, movement_matrix = move_matrix,
                                                         natural_mortality = natural_M)

  #MLE_report$equilibrium_natage_f[1,]
  #init_numbers_at_age[,1]
  #MLE_report$equilibrium_natage_f[29,]
  #init_numbers_at_age[,29]
  #MLE_report$equilibrium_natage_f[30,]
  #init_numbers_at_age[,30]
  #rowSums(init_numbers_at_age)

  start_Fs = apply(MLE_report$annual_F_fixed + MLE_report$annual_F_trwl, 1, mean)
  # Optimize for F that reaches percent Bzero
  optim_F = nlminb(start = log(start_Fs), objective = trail_Fspr, fixed_F_proportion = fixed_F_proportion, n_future_years = n_future_years, init_numbers_at_age = init_numbers_at_age, R0 = MLE_report$mean_rec * 0.5, natural_mortality = natural_M, maturity_weight_at_age = weight_maturity_age, fixed_sel = female_fixed_sel, trwl_sel = female_trwl_sel,
                   movement_matrix = move_matrix,
                   target_percent_Bzero = percent_Bzero, do_recruits_move = data$do_recruits_move,
                   verbose = verbose)


  terminal_ssb = get_terminal_ssbs(trial_Fs = exp(optim_F$par), fixed_F_proportion = fixed_F_proportion, n_future_years = n_future_years, init_numbers_at_age = init_numbers_at_age, R0 = MLE_report$mean_rec * 0.5, natural_mortality = natural_M, maturity_weight_at_age = weight_maturity_age, fixed_sel = female_fixed_sel, trwl_sel = female_trwl_sel,
                    movement_matrix = move_matrix, target_percent_Bzero = percent_Bzero, prop_Z_in_ssb = mean(data$spawning_time_proportion), do_recruits_move = data$do_recruits_move)

  #(rep(percent_Bzero, n_regions) - terminal_ssb$terminal_depletion)^2

  #terminal_ssb_start = get_terminal_ssbs(start_Fs, fixed_F_proportion = fixed_F_proportion, n_future_years = n_future_years, init_numbers_at_age = init_numbers_at_age, R0 = MLE_report$mean_rec * 0.5, natural_mortality = natural_M, maturity_weight_at_age = weight_maturity_age, fixed_sel = female_fixed_sel, trwl_sel = female_trwl_sel,
   #                                movement_matrix = move_matrix, target_percent_Bzero = percent_Bzero)
  #(rep(percent_Bzero, n_regions) - terminal_ssb_start$terminal_depletion)^2


  return(list(terminal_ssb = terminal_ssb, Fspr = exp(optim_F$par), fixed_gear_F_proportion = fixed_F_proportion))
}



