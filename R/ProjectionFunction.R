#' run_this_F
#' This function is run by the `find_fref` method
#' @param trial_F F value to F that is being trialed
#' @param target_percent_Bzero global percent B0 that we are trying to solve for
#' @param mle_obj an optimised obj that has been created by `TMB::MakeADFun`
#' @param fixed_F_allocation_by_region vector of proportions to allocated to each region for the fixed gear fishery
#' @param trwl_F_allocation_by_region vector of proportions to allocated to each region for the trawl gear fishery
#' @details the `trial_F` is allocated among both `fixed_F_allocation_by_region` and `trwl_F_allocation_by_region` which must sum = 0
#' @return sum of squares difference between the the trial F's Bzero and the value we are trying to solve for target_percent_Bzero
#' @export
run_this_F = function(trial_F, target_percent_Bzero, mle_obj, fixed_F_allocation_by_region, trwl_F_allocation_by_region) {
  ## set future F's
  mle_obj$env$data$future_fishing_inputs_fixed = matrix(trial_F * fixed_F_allocation_by_region, nrow = mle_obj$env$data$n_regions, ncol = mle_obj$env$data$n_projections_years, byrow = F)
  mle_obj$env$data$future_fishing_inputs_trwl = matrix(trial_F * trwl_F_allocation_by_region, nrow = mle_obj$env$data$n_regions, ncol = mle_obj$env$data$n_projections_years, byrow = F)
  ## run model
  proj_rep = mle_obj$report(get_tmb_fixed_effects(mle_obj))
  ## get SSBS
  this_ssb = proj_rep$SSB_yr
  ## get global depletion
  this_depletion = rowSums(this_ssb) / sum(proj_rep$Bzero)
  ## return sum of squares
  return((target_percent_Bzero - (this_depletion[length(this_depletion)] * 100))^2)
}


#' find_fref
#' deterministically search for an F that achieves some future Bzero value
#' @param mle_obj an optimised obj that has been created by `TMB::MakeADFun`
#' @param percent_Bzero This is the target biomass that will trigger the Fref value reported.
#' @details This algorithm will search over a range of F values to find one that closely results in percent B0 specified by percent_Bzero in the terminal projection year. .It is best practice to include future uncertainty when searching for Fmsy, it comes with considerable computational cost.
#' @return list with an F that has solved for a global F to achieve some percent_Bzero
#' @export

find_fref <- function(mle_obj, percent_Bzero = 40) {
  if(percent_Bzero < 0 | percent_Bzero > 100)
    stop("percent_Bzero is a percentage and needs to between 0 and 100")
  tmp_data = mle_obj$env$data
  mle_rep = mle_obj$report(get_tmb_fixed_effects(mle_obj))
  if(tmp_data$do_projection == 0) {
    ## set the do_projection = 1
    mle_obj$env$data$do_projection = 1
  }
  if(tmp_data$n_projections_years <= 0)
    stop("mle_obj: is not set up for projections. Found data$n_projections_years <= 0")

  ## set future fishing to F mode
  mle_obj$env$data$future_fishing_type = 0
  ## set future recruitment to mean recruitment no stochasticity
  mle_obj$env$data$future_recruitment_type = 2

  ## find F ratio between fisheries and regions
  terminal_F_by_fishery_and_region = c(mle_rep$annual_F_fixed[,length(tmp_data$years)], mle_rep$annual_F_trwl[,length(tmp_data$years)])

  fixed_F_allocation_by_region = mle_rep$annual_F_fixed[,length(tmp_data$years)] / sum(terminal_F_by_fishery_and_region)
  trwl_F_allocation_by_region = mle_rep$annual_F_trwl[,length(tmp_data$years)] / sum(terminal_F_by_fishery_and_region)
  if(abs(sum(c(fixed_F_allocation_by_region, trwl_F_allocation_by_region)) - 1.0) > 0.0001)
    stop("code error allocation should be close to 1.0. This was not the case")
  # Optimize for F that reaches percent Bzero
  optim_F = nlminb(start = 0.1, objective = run_this_F, target_percent_Bzero = percent_Bzero, mle_obj = mle_obj, fixed_F_allocation_by_region = fixed_F_allocation_by_region, trwl_F_allocation_by_region = trwl_F_allocation_by_region)
  # Run with reference F
  ## set future F's
  mle_obj$env$data$future_fishing_inputs_fixed = matrix(optim_F$par * fixed_F_allocation_by_region, nrow = mle_obj$env$data$n_regions, ncol = mle_obj$env$data$n_projections_years, byrow = F)
  mle_obj$env$data$future_fishing_inputs_trwl = matrix(optim_F$par * trwl_F_allocation_by_region, nrow = mle_obj$env$data$n_regions, ncol = mle_obj$env$data$n_projections_years, byrow = F)
  ## run model
  proj_rep = mle_obj$report(get_tmb_fixed_effects(mle_obj))
  ## get SSBS
  this_ssb = proj_rep$SSB_yr
  ## get global depletion
  this_depletion = rowSums(this_ssb) / sum(proj_rep$Bzero)
  ##
  return(list(proj_data = mle_obj$env$data, F_ref = optim_F$par, global_depletion = this_depletion))
}
