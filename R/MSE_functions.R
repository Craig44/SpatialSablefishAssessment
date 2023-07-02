#' simulate_future_data
#' @param data TMB data list which is up for current year
#' @param parameters TMB parameter which is up for current year
#' @param n_future_years number of years to run into the future
#' @param future_ll_Fs Annual fishing mortalities for the Longline fishery
#' @param future_trwl_Fs Annual fishing mortalities for the Trawl fishery
#' @param future_rec_devs future recruitment devs these are in logspace. If NULL we simulate from the prior
#' @details
#' Only used when `data$model == "Assessment"`
#' Any year based input quantity will be extended into the future using the last years values i.e., mean weight, M, age-length transition etc.
#'
#' @return list with the following elements sim_data, data, sim_parameters
#' @export
#' @examples
#' \dontrun{
#' # this loads in mock data and parameters
#' load(system.file("testdata", "MockAssessmentModel.RData",package="SpatialSablefishAssessment"))
#' n_future_years = 3
#' future_ll_Fs = rlnorm(n_future_years, log(0.3, 0.1))
#' future_trwl_Fs = rlnorm(n_future_years, log(0.3, 0.1))
#' future_rec_devs = rnorm(n_future_years, -0.5 * data$sigma_R^2, data$sigma_R)
#' # run function get simulated data
#' new_sim_data = simulate_future_data(data, parameters, n_future_years, future_ll_Fs, future_trwl_Fs,future_rec_devs)
#' }
simulate_future_data <-function(data, parameters, n_future_years, future_ll_Fs, future_trwl_Fs, future_rec_devs = NULL) {
  if(data$model != "Assessment")
    return("This function only works with 'data$model == Assessment'")
  ## do some sensible checks
  if(length(future_ll_Fs) != n_future_years)
    stop(paste0("'future_ll_Fs' needs to have length = ", n_future_years, " ('n_future_years'). Its current length is ", length(future_ll_Fs)))
  if(length(future_trwl_Fs) != n_future_years)
    stop(paste0("'future_trwl_Fs' needs to have length = ", n_future_years, " ('n_future_years'). Its current length is ", length(future_trwl_Fs)))

  ## simulate future recruitment
  if(is.null(future_rec_devs))
    future_rec_devs = rnorm(n_future_years, -0.5 * data$sigma_R^2, data$sigma_R)

  ## start with the parameters
  future_parameters = parameters
  future_parameters$ln_rec_dev = c(future_parameters$ln_rec_dev, future_rec_devs)
  # calculate F-dev future_parameters
  ln_ll_F_dev = log(future_trwl_Fs) - future_parameters$ln_ll_F_avg
  ln_trwl_F_dev = log(future_trwl_Fs) - future_parameters$ln_trwl_F_avg
  future_parameters$ln_ll_F_devs = c(future_parameters$ln_ll_F_devs, ln_ll_F_dev)
  future_parameters$ln_trwl_F_devs = c(future_parameters$ln_trwl_F_devs, ln_trwl_F_dev)
  # all other future_parameters are assumed to be the same in future years
  ## Now manipulate data objects
  future_data = data

  future_data$years = as.numeric(min(future_data$years):(max(future_data$years) + n_future_years))
  n_years = length(future_data$years)
  n_ages = length(data$ages)
  n_lens = length(data$length_bins)
  n_projyears = n_years +  future_data$n_projections_years
  # extend M
  future_data$M = cbind(future_data$M, matrix(future_data$M[,ncol(future_data$M)], nrow = n_ages, ncol = n_future_years, byrow = F))
  # extend maturity
  future_data$maturity = cbind(future_data$maturity, matrix(future_data$maturity[,ncol(future_data$maturity)], nrow = n_ages, ncol = n_future_years, byrow = F))
  # extend mean weight at age
  future_data$male_mean_weight_by_age = cbind(future_data$male_mean_weight_by_age, matrix(future_data$male_mean_weight_by_age[,ncol(future_data$male_mean_weight_by_age)], nrow = n_ages, ncol = n_future_years, byrow = F))
  future_data$female_mean_weight_by_age = cbind(future_data$female_mean_weight_by_age, matrix(future_data$female_mean_weight_by_age[,ncol(future_data$female_mean_weight_by_age)], nrow = n_ages, ncol = n_future_years, byrow = F))
  # proportion_male
  future_data$proportion_male = c(future_data$proportion_male, rep(mean(future_data$proportion_male), n_future_years))
  future_data$proportion_male2 = c(future_data$proportion_male2, rep(mean(future_data$proportion_male2), n_future_years))
  # spawning time
  future_data$spawning_time_proportion = c(future_data$spawning_time_proportion, rep(mean(future_data$spawning_time_proportion), n_future_years))
  # age-length-transition matrix
  AL_dim = dim(future_data$male_age_length_transition)
  tmp_male_AL = array(future_data$male_age_length_transition[,,AL_dim[3]], dim = c(AL_dim[1],AL_dim[2], n_future_years))
  tmp_female_AL = array(future_data$female_age_length_transition[,,AL_dim[3]], dim = c(AL_dim[1],AL_dim[2], n_future_years))
  future_data$male_age_length_transition = abind::abind(future_data$male_age_length_transition, tmp_male_AL)
  future_data$female_age_length_transition = abind::abind(future_data$female_age_length_transition, tmp_female_AL)
  # Catch - fill with empty values these will be simualted over based on input F's
  future_data$ll_fishery_catch = c(future_data$ll_fishery_catch, rep(1, n_future_years))
  future_data$trwl_fishery_catch = c(future_data$trwl_fishery_catch, rep(1, n_future_years))
  # selectivity and Q indicator variables
  future_data$ll_sel_by_year_indicator = c(future_data$ll_sel_by_year_indicator, rep(future_data$ll_sel_by_year_indicator[length(future_data$ll_sel_by_year_indicator)], n_future_years))
  future_data$trwl_sel_by_year_indicator = c(future_data$trwl_sel_by_year_indicator, rep(future_data$trwl_sel_by_year_indicator[length(future_data$trwl_sel_by_year_indicator)], n_future_years))
  future_data$srv_dom_ll_sel_by_year_indicator = c(future_data$srv_dom_ll_sel_by_year_indicator, rep(future_data$srv_dom_ll_sel_by_year_indicator[length(future_data$srv_dom_ll_sel_by_year_indicator)], n_future_years))
  future_data$srv_dom_ll_q_by_year_indicator = c(future_data$srv_dom_ll_q_by_year_indicator, rep(future_data$srv_dom_ll_q_by_year_indicator[length(future_data$srv_dom_ll_q_by_year_indicator)], n_future_years))
  future_data$srv_jap_ll_sel_by_year_indicator = c(future_data$srv_jap_ll_sel_by_year_indicator, rep(future_data$srv_jap_ll_sel_by_year_indicator[length(future_data$srv_jap_ll_sel_by_year_indicator)], n_future_years))
  future_data$srv_jap_ll_q_by_year_indicator = c(future_data$srv_jap_ll_q_by_year_indicator, rep(future_data$srv_jap_ll_q_by_year_indicator[length(future_data$srv_jap_ll_q_by_year_indicator)], n_future_years))
  future_data$srv_nmfs_trwl_sel_by_year_indicator = c(future_data$srv_nmfs_trwl_sel_by_year_indicator, rep(future_data$srv_nmfs_trwl_sel_by_year_indicator[length(future_data$srv_nmfs_trwl_sel_by_year_indicator)], n_future_years))
  future_data$srv_nmfs_trwl_q_by_year_indicator = c(future_data$srv_nmfs_trwl_q_by_year_indicator, rep(future_data$srv_nmfs_trwl_q_by_year_indicator[length(future_data$srv_nmfs_trwl_q_by_year_indicator)], n_future_years))
  future_data$srv_jap_fishery_ll_sel_by_year_indicator = c(future_data$srv_jap_fishery_ll_sel_by_year_indicator, rep(future_data$srv_jap_fishery_ll_sel_by_year_indicator[length(future_data$srv_jap_fishery_ll_sel_by_year_indicator)], n_future_years))
  future_data$srv_jap_fishery_ll_q_by_year_indicator = c(future_data$srv_jap_fishery_ll_q_by_year_indicator, rep(future_data$srv_jap_fishery_ll_q_by_year_indicator[length(future_data$srv_jap_fishery_ll_q_by_year_indicator)], n_future_years))
  future_data$ll_cpue_q_by_year_indicator = c(future_data$ll_cpue_q_by_year_indicator, rep(future_data$ll_cpue_q_by_year_indicator[length(future_data$ll_cpue_q_by_year_indicator)], n_future_years))

  # Observation indicator and input containers need changing
  # observations containers will be filled with the first value.
  # as the effective sample size from the last year.
  # these will be over written by simulated values so we just want to make sure
  # we specify the correct observation error.
  # LL age
  future_data$ll_catchatage_indicator = c(future_data$ll_catchatage_indicator, rep(1, n_future_years))
  N_eff = sum(future_data$obs_ll_catchatage[,ncol(future_data$obs_ll_catchatage)])
  future_data$obs_ll_catchatage = matrix(N_eff / n_ages, nrow = n_ages, ncol = sum(future_data$ll_catchatage_indicator))
  # LL Length male
  future_data$ll_catchatlgth_indicator = c(future_data$ll_catchatlgth_indicator, rep(1, n_future_years))
  N_eff = sum(future_data$obs_ll_catchatlgth_m[,ncol(future_data$obs_ll_catchatlgth_m)])
  future_data$obs_ll_catchatlgth_m = matrix(N_eff / n_lens, nrow = n_lens, ncol = sum(future_data$ll_catchatlgth_indicator))
  # LL Length female
  N_eff = sum(future_data$obs_ll_catchatlgth_f[,ncol(future_data$obs_ll_catchatlgth_f)])
  future_data$obs_ll_catchatlgth_f = matrix(N_eff / n_lens, nrow = n_lens, ncol = sum(future_data$ll_catchatlgth_indicator))
  # Trawl length
  future_data$trwl_catchatlgth_indicator = c(future_data$trwl_catchatlgth_indicator, rep(1, n_future_years))
  N_eff = sum(future_data$obs_trwl_catchatlgth_m[,ncol(future_data$obs_trwl_catchatlgth_m)])
  future_data$obs_trwl_catchatlgth_m = matrix(N_eff / n_lens, nrow = n_lens, ncol = sum(future_data$trwl_catchatlgth_indicator))
  # LL Length female
  N_eff = sum(future_data$obs_trwl_catchatlgth_f[,ncol(future_data$obs_trwl_catchatlgth_f)])
  future_data$obs_trwl_catchatlgth_f = matrix(N_eff / n_lens, nrow = n_lens, ncol = sum(future_data$trwl_catchatlgth_indicator))
  # Survey biomass for the Domestic Longline survey
  future_data$srv_dom_ll_bio_indicator = c(future_data$srv_dom_ll_bio_indicator, rep(1, n_future_years))
  future_data$obs_dom_ll_bio = c(future_data$obs_dom_ll_bio, rep(1.0, n_future_years))
  future_data$se_dom_ll_bio = c(future_data$se_dom_ll_bio, rep(future_data$se_dom_ll_bio[length(future_data$se_dom_ll_bio)], n_future_years))
  # Survey biomass for the Japanease Longline survey
  future_data$srv_jap_ll_bio_indicator = c(future_data$srv_jap_ll_bio_indicator, rep(1, n_future_years))
  future_data$obs_jap_ll_bio = c(future_data$obs_jap_ll_bio, rep(1.0, n_future_years))
  future_data$se_jap_ll_bio = c(future_data$se_jap_ll_bio, rep(future_data$se_jap_ll_bio[length(future_data$se_jap_ll_bio)], n_future_years))
  # Survey NMFS GOA trawl survey
  future_data$srv_nmfs_trwl_bio_indicator = c(future_data$srv_nmfs_trwl_bio_indicator, rep(1, n_future_years))
  future_data$obs_nmfs_trwl_bio = c(future_data$obs_nmfs_trwl_bio, rep(1.0, n_future_years))
  future_data$se_nmfs_trwl_bio = c(future_data$se_nmfs_trwl_bio, rep(future_data$se_nmfs_trwl_bio[length(future_data$se_nmfs_trwl_bio)], n_future_years))
  # Longline CPUE
  future_data$ll_cpue_indicator = c(future_data$ll_cpue_indicator, rep(1, n_future_years))
  future_data$obs_ll_cpue = c(future_data$obs_ll_cpue, rep(1.0, n_future_years))
  future_data$se_ll_cpue = c(future_data$se_ll_cpue, rep(future_data$se_ll_cpue[length(future_data$se_ll_cpue)], n_future_years))
  # Japanease Longline Fishery
  future_data$srv_jap_fishery_ll_bio_indicator = c(future_data$srv_jap_fishery_ll_bio_indicator, rep(1, n_future_years))
  future_data$obs_jap_fishery_ll_bio = c(future_data$obs_jap_fishery_ll_bio, rep(1.0, n_future_years))
  future_data$se_jap_fishery_ll_bio = c(future_data$se_jap_fishery_ll_bio, rep(future_data$se_jap_fishery_ll_bio[length(future_data$se_jap_fishery_ll_bio)], n_future_years))
  # Survey age for the Domestic Longline
  future_data$srv_dom_ll_age_indicator = c(future_data$srv_dom_ll_age_indicator, rep(1, n_future_years))
  N_eff = sum(future_data$obs_srv_dom_ll_age[,ncol(future_data$obs_srv_dom_ll_age)])
  future_data$obs_srv_dom_ll_age = matrix(N_eff / n_ages, nrow = n_ages, ncol = sum(future_data$srv_dom_ll_age_indicator))
  # Survey length Comp from the Domestic Longline survey
  future_data$srv_dom_ll_lgth_indicator = c(future_data$srv_dom_ll_lgth_indicator, rep(1, n_future_years))
  N_eff = sum(future_data$obs_srv_dom_ll_lgth_m[,ncol(future_data$obs_srv_dom_ll_lgth_m)])
  future_data$obs_srv_dom_ll_lgth_m = matrix(N_eff / n_lens, nrow = n_lens, ncol = sum(future_data$srv_dom_ll_lgth_indicator))

  N_eff = sum(future_data$obs_srv_dom_ll_lgth_f[,ncol(future_data$obs_srv_dom_ll_lgth_f)])
  future_data$obs_srv_dom_ll_lgth_f = matrix(N_eff / n_lens, nrow = n_lens, ncol = sum(future_data$srv_dom_ll_lgth_indicator))
  # Japanese LL age early survey
  future_data$srv_jap_ll_age_indicator = c(future_data$srv_jap_ll_age_indicator, rep(1, n_future_years))
  N_eff = sum(future_data$obs_srv_jap_ll_age[,ncol(future_data$obs_srv_jap_ll_age)])
  future_data$obs_srv_jap_ll_age = matrix(N_eff / n_ages, nrow = n_ages, ncol = sum(future_data$srv_jap_ll_age_indicator))
  # Japanese LL length Comp early survey
  future_data$srv_jap_ll_lgth_indicator = c(future_data$srv_jap_ll_lgth_indicator, rep(1, n_future_years))
  N_eff = sum(future_data$obs_srv_jap_ll_lgth_m[,ncol(future_data$obs_srv_jap_ll_lgth_m)])
  future_data$obs_srv_jap_ll_lgth_m = matrix(N_eff / n_lens, nrow = n_lens, ncol = sum(future_data$srv_jap_ll_lgth_indicator))
  N_eff = sum(future_data$obs_srv_jap_ll_lgth_f[,ncol(future_data$obs_srv_jap_ll_lgth_f)])
  future_data$obs_srv_jap_ll_lgth_f = matrix(N_eff / n_lens, nrow = n_lens, ncol = sum(future_data$srv_jap_ll_lgth_indicator))
  # Japanese LL Fishery length Comp early survey (Sex aggregated)
  future_data$srv_jap_fishery_ll_lgth_indicator = c(future_data$srv_jap_fishery_ll_lgth_indicator, rep(1, n_future_years))
  N_eff = sum(future_data$obs_srv_jap_fishery_ll_lgth[,ncol(future_data$obs_srv_jap_fishery_ll_lgth)])
  future_data$obs_srv_jap_fishery_ll_lgth = matrix(N_eff / n_lens, nrow = n_lens, ncol = sum(future_data$srv_jap_fishery_ll_lgth_indicator))
  # NMFS bottom trawl age frequency
  future_data$srv_nmfs_trwl_age_indicator = c(future_data$srv_nmfs_trwl_age_indicator, rep(1, n_future_years))
  N_eff = sum(future_data$obs_srv_nmfs_trwl_age[,ncol(future_data$obs_srv_nmfs_trwl_age)])
  future_data$obs_srv_nmfs_trwl_age = matrix(N_eff / n_ages, nrow = n_ages, ncol = sum(future_data$srv_nmfs_trwl_age_indicator))
  #  NMFS bottom trawl length Comp male
  future_data$srv_nmfs_trwl_lgth_indicator = c(future_data$srv_nmfs_trwl_lgth_indicator, rep(1, n_future_years))
  N_eff = sum(future_data$obs_srv_nmfs_trwl_lgth_m[,ncol(future_data$obs_srv_nmfs_trwl_lgth_m)])
  future_data$obs_srv_nmfs_trwl_lgth_m = matrix(N_eff / n_lens, nrow = n_lens, ncol = sum(future_data$srv_nmfs_trwl_lgth_indicator))
  # NMFS bottom trawl length Comp female
  N_eff = sum(future_data$obs_srv_nmfs_trwl_lgth_f[,ncol(future_data$obs_srv_nmfs_trwl_lgth_f)])
  future_data$obs_srv_nmfs_trwl_lgth_f = matrix(N_eff / n_lens, nrow = n_lens, ncol = sum(future_data$srv_nmfs_trwl_lgth_indicator))

  ## validate future_data and future_parameters
  result = validate_input_data_and_parameters(future_data, future_parameters)
  if(!result)
    stop("This function needs adapting. We failed the 'validate_input_data_and_future_parameters' check. If we progress this will likely cause an system crash.")
  ## Build OM
  OM <- TMB::MakeADFun(data = future_data,
                       parameters = future_parameters,
                       DLL = "SpatialSablefishAssessment_TMBExports", silent  = T, checkParameterOrder = T)

  sim_data = OM$simulate(complete = T)
  sim_data = convert_simdata_integers(sim_data, future_data)

  return(list(sim_data = sim_data, future_data = future_data, future_parameters = future_parameters))
}
