/* @file TagIntegrated.hpp
 * Statistical, separable space, sex and age-structured population model for sablefish
 * Alaska Fisheries Science Center, October 2022
 * Written by C Marsh craig.marsh10@gmail.com
 * This follows on from the current assessment that was originally written by
 * D. Hanselman:dana.hanselman@noaa.gov
 * UPDATED (and commented)  by D. Goethel: daniel.goethel@noaa.gov   (10/15/20)
 * Tips
 * - TMB indicies start at 0 (i.e., like C++) where as ADMB starts at 1 (i.e., like R)
 * - parameter labels should start with the transformation that is assumed for example natural log of F for longline should follow ln_F_fixed
 *   and for logistic proportion something like logis_prop. This is to aid readability and keep syntax consistent
 *
 * - Key when reading object names
 *      - dom = domestic
 *      - fixed = fixed gear LL + Pot
 *      - trwl = Trawl
 *      - srv = survey
 *      - lgth = Length
 *
 *
 */

#ifndef TagIntegrated_hpp
#define TagIntegrated_hpp

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

//
template<class Type>
Type TagIntegrated(objective_function<Type>* obj) {
  using namespace density;

  // Input parameters
  // model dimensions
  DATA_VECTOR(ages);                                // assumes min(ages) >= 1, also assumes the last age is a plus group
  DATA_VECTOR(years);                               // annual years
  DATA_VECTOR(length_bins);                         // Length bins, the last length bin value is the minimum for a length plus group
  DATA_INTEGER(n_projections_years);                // number of years to project the model beyond max(years)
  DATA_INTEGER(do_projection);                      // Should we project the model to last_projection_year. 1 = yes, 0 = no
  DATA_INTEGER(n_regions);                          // number of regions in the model

  int n_years = years.size();
  int n_projyears = n_years + n_projections_years;
  int n_ages = ages.size();
  int n_length_bins = length_bins.size();


  // Biology parameters
  DATA_INTEGER(global_rec_devs);                    // Are there recruit devs parameters for each region (= 0), or do all regions have the same rec devs (=1)
  DATA_INTEGER(rec_devs_sum_to_zero);               // Should the recruit devs in each region sum to zero? yes = 1, no = 0. I yes then this the parameter trans_rec_dev has one less parameter
  DATA_VECTOR(Q_r_for_sum_to_zero);                 // A vector that has length (trans_rec_dev.dim(1) + 1) * 2. Only used if rec_devs_sum_to_zero = 1. Should have been created by the R function Q_sum_to_zero_QR

  DATA_INTEGER(n_init_rec_devs);                    // Number of initial recruitment devs parameters "init_ln_rec_dev" Note: should cannot be greater than n_ages - 2 (we don't apply it to first age or plus group)

  // this will effect the expected size of the parameter 'ln_rec_dev', if global_rec_devs = 1. then ln_rec_dev.size() = n_years + n_ages + 1 else n_regions * (n_years + n_ages + 1). with the first (n_years + n_ages + 1) corresponding to region 1 and so in block
  DATA_ARRAY(M);                                    // Natural Mortality: dim = n_ages x n_projyears
  DATA_ARRAY(maturity);                             // Proportion ages mature: dim = n_ages x n_projyears
  DATA_ARRAY(male_mean_weight_by_age);              // male_mean_weight_by_age (tonnes): dim = n_ages x n_projyears
  DATA_ARRAY(female_mean_weight_by_age);            // female_mean_weight_by_age (tonnes): dim = n_ages x n_projyears

  DATA_ARRAY(male_age_length_transition);           // Proportion at among length bins for each age for male: dim = n_ages x n_lengths x n_years
  DATA_ARRAY(female_age_length_transition);         // Proportion at among length bins for each age for female: dim = n_ages x n_lengths x n_years


  DATA_INTEGER(SrType);                             // Stock recruitment type 3=average, 2=Bholt, 1=Ricker
  DATA_VECTOR(spawning_time_proportion);            // proportion of time within a year that spawning occurs needed for each year length = n_projyears, bound between 0 and 1

  // the reason I added this fixed fixed movement switch is because we use the simplex to transform parameters for the estimated movement
  // matrix. This parameterisation will cause NaNs or Inf when there is a value = 1 and the rest zeros i.e., no movement. For this scenario you
  // you should use this input functionality.
  DATA_INTEGER(apply_fixed_movement);               // 0 means will use estimated movement matrix, 1 means will use input movement matrix.
  DATA_MATRIX(fixed_movement_matrix);               // n_regions x n_regions. only used if apply_fixed_movement = 1

  // Fishing stuff
  DATA_SCALAR(prop_F_hist);                         // Proportion of fixed_F_avg that is applied during initialization
  DATA_INTEGER(F_method );                          // 0 = estimate F's as free parameters, 1 = Hybrid method
  DATA_SCALAR(F_max);                               // max F = 2.56
  DATA_INTEGER(F_iterations);                       // should be between 2-5

  DATA_ARRAY(fixed_fishery_catch);                  // Observed catch for Longline fishery. dim: n_region x n_years
  DATA_ARRAY(trwl_fishery_catch);                   // Observed catch for Trawl fishery. dim: n_region x n_years

  // Selectivity indicator switches
  DATA_IVECTOR(fixed_sel_type);                        // Selectivity type for each row of ln_fixed_sel_m_pars and ln_fixed_sel_f_pars
  DATA_IVECTOR(fixed_sel_by_year_indicator);           // Selectivity time-block to apply in each model year
  DATA_IVECTOR(trwl_sel_type);                      // Selectivity type for each row of ln_trwl_sel_m_pars and ln_fixed_sel_f_pars
  DATA_IVECTOR(trwl_sel_by_year_indicator);         // Selectivity time-block to apply in each model year

  // Survey stuff
  DATA_IVECTOR(srv_dom_ll_sel_type);                // Selectivity type for each row of ln_fixed_sel_m_pars and ln_fixed_sel_f_pars
  DATA_IVECTOR(srv_dom_ll_sel_by_year_indicator);   // Selectivity time-block to apply in each model year


  // Tag-release information
  DATA_IVECTOR(tag_release_event_this_year);                 // dim: n_years.  1 = release events, 0 skip
  int n_years_with_tag_releases = sum(tag_release_event_this_year);
  DATA_ARRAY(male_tagged_cohorts_by_age);                    // numbers at age. dim: n_ages x n_region x n_years_with_tag_releases
  DATA_ARRAY(female_tagged_cohorts_by_age);                  // numbers at age. dim: n_ages x n_region x n_years_with_tag_releases
  DATA_INTEGER(n_years_to_retain_tagged_cohorts_for);        // How many years are tagged cohorts tracked in the tagged partition before they merged?
  DATA_VECTOR(initial_tag_induced_mortality);                // dim: n_years_with_tag_releases
  DATA_SCALAR(annual_tag_shedding_rate);                     // same for all regions and years -

  // Observational stuff
  DATA_MATRIX(ageing_error_matrix);                 // Ageing error/missclassification matrix n_ages x n_ages

  // Longline fishery catch at age (sex dis aggregated)
  DATA_IARRAY(fixed_catchatage_indicator);            // dim: n_regions x n_years.  1 = calculate catch at age in this year and region, 0 = don't calculate catch at age
  DATA_ARRAY(obs_fixed_catchatage);                   // Longline fishery composition observations dim = 2*n_ages x n_regions x n_years. NOTE: male first then female
  DATA_ARRAY_INDICATOR(keep_fixed_catchatage_comp, obs_fixed_catchatage); // Used for OSA residuals, when not using the multinomial likelihood
  DATA_INTEGER(fixed_catchatage_covar_structure);             // 0 = iid, 5 = AR(1), 2 = Unstructured.
  DATA_INTEGER(fixed_catchatage_comp_likelihood);             // 0 = old multinomial, 1 = TMB's multinomial. //0 = MVN (can be applied to both comp_type), 1 = Multinomial, 2 = dirichlet-multinomial
  array<Type> pred_fixed_catchatage(obs_fixed_catchatage.dim); // Sex disaggregated predicted catch at age

  // Trawl fishery catch at length (sex dis aggregated)
  DATA_IARRAY(trwl_catchatlgth_indicator);            // dim: n_regions x n_years.  1 = calculate catch at age in this year and region, 0 = don't calculate catch at age
  DATA_ARRAY(obs_trwl_catchatlgth);                   // Longline fishery composition observations dim = 2*n_length_bins x n_regions x n_years. NOTE: male first then female
  DATA_ARRAY_INDICATOR(keep_trwl_catchatlgth_comp, obs_trwl_catchatlgth); // Used for OSA residuals, when not using the multinomial likelihood
  DATA_INTEGER(trwl_catchatlgth_covar_structure);             // 0 = iid, 5 = AR(1), 2 = Unstructured.
  DATA_INTEGER(trwl_catchatlgth_comp_likelihood);             // 0 = old multinomial, 1 = TMB's multinomial. //0 = MVN (can be applied to both comp_type), 1 = Multinomial, 2 = dirichlet-multinomial
  array<Type> pred_trwl_catchatlgth(obs_trwl_catchatlgth.dim); // Sex disaggregated predicted catch at age

  // Fixed gear fishery catch at length (sex dis aggregated)
  DATA_IARRAY(fixed_catchatlgth_indicator);            // dim: n_regions x n_years.  1 = calculate catch at age in this year and region, 0 = don't calculate catch at age
  DATA_ARRAY(obs_fixed_catchatlgth);                   // Longline fishery composition observations dim = 2*n_length_bins x n_regions x n_years. NOTE: male first then female
  DATA_ARRAY_INDICATOR(keep_fixed_catchatlgth_comp, obs_fixed_catchatlgth); // Used for OSA residuals, when not using the multinomial likelihood
  DATA_INTEGER(fixed_catchatlgth_covar_structure);             // 0 = iid, 5 = AR(1), 2 = Unstructured.
  DATA_INTEGER(fixed_catchatlgth_comp_likelihood);             // 0 = old multinomial, 1 = TMB's multinomial. //0 = MVN (can be applied to both comp_type), 1 = Multinomial, 2 = dirichlet-multinomial
  array<Type> pred_fixed_catchatlgth(obs_fixed_catchatlgth.dim); // Sex disaggregated predicted catch at age

  // Longline survey catch at age
  DATA_IARRAY(srv_dom_ll_catchatage_indicator);            // dim: n_regions x n_years.  1 = calculate catch at age in this year and region, 0 = don't calculate catch at age
  DATA_ARRAY(obs_srv_dom_ll_catchatage);                   // Longline domestic survey composition observations dim = 2*n_ages x n_regions x n_years. NOTE: male first then female
  DATA_ARRAY_INDICATOR(keep_srv_dom_ll_catchatage_comp, obs_srv_dom_ll_catchatage); // Used for OSA residuals, when not using the multinomial likelihood
  DATA_INTEGER(srv_dom_ll_catchatage_covar_structure);             // 0 = iid, 5 = AR(1), 2 = Unstructured.
  DATA_INTEGER(srv_dom_ll_catchatage_comp_likelihood);             // 0 = old multinomial, 1 = TMB's multinomial. //0 = MVN (can be applied to both comp_type), 1 = Multinomial, 2 = dirichlet-multinomial
  array<Type> pred_srv_dom_ll_catchatage(obs_srv_dom_ll_catchatage.dim); // Sex disaggregated predicted catch at age

  // Longline survey biomass
  DATA_IARRAY(srv_dom_ll_bio_indicator);                      // dim: n_regions x n_years.  1 = calculate catch at age in this year and region, 0 = don't calculate observation
  DATA_ARRAY(obs_srv_dom_ll_bio);                             // Longline domestic survey biomass observations dim = n_regions x n_years.
  DATA_ARRAY(obs_srv_dom_ll_se);                              // Longline domestic survey biomass standard errors
  DATA_ARRAY_INDICATOR(keep_srv_dom_ll_bio_comp, obs_srv_dom_ll_bio); // Used for OSA residuals, when not using the multinomial likelihood
  DATA_INTEGER(srv_dom_ll_bio_comp_likelihood);               // 0 is old lognormal, 1 = dlnorm call
  DATA_INTEGER(srv_dom_ll_obs_is_abundance);                     // 1 = Abundance (Numbers), 0 = Biomass (Weight)
  array<Type> pred_srv_dom_ll_bio(obs_srv_dom_ll_bio.dim);    // Sex disaggregated predicted catch at age
  DATA_IVECTOR(srv_dom_ll_q_by_year_indicator);               // Catchability time-block to apply when deriving model predictions each year
  DATA_INTEGER(srv_dom_ll_q_transformation);                  // 0 = log, 1 = logistic (bound between 0-1)


  // Tag recovery observations Sexually disaggregated?
  // Note: All tag recoveries are assumed from the fixed fishery!!!
  DATA_IVECTOR(tag_recovery_indicator_by_year);                         // dim: n_years.  1 = calculate fitted values for tag-recoveries this year and region, 0 = don't calculate tag recovery observation for this year and region

  int n_tag_recovery_years = tag_recovery_indicator_by_year.sum();
  DATA_IARRAY(tag_recovery_indicator);
  // if tag_likelihood %in% c(0,1)
  // dim: n_tag_release_events x n_regions x n_tag_recovery_years.  1 = calculate fitted values for tag-recoveries this year and region, 0 = don't calculate tag recovery observation for this year and region
  // if tag_likelihood == 2
  // dim: n_years (indicates release year) x n_release_regions

  DATA_ARRAY(obs_tag_recovery);
  // if tag_likelihood %in% c(0,1)
  // dim: n_tag_release_events x n_regions x n_tag_recovery_years.
  // if tag_likelihood == 2
  // dim:n_recapture_events + 1 (plus one for the NC group) x n_release_regions x n_years
  DATA_INTEGER(tag_likelihood);                                 // likelihood type 0 = Poisson, 1 = Negative Binomial, 2 = Multinomial release conditioned
  array<Type> pred_tag_recovery(obs_tag_recovery.dim);


  /*
   *  Estimable parameters
   *
   */
  PARAMETER_VECTOR(ln_mean_rec);                        // Unfish equil recruitment (logged) (estimated) for each spatial region
  PARAMETER_ARRAY(trans_rec_dev);                       // Transformed Recruitment deviations. If rec_devs_sum_to_zero = 0 these are logged values, otherwise they are sum to zero values of which there are n_years-1 for each region
  PARAMETER_VECTOR(ln_init_rec_dev);                    // Initial age deviations to apply during initialization they include years before the assessment starts: length = n_init_rec_devs

  // Fishery selectivities
  PARAMETER_ARRAY(ln_fixed_sel_pars);                       // log selectivity parameters for fixed gear, dim: time-blocks:  max(sel parameters): sex
  PARAMETER_ARRAY(ln_trwl_sel_pars);                        // log selectivity parameters for Trawl gear, dim: time-blocks:  max(sel parameters): sex


  PARAMETER_ARRAY(transformed_movement_pars);               // transformed parameters for movmenet (consider both simplex and logistic? or what ever it is). dimension:  (n_regions - 1) x n_regions


  // Estimated if F_method == 0, otherwise these are derived.
  PARAMETER(ln_fixed_F_avg);                                // log average longline Fishing mortality
  PARAMETER_ARRAY(ln_fixed_F_devs);                         // Annual fishing mortality deviation dim: n_regions x n_years
  PARAMETER(ln_trwl_F_avg);                                 // log average trawl Fishing mortality
  PARAMETER_ARRAY(ln_trwl_F_devs);                          // Annual fishing mortality deviations dim: n_regions x n_years

  PARAMETER(ln_init_F_avg);                                 // log average initial Fishing mortality used when F_method == 1, else should not be estimated
  PARAMETER(ln_catch_sd);                                   // Shared across all gears
  //
  PARAMETER_ARRAY(trans_srv_dom_ll_q);                   // logistic catchabilities parameters for srv_dom_ll n_regions x n_q_time-blocks
  PARAMETER_ARRAY(ln_srv_dom_ll_sel_pars);                  // log selectivity parameters for domestic longline surveyr, dim: time-blocks:  max(sel parameters): sex
  PARAMETER_ARRAY(logistic_tag_reporting_rate);             // logistic tag-reporting dim: n_regions x n_tag_recovery_years
  // nuisance parameters
  PARAMETER(ln_tag_phi);                                    // log variance for tag data likelihood- currently only used if tag_likelihood = 1 (Negative binomial)
  PARAMETER(ln_sigma_R);                                    // standard deviation for recruitment;
  PARAMETER(ln_sigma_init_devs);                            // standard deviation for recruitment;

  // composition parameters
  PARAMETER_VECTOR(trans_trwl_catchatlgth_error);     //
  PARAMETER_VECTOR(trans_fixed_catchatlgth_error);     //
  PARAMETER_VECTOR(trans_fixed_catchatage_error);     //
  PARAMETER_VECTOR(trans_srv_dom_ll_catchatage_error);    //

  PARAMETER_VECTOR(logistic_prop_recruit_male);             //  logistic_prop_recruit_male. length: n_years

  // Initialise consistently used variables throughout the code
  int year_ndx;
  int proj_year_ndx;
  int age_ndx;
  //int len_ndx;
  int region_ndx;
  int release_region_ndx;
  int fishery_ndx;
  int tag_ndx;
  int tag_release_event_ndx = 0;
  int tag_recovery_event_ndx = 0;
  int release_year_ndx = 0;

  Type m_plus_group = 0.0;
  Type f_plus_group = 0.0;
  Type effective_sample_size = 0.0;
  Type predicted_tags;
  Type pen_posfun = 0; // this is passed to the utility posfun function and added to the likelihood as apenalty
  Type eps_for_posfun = 0.00001; // used for the posfun object to scale values above zero

  Type s1; // used in the negative binomial likelihood
  Type s2; // used in the negative binomial likelihood

  // Note: about ln_init_rec_dev - it is the opposite order to how ADMB model is formulated. I am sorry but it was easier to code.
  //       the first init_rec_dev corresponds to recruitment for age class 2 which would have arrived in styr - 1, and so on.
  // Untransform parameters
  vector<Type> mean_rec = exp(ln_mean_rec);
  Type sigma_R = exp(ln_sigma_R);
  Type sigma_init_devs = exp(ln_sigma_init_devs);
  Type sigma_init_devs_sq = sigma_init_devs * sigma_init_devs;
  Type sigma_R_sq = sigma_R * sigma_R;
  vector<Type> init_rec_dev = exp(ln_init_rec_dev);
  array<Type> recruitment_multipliers(n_regions, n_projyears);
  recruitment_multipliers.fill(1.0);
  array<Type> recruitment_devs(n_regions, n_projyears);
  recruitment_devs.fill(0.0);
  if(rec_devs_sum_to_zero == 0) {
    if(global_rec_devs == 1) {
      for(year_ndx = 0; year_ndx < trans_rec_dev.dim(1); ++year_ndx) {
        for(region_ndx = 0; region_ndx < n_regions; ++region_ndx) {
          recruitment_devs(region_ndx, year_ndx) = trans_rec_dev(0, year_ndx);
          recruitment_multipliers(region_ndx, year_ndx) = exp(recruitment_devs(region_ndx, year_ndx) - sigma_R_sq/2);
        }
      }
    } else {
      for(year_ndx = 0; year_ndx < trans_rec_dev.dim(1); ++year_ndx) {
        for(region_ndx = 0; region_ndx < n_regions; ++region_ndx) {
          recruitment_devs(region_ndx, year_ndx) = trans_rec_dev(region_ndx, year_ndx);
          recruitment_multipliers(region_ndx, year_ndx) = exp(recruitment_devs(region_ndx, year_ndx) - sigma_R_sq/2);
        }
      }
    }
  } else if(rec_devs_sum_to_zero == 1) {
    int N = trans_rec_dev.dim(1) + 1;

    if(global_rec_devs == 1) {
      Type rec_aux = 0.0;
      Type rec_multi_temp = 0.0;
      for(year_ndx = 0; year_ndx < trans_rec_dev.dim(1); ++year_ndx) {
        rec_multi_temp = rec_aux + trans_rec_dev(0, year_ndx) * Q_r_for_sum_to_zero(year_ndx);
        rec_aux += trans_rec_dev(0, year_ndx) * Q_r_for_sum_to_zero(year_ndx + N);
        for(region_ndx = 0; region_ndx < n_regions; ++region_ndx) {
          recruitment_devs(region_ndx, year_ndx) = rec_multi_temp;
          recruitment_multipliers(region_ndx, year_ndx) = exp(recruitment_devs(region_ndx, year_ndx) - sigma_R_sq/2.0); // do we need to add the -sigma^2
        }
      }
      // the last group
      for(region_ndx = 0; region_ndx < n_regions; ++region_ndx) {
        recruitment_devs(region_ndx, N - 1) = rec_aux;
        recruitment_multipliers(region_ndx, N - 1) = exp(recruitment_devs(region_ndx, N - 1) - sigma_R_sq/2.0);
      }
    } else {
      for(region_ndx = 0; region_ndx < n_regions; ++region_ndx) {
        Type rec_aux = 0.0;
        Type rec_multi_temp = 0.0;
        for(year_ndx = 0; year_ndx < trans_rec_dev.dim(1); ++year_ndx) {
          rec_multi_temp = rec_aux + trans_rec_dev(region_ndx, year_ndx) * Q_r_for_sum_to_zero(year_ndx);
          rec_aux += trans_rec_dev(region_ndx, year_ndx) * Q_r_for_sum_to_zero(year_ndx + N);
          recruitment_devs(region_ndx, year_ndx) = rec_multi_temp;
          recruitment_multipliers(region_ndx, year_ndx) = exp(recruitment_devs(region_ndx, year_ndx) - sigma_R_sq/2.0); // do we need to add the -sigma^2
        }
        // the last group
        recruitment_devs(region_ndx, N - 1) = rec_aux;
        recruitment_multipliers(region_ndx, N - 1) = exp(recruitment_devs(region_ndx, N - 1) - sigma_R_sq/2.0);
      }
    }
  }

  Type F_hist;
  if(F_method == 0) {
    F_hist = exp(ln_fixed_F_avg);
  } else {
    F_hist = exp(ln_init_F_avg);
  }
  array<Type> tag_reporting_rate(logistic_tag_reporting_rate.dim);
  for(int i = 0; i < logistic_tag_reporting_rate.dim[0]; ++i) {
    for(int j = 0; j < logistic_tag_reporting_rate.dim[1]; ++j) {
      tag_reporting_rate(i, j) = invlogit(logistic_tag_reporting_rate(i,j));
    }
  }
  // recruitment sex proportion
  vector<Type> prop_recruit_male(n_projyears);
  prop_recruit_male.fill(0.5);
  vector<Type> prop_recruit_female(n_projyears);
  prop_recruit_female.fill(0.5);
  for(int i = 0; i < logistic_prop_recruit_male.size(); ++i) {
    prop_recruit_male(i) = invlogit(logistic_prop_recruit_male(i));
    prop_recruit_female(i) = 1.0 - prop_recruit_male(i);
  }
  Type init_F_hist = F_hist * prop_F_hist;
  Type catch_sd = exp(ln_catch_sd);

  array<Type> fixed_sel_pars(ln_fixed_sel_pars.dim);
  fixed_sel_pars = exp(ln_fixed_sel_pars);

  array<Type> trwl_sel_pars(ln_trwl_sel_pars.dim);
  trwl_sel_pars = exp(ln_trwl_sel_pars);

  array<Type> srv_dom_ll_sel_pars(ln_srv_dom_ll_sel_pars.dim);
  srv_dom_ll_sel_pars = exp(ln_srv_dom_ll_sel_pars);

  array<Type> srv_dom_ll_q(trans_srv_dom_ll_q.dim);
  if(srv_dom_ll_q_transformation == 0) {
    for(int i = 0; i < srv_dom_ll_q.dim(0); ++i) { // region
      for(int j = 0; j < srv_dom_ll_q.dim(1); ++j) // time-blocks
        srv_dom_ll_q(i,j) = exp(trans_srv_dom_ll_q(i,j));
    }
  } else if(srv_dom_ll_q_transformation == 1) {
    for(int i = 0; i < srv_dom_ll_q.dim(0); ++i) { // region
      for(int j = 0; j < srv_dom_ll_q.dim(1); ++j) // time-blocks
        srv_dom_ll_q(i,j) = invlogit(trans_srv_dom_ll_q(i,j));
    }
  }


  // deal with movement
  array<Type> movement_matrix(n_regions,n_regions);                  // n_regions x n_regions. Rows sum = 1 (aka source)
  movement_matrix.fill(1.0);
  if(n_regions > 1) {
    vector<Type> cache_log_k_value(n_regions - 1);
    for(int k = 0; k < (n_regions - 1); k++)
      cache_log_k_value[k] = log(n_regions - 1 - k);

    for(region_ndx = 0; region_ndx < n_regions; ++region_ndx) {
      Type stick_length = 1.0;
      for (int k = 0; k < (n_regions - 1); ++k) {
        movement_matrix(region_ndx, k) = stick_length * invlogit(transformed_movement_pars(k, region_ndx) - cache_log_k_value(k));
        stick_length -= movement_matrix(region_ndx, k);
      }
      // plus group
      movement_matrix(region_ndx, n_regions - 1) = stick_length;
    }
  }
  Type tag_phi = exp(ln_tag_phi);

  // Parameters for dirichlet-multinomial composition
  Type theta_fixed_catchatage;
  Type theta_fixed_catchatlgth;
  Type theta_trwl_catchatlgth;
  Type theta_srv_dom_ll_catchatage;

  if(fixed_catchatage_comp_likelihood == 1)
    theta_fixed_catchatage = exp(trans_fixed_catchatage_error(0));

  if(fixed_catchatlgth_comp_likelihood == 1)
    theta_fixed_catchatlgth = exp(trans_fixed_catchatlgth_error(0));

  if(trwl_catchatlgth_comp_likelihood == 1)
    theta_trwl_catchatlgth = exp(trans_trwl_catchatlgth_error(0));

  if(srv_dom_ll_catchatage_comp_likelihood == 1)
    theta_srv_dom_ll_catchatage = exp(trans_srv_dom_ll_catchatage_error(0));

  // Declare Derived quantities
  array<Type>  SSB_yr(n_projyears, n_regions);
  vector<Type>  SSB_all_areas(n_projyears);
  array<Type>  recruitment_yr(n_projyears, n_regions);
  array<Type> init_natage_m(n_ages, n_regions);                     // Initial numbers at age Males
  array<Type> init_natage_f(n_ages, n_regions);                     // Initial numbers at age Females

  array<Type> cache_natage_m(n_ages, n_regions);                     // Initial numbers at age Males
  array<Type> cache_natage_f(n_ages, n_regions);                     // Initial numbers at age Females

  array<Type> weight_maturity_prod_f(n_ages, n_projyears);// Female weight and proportion mature, used to calcualte SSBs etc
  array<Type> natage_m(n_ages, n_regions, n_projyears + 1);          // Male numbers at age at the beginning of the year from start year to end of last projection year
  array<Type> natage_f(n_ages, n_regions, n_projyears + 1);          // Female numbers at age at the beginning of the year from start year to end of last projection year

  /*
   * the tagging partition has a slightly complex structure. For each sex we will track (n_years_to_retain_tagged_cohorts_for + 1) * n_regions release events at any point time.
   * This means that when a tagged fish are in the partition for longer than n_years_to_retain_tagged_cohorts_for we lose release information.
   * The third dimesion of tagged_natage_m and tagged_natage_f is the tag-release event id which is actually two dimensions (release year and release region)
   * To access a specific release event you need to map those two dimesions to get the correct release event. the index representation is conditional on
   * the current year. using the following notation
   *
   * t1 is tag fished released this year, t2 is tag fished released last year, t3 is tag fished released 2 years age etc
   * and
   * r1 is region 1 .... rn is the n_region
   * Then the index representation of the tagged partition follows
   * (t1-r1, t1-r2, ..., t1-rn, t2-r1, t2-r2, ..., t2-rn , ...., tn-r1, ...., tn-rn)
   * this is important when ageing and applying mortality
   * use the get_tag_release_event_ndx() function defined in AuxillaryFuns.h to retrieve the correct ndx for a given release region and release year.
   *
   * Additional Note: tagged fish are in actual numbers this is in contrast to the rest of the partition where 1 = a 1000 fish
   * you will see everytime tagged fish contribute to catch or ssb there is divide by 1000 to account for this.
   */

  array<Type> tagged_natage_m(n_ages, n_regions, (n_years_to_retain_tagged_cohorts_for + 1) * n_regions); // tagged Male number partition at age, we only retain tagged fish for n_years_to_retain_tagged_cohorts_for years then they move back to a pooled tagged group. This is to reduce computational burden
  array<Type> tagged_natage_f(n_ages, n_regions, (n_years_to_retain_tagged_cohorts_for + 1) * n_regions); // tagged Female numbers partition at age, we only retain tagged fish for n_years_to_retain_tagged_cohorts_for years then they move back to a pooled tagged group. This is to reduce computational burden
  vector<Type> temp_numbers_at_age_m(n_ages);                          // used during interim calculations for observations
  vector<Type> temp_numbers_at_age_f(n_ages);                          // used during interim calculations for observations

  vector<Type> numbers_at_age_and_sex(n_ages * 2);                   // used when calculating predicted proportions for AF observations.
  vector<Type> temp_numbers_at_age(n_ages);                          // used during interim calculations for observations
  vector<Type> temp_numbers_at_age_after_ageing_error(n_ages);       // used during interim calculations for observations
  vector<Type> temp_observed_age_and_sex(n_ages * 2);                // used to pull out observations
  vector<Type> temp_numbers_at_lgth(n_length_bins);       // used during interim calculations for observations
  vector<Type> temp_numbers_at_lgth_and_sex(n_length_bins * 2);       // used during interim calculations for observations
  vector<Type> temp_observed_lgth_and_sex(n_length_bins * 2);       // used during interim calculations for observations

  array<Type> Z_m(n_ages, n_regions, n_projyears);                   // Male total mortality at age from start year to end year
  array<Type> Z_f(n_ages, n_regions, n_projyears);                   // Female total mortality at age from start year to end year
  array<Type> S_m(n_ages, n_regions, n_projyears);                   // Male Survival at age from start year to end year
  array<Type> S_f(n_ages, n_regions, n_projyears);                   // Female Survival at age from start year to end year
  array<Type> S_m_mid(n_ages, n_regions, n_projyears);               // Male Survival at age from start year to end year
  array<Type> S_f_mid(n_ages, n_regions, n_projyears);               // Female Survival at age from start year to end year
  array<Type> F_fixed_m(n_ages, n_regions, n_projyears);                // Male Fishing mortality Longline at age from start year to end year
  array<Type> F_fixed_f(n_ages, n_regions, n_projyears);               // Female Fishing mortality Longline at age from start year to end year
  array<Type> F_trwl_m(n_ages, n_regions, n_projyears);              // Male Fishing mortality Trawl at age from start year to end year
  array<Type> F_trwl_f(n_ages, n_regions, n_projyears);              // Female Fishing mortality Trawl at age from start year to end year

  array<Type> catchatage_fixed_m(n_ages, n_regions, n_projyears);           // Male Catch at age Longline at age from start year to end year
  array<Type> catchatage_fixed_f(n_ages, n_regions, n_projyears);           // Female Catch at age Longline at age from start year to end year
  array<Type> catchatage_trwl_m(n_ages, n_regions, n_projyears);         // Male Catch at age Trawl at age from start year to end year
  array<Type> catchatage_trwl_f(n_ages, n_regions, n_projyears);         // Female Catch at age Trawl at age from start year to end year

  array<Type> annual_F_fixed(n_regions, n_projyears);                      // Fishing mortality for longline gear for each model year
  array<Type> annual_fixed_catch_pred(n_regions, n_projyears);             // Fishing mortality for longline gear for each model year
  array<Type> annual_F_trwl(n_regions, n_projyears);                    // Fishing mortality for trawl gear for each model year
  array<Type> annual_trwl_catch_pred(n_regions, n_projyears);           // Fishing mortality for trawl gear for each model year
  annual_fixed_catch_pred.fill(0.0);
  annual_trwl_catch_pred.fill(0.0);
  SSB_all_areas.setZero();

  array<Type> sel_fixed_f(n_ages, ln_fixed_sel_pars.dim(0));                  // Longline selectivity Female. dim: n_ages x n_time_blocks
  array<Type> sel_fixed_m(n_ages, ln_fixed_sel_pars.dim(0));                  // Longline selectivity Male. dim: n_age x n_time_blocks
  array<Type> sel_trwl_f(n_ages, ln_trwl_sel_pars.dim(0));                    // Trawl selectivity Female. dim: n_ages x n_projyears
  array<Type> sel_trwl_m(n_ages, ln_trwl_sel_pars.dim(0));                    // Trawl selectivity Male. dim: n_ages x n_projyears
  array<Type> sel_srv_dom_ll_f(n_ages, ln_srv_dom_ll_sel_pars.dim(0));        // Longline survey selectivity Female. dim: n_ages x n_projyears
  array<Type> sel_srv_dom_ll_m(n_ages, ln_srv_dom_ll_sel_pars.dim(0));        // Longline survey selectivity Male. dim: n_ages x n_projyears

  vector<Type> pred_recoveries_multinomial_release(n_regions * n_years_to_retain_tagged_cohorts_for + 1);
  vector<Type> obs_recoveries_multinomial_release(n_regions * n_years_to_retain_tagged_cohorts_for + 1);


  Type alpha = 0.0;                                       // alpha for the stock recruit relationship
  Type beta = 0.0;                                        // beta for the stock recruit relationship
  Type N_eff = 0.0;

  vector<Type> Bzero(n_regions);
  // Initialise Derived quantities
  /*
   * Calculate some initial Derived quantities
   */
  weight_maturity_prod_f = maturity * female_mean_weight_by_age;
  /*
   * Build selectivity objects
   */
  // TODO: Change BuildFisherySelectivity() to be like BuildSelectivity(). This will collapse
  // the number of years needed for the sel_ll_f container.
  BuildSelectivity(fixed_sel_pars.col(0), fixed_sel_type, ages, sel_fixed_m);
  BuildSelectivity(fixed_sel_pars.col(1), fixed_sel_type, ages, sel_fixed_f);
  BuildSelectivity(trwl_sel_pars.col(0), trwl_sel_type, ages, sel_trwl_m);
  BuildSelectivity(trwl_sel_pars.col(1), trwl_sel_type, ages, sel_trwl_f);
  BuildSelectivity(srv_dom_ll_sel_pars.col(0), srv_dom_ll_sel_type, ages, sel_srv_dom_ll_m);
  BuildSelectivity(srv_dom_ll_sel_pars.col(1), srv_dom_ll_sel_type, ages, sel_srv_dom_ll_f);

  // Pre-calculate F, Z and survivorship only if F_method == 0

  annual_F_fixed = exp(ln_fixed_F_avg + ln_fixed_F_devs);
  annual_F_trwl = exp(ln_trwl_F_avg + ln_trwl_F_devs);
  if(F_method == 0) {
    for(year_ndx = 0; year_ndx < n_years; ++year_ndx) {
      for(region_ndx = 0; region_ndx < n_regions; ++region_ndx) {
        for(age_ndx = 0; age_ndx < n_ages; ++age_ndx) {
          F_fixed_m(age_ndx, region_ndx, year_ndx) = annual_F_fixed(region_ndx, year_ndx) * sel_fixed_m(age_ndx, fixed_sel_by_year_indicator(year_ndx));
          F_fixed_f(age_ndx, region_ndx, year_ndx) = annual_F_fixed(region_ndx, year_ndx) * sel_fixed_f(age_ndx, fixed_sel_by_year_indicator(year_ndx));
          F_trwl_m(age_ndx, region_ndx, year_ndx) = annual_F_trwl(region_ndx, year_ndx) * sel_trwl_m(age_ndx, trwl_sel_by_year_indicator(year_ndx));
          F_trwl_f(age_ndx, region_ndx, year_ndx) = annual_F_trwl(region_ndx, year_ndx) * sel_trwl_f(age_ndx, trwl_sel_by_year_indicator(year_ndx));
          Z_f(age_ndx, region_ndx, year_ndx) = M(age_ndx, year_ndx) + F_fixed_f(age_ndx, region_ndx, year_ndx) + F_trwl_f(age_ndx, region_ndx, year_ndx);
          Z_m(age_ndx, region_ndx, year_ndx) = M(age_ndx, year_ndx) + F_fixed_m(age_ndx, region_ndx, year_ndx) + F_trwl_m(age_ndx, region_ndx, year_ndx);
          S_f(age_ndx, region_ndx, year_ndx) = exp(-1.0 * Z_f(age_ndx, region_ndx, year_ndx));
          S_m(age_ndx, region_ndx, year_ndx) = exp(-1.0 * Z_m(age_ndx, region_ndx, year_ndx));
          S_f_mid(age_ndx, region_ndx, year_ndx) = exp(-0.5 * Z_f(age_ndx, region_ndx, year_ndx));
          S_m_mid(age_ndx, region_ndx, year_ndx) = exp(-0.5 * Z_m(age_ndx, region_ndx, year_ndx));
        }
      }
    }
  }
  // These containers are only used for the hybrid F calculation i.e., F_method == 1
  vector<Type> vulnerable_bio(2);  // WHy the number 2 equals number of fisheries i.e., fixed, trawl
  vector<Type> init_popes_rate(2);
  vector<Type> catch_this_year(2);
  vector<Type> steep_jointer(2);
  vector<Type> annual_Fs(2);
  vector<Type> init_F(2);
  vector<Type> temp_Z_vals_m(n_ages);
  vector<Type> temp_Z_vals_f(n_ages);
  vector<Type> survivorship_m(n_ages);
  vector<Type> survivorship_f(n_ages);
  vector<Type> exploitation_rate(2);
  array<Type> exp_half_natural_mortality(M.dim);
  exp_half_natural_mortality = exp(-0.5 * M); // called many times in hybrid F calculation, good to cache
  Type interim_total_catch = 0;// interim catch over all fisheries in current year
  Type total_catch_this_year = 0;
  Type z_adjustment;
  Type number_of_tag_releases;

  vector<Type> nll(11); // slots
  nll.setZero();
  /* nll components
   * 0 - fixed - fishery age comp
   * 1 - trwl - fishery length comp
   * 2 - fixed - fishery length comp
   * 3 - LL domestic survey catch at age
   * 4 - LL domestic survey biomass
   * 5 - fixed fishery catch contribution
   * 6 - Trawl fishery catch contribution
   * 7 - Tag-recovery
   * 8 - Recruitment penalty/hyper prior if model is hierachical
   * 9 - init dev penalty/hyper prior if model is hierachical
   * 10 - Posfun penalty for values that must be > 0 but aren't
   */

  /*
   * Initialize the partition (age structure)
   * by running the "annual cycle" n_ages times and approximating the plus group
   * this "should" account for age accumlation along with movement
   */
  Type plus_c = 0.0;
  for(int init_iter = 0; init_iter < n_ages; ++init_iter) {
    // Equilibrium Age-structure
    // TODO: consider how to include the init devs into this calculation.
    cache_natage_f = init_natage_f;
    cache_natage_m = init_natage_m;
    for(region_ndx = 0; region_ndx < n_regions; ++region_ndx) {
      init_natage_m(0, region_ndx) = mean_rec(region_ndx)/2; //mean_rec * exp((sigma_R*sigma_R)/2)/2;
      init_natage_f(0, region_ndx) = mean_rec(region_ndx)/2;
      // Ageing + Z
      m_plus_group = cache_natage_m(n_ages - 1, region_ndx);
      f_plus_group = cache_natage_f(n_ages - 1, region_ndx);
      for(age_ndx = 1; age_ndx < n_ages; age_ndx++) {
        init_natage_f(age_ndx, region_ndx) = cache_natage_f(age_ndx - 1, region_ndx) * exp(- (M(age_ndx - 1, 0) + init_F_hist * sel_fixed_f(age_ndx, 0)));
        init_natage_m(age_ndx, region_ndx) = cache_natage_m(age_ndx - 1, region_ndx) * exp(- (M(age_ndx - 1, 0) + init_F_hist * sel_fixed_m(age_ndx, 0)));
      }
      // plus group
      init_natage_f(n_ages - 1, region_ndx) += f_plus_group *  exp(- (M(n_ages - 1, 0) + init_F_hist * sel_fixed_f(n_ages - 1, 0)));
      init_natage_m(n_ages - 1, region_ndx) += m_plus_group *  exp(- (M(n_ages - 1, 0) + init_F_hist * sel_fixed_f(n_ages - 1, 0)));
    }
    // Movement
    if(apply_fixed_movement) {
      init_natage_f = (init_natage_f.matrix() * fixed_movement_matrix).array();
      init_natage_m = (init_natage_m.matrix() * fixed_movement_matrix).array();
    } else {
      init_natage_f = (init_natage_f.matrix() * movement_matrix.matrix()).array();
      init_natage_m = (init_natage_m.matrix() * movement_matrix.matrix()).array();
    }
  }
  // Cache age-structure
  cache_natage_f = init_natage_f;
  cache_natage_m = init_natage_m;
  // Run annual cycle one more time
  for(region_ndx = 0; region_ndx < n_regions; ++region_ndx) {
    init_natage_m(0, region_ndx) = mean_rec(region_ndx)/2; //mean_rec * exp((sigma_R*sigma_R)/2)/2;
    init_natage_f(0, region_ndx) = mean_rec(region_ndx)/2;
    m_plus_group = init_natage_m(n_ages - 1, region_ndx);
    f_plus_group = init_natage_f(n_ages - 1, region_ndx);
    // Ageing + Z
    for(age_ndx = 1; age_ndx < n_ages; age_ndx++) {
      init_natage_f(age_ndx, region_ndx) = cache_natage_f(age_ndx - 1, region_ndx) * exp(- (M(age_ndx - 1, 0) + init_F_hist * sel_fixed_f(age_ndx, 0)));
      init_natage_m(age_ndx, region_ndx) = cache_natage_m(age_ndx - 1, region_ndx) * exp(- (M(age_ndx - 1, 0) + init_F_hist * sel_fixed_m(age_ndx, 0)));
    }
    // plus group
    init_natage_f(n_ages - 1, region_ndx) += f_plus_group *  exp(- (M(age_ndx - 1, 0) + init_F_hist * sel_fixed_f(n_ages - 1, 0)));
    init_natage_m(n_ages - 1, region_ndx) += m_plus_group *  exp(- (M(age_ndx - 1, 0) + init_F_hist * sel_fixed_m(n_ages - 1, 0)));
  }
  // Movement
  // Movement
  if(apply_fixed_movement) {
    init_natage_f = (init_natage_f.matrix() * fixed_movement_matrix).array();
    init_natage_m = (init_natage_m.matrix() * fixed_movement_matrix).array();
  } else {
    init_natage_f = (init_natage_f.matrix() * movement_matrix.matrix()).array();
    init_natage_m = (init_natage_m.matrix() * movement_matrix.matrix()).array();
  }
  // Approximate plus group
  for(region_ndx = 0; region_ndx < n_regions; ++region_ndx) {
    plus_c = init_natage_f(n_ages - 1, region_ndx) / cache_natage_f(n_ages - 1, region_ndx) - 1.0;
    init_natage_f(n_ages - 1, region_ndx) = cache_natage_f(n_ages - 1, region_ndx) * 1 / (1 - plus_c);
    plus_c = init_natage_m(n_ages - 1, region_ndx) / cache_natage_m(n_ages - 1, region_ndx) - 1.0;
    init_natage_m(n_ages - 1, region_ndx) = cache_natage_m(n_ages - 1, region_ndx) * 1 / (1 - plus_c);
    for(age_ndx = 0; age_ndx < n_ages; age_ndx++)
      Bzero(region_ndx) += init_natage_f(age_ndx, region_ndx) * pow(exp(- 0.5 * (M(age_ndx, 0) + init_F_hist * sel_fixed_f(age_ndx, 0))), spawning_time_proportion(0)) * weight_maturity_prod_f(age_ndx, 0);

  }
  // apply init rec devs
  if(n_init_rec_devs > 0) {
    for(region_ndx = 0; region_ndx < n_regions; ++region_ndx) {
      for(age_ndx = 1; age_ndx < (n_ages - 1); ++age_ndx) {
        if(age_ndx >= n_init_rec_devs) {
          init_natage_m(age_ndx, region_ndx) *= init_rec_dev(n_init_rec_devs - 1);
          init_natage_f(age_ndx, region_ndx) *= init_rec_dev(n_init_rec_devs - 1);
        } else {
          init_natage_m(age_ndx, region_ndx) *= init_rec_dev(age_ndx - 1);
          init_natage_f(age_ndx, region_ndx) *= init_rec_dev(age_ndx - 1);
        }
      }
    }
  }

  natage_m.col(0) = init_natage_m;
  natage_f.col(0) = init_natage_f;
  /*
   *
   *
   *
   * Run the annual cycle - main process dynmaic code
   *
   *
   *
   */
  int tag_year_counter = 0;
  int tag_recovery_counter = 0;
  for(year_ndx = 0; year_ndx < n_years; ++year_ndx) {
    // Are we releasing tags this year?
    if(tag_release_event_this_year(year_ndx) == 1) {
      for(release_region_ndx = 0; release_region_ndx < n_regions; ++release_region_ndx) {

        tag_release_event_ndx = get_tag_release_event_ndx(release_region_ndx, 0, n_regions);
        tagged_natage_m.col(tag_release_event_ndx).col(release_region_ndx) = male_tagged_cohorts_by_age.col(tag_year_counter).col(release_region_ndx);
        tagged_natage_f.col(tag_release_event_ndx).col(release_region_ndx) = female_tagged_cohorts_by_age.col(tag_year_counter).col(release_region_ndx);
        // apply initial tag_loss - applied as a mortality event
        tagged_natage_m.col(tag_release_event_ndx).col(release_region_ndx) *= exp(-initial_tag_induced_mortality(tag_year_counter));
        tagged_natage_f.col(tag_release_event_ndx).col(release_region_ndx) *= exp(-initial_tag_induced_mortality(tag_year_counter));

      }
      ++tag_year_counter;
    }

    // in each region we want to calculate recruitment, Ageing and total mortality
    for(region_ndx = 0; region_ndx < n_regions; ++region_ndx) {
      // Recruitment
      //fill in recruitment in current year - a slight ineffieciency as we have already done this for year i.e., year_ndx = 0 above but meh!
      natage_m(0, region_ndx, year_ndx) = (mean_rec(region_ndx) * recruitment_multipliers(region_ndx, year_ndx)) * prop_recruit_male(year_ndx);
      natage_f(0, region_ndx, year_ndx) = (mean_rec(region_ndx) * recruitment_multipliers(region_ndx, year_ndx)) * prop_recruit_female(year_ndx);

      recruitment_yr(year_ndx, region_ndx) = natage_m(0, region_ndx, year_ndx) + natage_f(0, region_ndx, year_ndx);

      // Calculate F and Z for this year and region using the hybrid method
      if(F_method == 1) {
        // calculate vulnerable and initial calculations
        catch_this_year(0) = fixed_fishery_catch(region_ndx, year_ndx);
        catch_this_year(1) = trwl_fishery_catch(region_ndx, year_ndx);
        total_catch_this_year = catch_this_year.sum();
        for(age_ndx = 0; age_ndx < n_ages; ++age_ndx) {
          // possibly add switches to turn off fisheries in some years
          vulnerable_bio(0) = natage_m(age_ndx, region_ndx, year_ndx) * exp_half_natural_mortality(age_ndx, year_ndx) * sel_fixed_m(age_ndx, fixed_sel_by_year_indicator(year_ndx)) * male_mean_weight_by_age(age_ndx, year_ndx);
          vulnerable_bio(0) += natage_f(age_ndx, region_ndx, year_ndx) * exp_half_natural_mortality(age_ndx, year_ndx) * sel_fixed_f(age_ndx, fixed_sel_by_year_indicator(year_ndx)) * female_mean_weight_by_age(age_ndx, year_ndx);
          vulnerable_bio(1) = natage_m(age_ndx, region_ndx, year_ndx) * exp_half_natural_mortality(age_ndx, year_ndx) * sel_trwl_m(age_ndx, trwl_sel_by_year_indicator(year_ndx)) * male_mean_weight_by_age(age_ndx, year_ndx);
          vulnerable_bio(1) += natage_f(age_ndx, region_ndx, year_ndx) * exp_half_natural_mortality(age_ndx, year_ndx) * sel_trwl_f(age_ndx, trwl_sel_by_year_indicator(year_ndx)) * female_mean_weight_by_age(age_ndx, year_ndx);
        }
        for(fishery_ndx = 0; fishery_ndx < 2; ++fishery_ndx) {
          init_popes_rate(fishery_ndx) = catch_this_year(fishery_ndx) / (vulnerable_bio(fishery_ndx) + 0.1 * catch_this_year(fishery_ndx)); //  Pope's rate  robust A.1.22 of SS appendix
          steep_jointer(fishery_ndx) = 1.0 / (1.0 + exp(30.0 * (init_popes_rate(fishery_ndx) - 0.95))); // steep logistic joiner at harvest rate of 0.95 //steep logistic joiner at harvest rate of 0.95
          exploitation_rate(fishery_ndx) = steep_jointer(fishery_ndx)  * init_popes_rate(fishery_ndx) + (1.0 - steep_jointer(fishery_ndx) ) * 0.95;
          init_F(fishery_ndx) = -log(1.0 - exploitation_rate(fishery_ndx));
        }
        // Now solve;
        for(int f_iter = 0; f_iter < F_iterations; ++f_iter) {
          interim_total_catch = 0;
          // Use calculate an initial Z
          for(age_ndx = 0; age_ndx < n_ages; ++age_ndx) {
            temp_Z_vals_m(age_ndx) = M(age_ndx, year_ndx);
            temp_Z_vals_f(age_ndx) = M(age_ndx, year_ndx);
            temp_Z_vals_m(age_ndx) += init_F(0) * sel_fixed_m(age_ndx, fixed_sel_by_year_indicator(year_ndx));
            temp_Z_vals_f(age_ndx) += init_F(0) * sel_fixed_f(age_ndx, fixed_sel_by_year_indicator(year_ndx));
            temp_Z_vals_m(age_ndx) += init_F(1) * sel_trwl_m(age_ndx, trwl_sel_by_year_indicator(year_ndx));
            temp_Z_vals_f(age_ndx) += init_F(1) * sel_trwl_f(age_ndx, trwl_sel_by_year_indicator(year_ndx));
          }
          // The survivorship is calculated as:
          for(age_ndx = 0; age_ndx < n_ages; ++age_ndx) {
            survivorship_m(age_ndx) = (1.0 - exp(-temp_Z_vals_m(age_ndx))) / temp_Z_vals_m(age_ndx);
            survivorship_f(age_ndx) = (1.0 - exp(-temp_Z_vals_f(age_ndx))) / temp_Z_vals_f(age_ndx);
          }
          // Calculate the expected total catch that would occur with the current Hrates and Z
          interim_total_catch += (natage_m.col(year_ndx).col(region_ndx).vec() * init_F(0) * sel_fixed_m.col(fixed_sel_by_year_indicator(year_ndx)).vec() * survivorship_m * male_mean_weight_by_age.col(year_ndx).vec()).sum();
          interim_total_catch += (natage_f.col(year_ndx).col(region_ndx).vec() * init_F(0) * sel_fixed_f.col(fixed_sel_by_year_indicator(year_ndx)).vec() * survivorship_f * female_mean_weight_by_age.col(year_ndx).vec()).sum();
          interim_total_catch += (natage_m.col(year_ndx).col(region_ndx).vec() * init_F(1) * sel_trwl_m.col(trwl_sel_by_year_indicator(year_ndx)).vec() * survivorship_m * male_mean_weight_by_age.col(year_ndx).vec()).sum();
          interim_total_catch += (natage_f.col(year_ndx).col(region_ndx).vec() * init_F(1) * sel_trwl_f.col(trwl_sel_by_year_indicator(year_ndx)).vec() * survivorship_f * female_mean_weight_by_age.col(year_ndx).vec()).sum();
          // make Z adjustments
          z_adjustment = total_catch_this_year / (interim_total_catch + 0.0001);

          for(age_ndx = 0; age_ndx < n_ages; ++age_ndx) {
            temp_Z_vals_m(age_ndx) = M(age_ndx, year_ndx) + z_adjustment * (temp_Z_vals_m(age_ndx) - M(age_ndx, year_ndx));
            temp_Z_vals_f(age_ndx) = M(age_ndx, year_ndx) + z_adjustment * (temp_Z_vals_f(age_ndx) - M(age_ndx, year_ndx));
            survivorship_m(age_ndx) = (1.0 - exp(-temp_Z_vals_m(age_ndx))) / temp_Z_vals_m(age_ndx);
            survivorship_f(age_ndx) = (1.0 - exp(-temp_Z_vals_f(age_ndx))) / temp_Z_vals_f(age_ndx);
          }
          // Now re-calculate a new pope rate using a vulnerable biomass based
          // on the newly adjusted F
          vulnerable_bio(0) = (natage_m.col(year_ndx).col(region_ndx).vec() * sel_fixed_m.col(fixed_sel_by_year_indicator(year_ndx)).vec() * survivorship_m * male_mean_weight_by_age.col(year_ndx).vec()).sum();
          vulnerable_bio(0) += (natage_f.col(year_ndx).col(region_ndx).vec() * sel_fixed_f.col(fixed_sel_by_year_indicator(year_ndx)).vec() * survivorship_f * female_mean_weight_by_age.col(year_ndx).vec()).sum();
          vulnerable_bio(1) = (natage_m.col(year_ndx).col(region_ndx).vec() * sel_trwl_m.col(trwl_sel_by_year_indicator(year_ndx)).vec() * survivorship_m * male_mean_weight_by_age.col(year_ndx).vec()).sum();
          vulnerable_bio(1) += (natage_f.col(year_ndx).col(region_ndx).vec() * sel_trwl_f.col(trwl_sel_by_year_indicator(year_ndx)).vec() * survivorship_f * female_mean_weight_by_age.col(year_ndx).vec()).sum();

          for(fishery_ndx = 0; fishery_ndx < 2; ++fishery_ndx) {
            exploitation_rate(fishery_ndx) = catch_this_year(fishery_ndx) / (vulnerable_bio(fishery_ndx) + 0.0001); //  Pope's rate  robust A.1.22 of SS appendix
            steep_jointer(fishery_ndx) = 1.0 / (1.0 + exp(30.0 * (exploitation_rate(fishery_ndx) - 0.95 * F_max))); // steep logistic joiner at harvest rate of 0.95 //steep logistic joiner at harvest rate of 0.95
            init_F(fishery_ndx) = steep_jointer(fishery_ndx) * exploitation_rate(fishery_ndx) + (1.0 - steep_jointer(fishery_ndx)) * F_max;
            annual_Fs(fishery_ndx) = init_F(fishery_ndx);
          }
        } // for(int f_iter = 0; f_iter < F_iterations; ++f_iter) {
        annual_F_fixed(region_ndx, year_ndx) = annual_Fs(0);
        annual_F_trwl(region_ndx, year_ndx) = annual_Fs(1);
        // cache the values
        for(age_ndx = 0; age_ndx < n_ages; ++age_ndx) {
          F_fixed_m(age_ndx, region_ndx, year_ndx) = annual_F_fixed(region_ndx, year_ndx) * sel_fixed_m(age_ndx, fixed_sel_by_year_indicator(year_ndx));
          F_fixed_f(age_ndx, region_ndx, year_ndx) = annual_F_fixed(region_ndx, year_ndx) * sel_fixed_f(age_ndx, fixed_sel_by_year_indicator(year_ndx));
          F_trwl_m(age_ndx, region_ndx, year_ndx) = annual_F_trwl(region_ndx, year_ndx) * sel_trwl_m(age_ndx, trwl_sel_by_year_indicator(year_ndx));
          F_trwl_f(age_ndx, region_ndx, year_ndx) = annual_F_trwl(region_ndx, year_ndx) * sel_trwl_f(age_ndx, trwl_sel_by_year_indicator(year_ndx));
          Z_f(age_ndx, region_ndx, year_ndx) = M(age_ndx, year_ndx) + F_fixed_f(age_ndx, region_ndx, year_ndx) + F_trwl_f(age_ndx, region_ndx, year_ndx);
          Z_m(age_ndx, region_ndx, year_ndx) = M(age_ndx, year_ndx) + F_fixed_m(age_ndx, region_ndx, year_ndx) + F_trwl_m(age_ndx, region_ndx, year_ndx);
          S_f(age_ndx, region_ndx, year_ndx) = exp(-Z_f(age_ndx, region_ndx, year_ndx));
          S_m(age_ndx, region_ndx, year_ndx) = exp(-Z_m(age_ndx, region_ndx, year_ndx));
          S_f_mid(age_ndx, region_ndx, year_ndx) = exp(-0.5 * Z_f(age_ndx, region_ndx, year_ndx));
          S_m_mid(age_ndx, region_ndx, year_ndx) = exp(-0.5 * Z_m(age_ndx, region_ndx, year_ndx));
        }
      } // if(F_method == 1) {

      // Z + Ageing
      m_plus_group = natage_m(n_ages - 1, region_ndx, year_ndx);
      f_plus_group = natage_f(n_ages - 1, region_ndx, year_ndx);
      for(age_ndx = 0; age_ndx < (n_ages - 1); age_ndx++) {
        natage_m(age_ndx + 1, region_ndx, year_ndx + 1) =  natage_m(age_ndx, region_ndx, year_ndx) * S_m(age_ndx, region_ndx, year_ndx);
        natage_f(age_ndx + 1, region_ndx, year_ndx + 1) =  natage_f(age_ndx, region_ndx, year_ndx) * S_f(age_ndx, region_ndx, year_ndx);
        SSB_yr(year_ndx, region_ndx) += natage_f(age_ndx, region_ndx, year_ndx) * pow(S_f(age_ndx, region_ndx, year_ndx), spawning_time_proportion(year_ndx)) * weight_maturity_prod_f(age_ndx, year_ndx);
        SSB_all_areas(year_ndx) += SSB_yr(year_ndx, region_ndx);
      }
      // SSB for the plus group
      SSB_yr(year_ndx, region_ndx) += natage_f(n_ages - 1, region_ndx, year_ndx) * pow(S_f(n_ages - 1, region_ndx, year_ndx), spawning_time_proportion(year_ndx)) * weight_maturity_prod_f(n_ages - 1, year_ndx);
      SSB_all_areas(year_ndx) += SSB_yr(year_ndx, region_ndx);

      // plus group accumulation
      natage_m(n_ages - 1, region_ndx, year_ndx + 1) +=  m_plus_group * S_m(n_ages - 1, region_ndx, year_ndx);
      natage_f(n_ages - 1, region_ndx, year_ndx + 1) +=  f_plus_group * S_f(n_ages - 1, region_ndx, year_ndx);
      // Calculate Catch at age
      for(age_ndx = 0; age_ndx < n_ages; age_ndx++) {
        catchatage_fixed_m(age_ndx, region_ndx, year_ndx) = F_fixed_m(age_ndx, region_ndx, year_ndx) / Z_m(age_ndx, region_ndx, year_ndx) * natage_m(age_ndx, region_ndx, year_ndx) * (1.0 - S_m(age_ndx, region_ndx, year_ndx));
        catchatage_fixed_f(age_ndx, region_ndx, year_ndx) = F_fixed_f(age_ndx, region_ndx, year_ndx) / Z_f(age_ndx, region_ndx, year_ndx) * natage_f(age_ndx, region_ndx, year_ndx) * (1.0 - S_f(age_ndx, region_ndx, year_ndx));
        catchatage_trwl_m(age_ndx, region_ndx, year_ndx) = F_trwl_m(age_ndx, region_ndx, year_ndx) / Z_m(age_ndx, region_ndx, year_ndx) * natage_m(age_ndx, region_ndx, year_ndx) * (1.0 - S_m(age_ndx, region_ndx, year_ndx));
        catchatage_trwl_f(age_ndx, region_ndx, year_ndx) = F_trwl_f(age_ndx, region_ndx, year_ndx) / Z_f(age_ndx, region_ndx, year_ndx) * natage_f(age_ndx, region_ndx, year_ndx) * (1.0 - S_f(age_ndx, region_ndx, year_ndx));
        // Account for tagged fish in Catch at age


        if(tag_year_counter > 0) {
          for(tag_ndx = 0; tag_ndx <= n_years_to_retain_tagged_cohorts_for; ++tag_ndx) {
            for(release_region_ndx = 0; release_region_ndx < n_regions; ++release_region_ndx) {
              tag_release_event_ndx = get_tag_release_event_ndx(release_region_ndx, tag_ndx, n_regions);
              catchatage_fixed_m(age_ndx, region_ndx, year_ndx) += F_fixed_m(age_ndx, region_ndx, year_ndx) / Z_m(age_ndx, region_ndx, year_ndx) * tagged_natage_m(age_ndx, region_ndx, tag_release_event_ndx) / 1000 * (1.0 - S_m(age_ndx, region_ndx, year_ndx));
              catchatage_fixed_f(age_ndx, region_ndx, year_ndx) += F_fixed_f(age_ndx, region_ndx, year_ndx) / Z_f(age_ndx, region_ndx, year_ndx) * tagged_natage_f(age_ndx, region_ndx, tag_release_event_ndx) / 1000 * (1.0 - S_f(age_ndx, region_ndx, year_ndx));
              catchatage_trwl_m(age_ndx, region_ndx, year_ndx) += F_trwl_m(age_ndx, region_ndx, year_ndx) / Z_m(age_ndx, region_ndx, year_ndx) * tagged_natage_m(age_ndx, region_ndx, tag_release_event_ndx) / 1000 * (1.0 - S_m(age_ndx, region_ndx, year_ndx));
              catchatage_trwl_f(age_ndx, region_ndx, year_ndx) += F_trwl_f(age_ndx, region_ndx, year_ndx) / Z_f(age_ndx, region_ndx, year_ndx) * tagged_natage_f(age_ndx, region_ndx, tag_release_event_ndx) / 1000 * (1.0 - S_f(age_ndx, region_ndx, year_ndx));
            }
          }
        }
        // Calculate expected catch per method.
        annual_fixed_catch_pred(region_ndx, year_ndx) += catchatage_fixed_f(age_ndx, region_ndx, year_ndx) * female_mean_weight_by_age(age_ndx, year_ndx) + catchatage_fixed_m(age_ndx, region_ndx, year_ndx) * male_mean_weight_by_age(age_ndx, year_ndx);
        annual_trwl_catch_pred(region_ndx, year_ndx) += catchatage_trwl_f(age_ndx, region_ndx, year_ndx) * female_mean_weight_by_age(age_ndx, year_ndx) + catchatage_trwl_m(age_ndx, region_ndx, year_ndx) * male_mean_weight_by_age(age_ndx, year_ndx);
      }
      // account for SSB contribution from Tagged partition

      if(tag_year_counter > 0) { // there has to be at least one tag cohort in the partition before we waste computational resources
        // add tagging components to SSB
        for(tag_ndx = 0; tag_ndx <= n_years_to_retain_tagged_cohorts_for; ++tag_ndx) {
          for(release_region_ndx = 0; release_region_ndx < n_regions; ++release_region_ndx) {
            tag_release_event_ndx = get_tag_release_event_ndx(release_region_ndx, tag_ndx, n_regions);
            for(age_ndx = 0; age_ndx < n_ages; age_ndx++) {
              SSB_yr(year_ndx, region_ndx) += tagged_natage_f(age_ndx, region_ndx, tag_release_event_ndx) / 1000 * pow(S_f(age_ndx, region_ndx, year_ndx), spawning_time_proportion(year_ndx)) * weight_maturity_prod_f(age_ndx, year_ndx);
              SSB_all_areas(year_ndx) += SSB_yr(year_ndx, region_ndx);
            }
          }
        }

        // are there any Tag-recovery observations this year and region - unfortunately they need to be calculated during the annual cycle because
        // of the way we designed the tag partition i.e., we loose release year information after n_years_to_retain_tagged_cohorts_for years
        if(tag_likelihood == 2) {
          for(tag_ndx = 0; tag_ndx < n_years_to_retain_tagged_cohorts_for; ++tag_ndx) {
            release_year_ndx = year_ndx - tag_ndx;
            // in the early years we don't want to look beyond the first year will cause a crash
            if(release_year_ndx < 0)
              continue;
            //std::cerr << "year_ndx " << year_ndx << " tag_ndx "<<  tag_ndx << " release_year_ndx " << release_year_ndx << " recovery region = " <<region_ndx << "\n";
            for(release_region_ndx = 0; release_region_ndx < n_regions; ++release_region_ndx) {
              if(tag_recovery_indicator(release_year_ndx, release_region_ndx) == 1) {
                // we have an observation for this release event calculate predicted values
                tag_recovery_event_ndx = get_tag_release_event_ndx(region_ndx, tag_ndx, n_regions);
                tag_release_event_ndx = get_tag_release_event_ndx(release_region_ndx, tag_ndx, n_regions);

                temp_numbers_at_age_m = tagged_natage_m.col(tag_release_event_ndx).col(region_ndx).vec() * F_fixed_m.col(year_ndx).col(region_ndx).vec() / Z_m.col(year_ndx).col(region_ndx).vec() * (1.0 - S_m.col(year_ndx).col(region_ndx).vec());
                numbers_at_age_and_sex.segment(0,n_ages) = temp_numbers_at_age_m;
                temp_numbers_at_age_f = tagged_natage_f.col(tag_release_event_ndx).col(region_ndx).vec() * F_fixed_f.col(year_ndx).col(region_ndx).vec() / Z_f.col(year_ndx).col(region_ndx).vec() * (1.0 - S_f.col(year_ndx).col(region_ndx).vec());
                numbers_at_age_and_sex.segment(n_ages,n_ages) = temp_numbers_at_age_f;

                numbers_at_age_and_sex *= tag_reporting_rate(region_ndx, tag_recovery_counter);

                predicted_tags = numbers_at_age_and_sex.sum();
                // posfun to make sure no zero probabilities
                predicted_tags = posfun(predicted_tags, eps_for_posfun, pen_posfun);
                pred_tag_recovery(tag_recovery_event_ndx, release_region_ndx, release_year_ndx) = predicted_tags;
              }
            }
          }
        } else {
          if(tag_recovery_indicator_by_year(year_ndx) == 1) {
            for(tag_ndx = 0; tag_ndx <= n_years_to_retain_tagged_cohorts_for; ++tag_ndx) {
              for(release_region_ndx = 0; release_region_ndx < n_regions; ++release_region_ndx) {
                tag_release_event_ndx = get_tag_release_event_ndx(release_region_ndx, tag_ndx, n_regions);
                if(tag_recovery_indicator(tag_release_event_ndx, region_ndx, tag_recovery_counter) == 1) {
                  //pred_tag_recovery
                  temp_numbers_at_age_m = tagged_natage_m.col(tag_release_event_ndx).col(region_ndx).vec() * F_fixed_m.col(year_ndx).col(region_ndx).vec() / Z_m.col(year_ndx).col(region_ndx).vec() * (1.0 - S_m.col(year_ndx).col(region_ndx).vec());
                  numbers_at_age_and_sex.segment(0,n_ages) = temp_numbers_at_age_m;
                  temp_numbers_at_age_f = tagged_natage_f.col(tag_release_event_ndx).col(region_ndx).vec() * F_fixed_f.col(year_ndx).col(region_ndx).vec() / Z_f.col(year_ndx).col(region_ndx).vec() * (1.0 - S_f.col(year_ndx).col(region_ndx).vec());
                  numbers_at_age_and_sex.segment(n_ages,n_ages) = temp_numbers_at_age_f;

                  // apply reporting rate
                  numbers_at_age_and_sex *= tag_reporting_rate(region_ndx, tag_recovery_counter);
                  predicted_tags = numbers_at_age_and_sex.sum();
                  predicted_tags = posfun(predicted_tags, eps_for_posfun, pen_posfun);

                  // likelihood contribution
                  if(tag_likelihood == 0) {
                    pred_tag_recovery(tag_release_event_ndx, region_ndx, tag_recovery_counter) = predicted_tags;
                    nll(7) -= dpois(obs_tag_recovery(tag_release_event_ndx, region_ndx, tag_recovery_counter), predicted_tags, true);
                    SIMULATE {
                      // store the simulated tag-observation in the first age-sex bin of obs_tag_recovery
                      obs_tag_recovery(tag_release_event_ndx, region_ndx, tag_recovery_counter) = rpois(predicted_tags);
                    }
                  } else if(tag_likelihood == 1) {
                    pred_tag_recovery(tag_release_event_ndx, region_ndx, tag_recovery_counter) = predicted_tags;
                    s1 = log(predicted_tags);                          // log(mu)
                    s2 = 2. * s1 - ln_tag_phi;                         // log(var - mu)
                    nll(7) -= dnbinom_robust(obs_tag_recovery(tag_release_event_ndx, region_ndx, tag_recovery_counter), s1, s2, true);
                    SIMULATE{
                      s1 = predicted_tags;
                      s2 = predicted_tags * (1.0 + tag_phi);  // (1+phi) guarantees that var >= mu
                      obs_tag_recovery(tag_release_event_ndx, region_ndx, tag_recovery_counter) = rnbinom2(s1, s2);
                    }
                  }
                }
              }
            }
          }
        } // if(tag_recovery_indicator_by_year(year_ndx) == 1) {


        // now do Z, ageing for the tagged partition
        // Sorry for anyone trying to wrap their head around this. It is complicated because
        // of how we designed the tag partition. Probably the most complex piece of code in the model.
        // Ageing means tagged fish move across an age dimension and tag-release event dimension!!!
        // This is the general algorithm
        // - apply Z and ageing in the last tagged cohort (it is a pooled tag group, by release region)
        // - apply Z and ageing in the second to last tagged cohort and add add it to the pooled tag group
        // - apply Z and ageing in the remaining tagged cohorts.
        for(release_region_ndx = 0; release_region_ndx < n_regions; ++release_region_ndx) {
          // Z + ageing to plus tag group
          tag_release_event_ndx = get_tag_release_event_ndx(release_region_ndx, n_years_to_retain_tagged_cohorts_for, n_regions);
          //std::cout << "year = " << year_ndx << " release region = " << release_region_ndx << " plus group release event ndx " << tag_release_event_ndx << "\n";

          temp_numbers_at_age_m = tagged_natage_m.col(tag_release_event_ndx).col(region_ndx);
          temp_numbers_at_age_f = tagged_natage_f.col(tag_release_event_ndx).col(region_ndx);
          // clear the partition then populate with temp container- this is to remove lingering early ages that need to get wiped as tagged fish move along the age dimension
          tagged_natage_m.col(tag_release_event_ndx).col(region_ndx).fill(0.0);
          tagged_natage_f.col(tag_release_event_ndx).col(region_ndx).fill(0.0);

          for(age_ndx = 0; age_ndx < (n_ages - 2); age_ndx++) {
            tagged_natage_m(age_ndx + 1, region_ndx, tag_release_event_ndx) = temp_numbers_at_age_m(age_ndx) * S_m(age_ndx, region_ndx, year_ndx);
            tagged_natage_f(age_ndx + 1, region_ndx, tag_release_event_ndx) = temp_numbers_at_age_f(age_ndx) * S_f(age_ndx, region_ndx, year_ndx);
          }
          tagged_natage_m(n_ages - 1, region_ndx, tag_release_event_ndx) =  temp_numbers_at_age_m(n_ages - 1) * S_m(n_ages - 1, region_ndx, year_ndx) + temp_numbers_at_age_m(n_ages - 2) * S_m(n_ages - 2, region_ndx, year_ndx);
          tagged_natage_f(n_ages - 1, region_ndx, tag_release_event_ndx) =  temp_numbers_at_age_f(n_ages - 1) * S_f(n_ages - 1, region_ndx, year_ndx) + temp_numbers_at_age_f(n_ages - 2) * S_f(n_ages - 2, region_ndx, year_ndx);
          // Now bring the second to last tag-release-cohort into the accumulation tag-cohort
          temp_numbers_at_age_m = tagged_natage_m.col(get_tag_release_event_ndx(release_region_ndx, n_years_to_retain_tagged_cohorts_for - 1, n_regions)).col(region_ndx);
          temp_numbers_at_age_f = tagged_natage_f.col(get_tag_release_event_ndx(release_region_ndx, n_years_to_retain_tagged_cohorts_for - 1, n_regions)).col(region_ndx);

          for(age_ndx = 0; age_ndx < (n_ages - 2); age_ndx++) {
            tagged_natage_m(age_ndx + 1, region_ndx, tag_release_event_ndx) += temp_numbers_at_age_m(age_ndx) * S_m(age_ndx, region_ndx, year_ndx);
            tagged_natage_f(age_ndx + 1, region_ndx, tag_release_event_ndx) += temp_numbers_at_age_f(age_ndx) * S_f(age_ndx, region_ndx, year_ndx);
          }
          tagged_natage_m(n_ages - 1, region_ndx, tag_release_event_ndx) +=  temp_numbers_at_age_m(n_ages - 1) * S_m(n_ages - 1, region_ndx, year_ndx) + temp_numbers_at_age_m(n_ages - 2) * S_m(n_ages - 2, region_ndx, year_ndx);
          tagged_natage_f(n_ages - 1, region_ndx, tag_release_event_ndx) +=  temp_numbers_at_age_f(n_ages - 1) * S_f(n_ages - 1, region_ndx, year_ndx) + temp_numbers_at_age_f(n_ages - 2) * S_f(n_ages - 2, region_ndx, year_ndx);

          // Now do the rest of the tagged cohorts going backwards
          for(tag_ndx = n_years_to_retain_tagged_cohorts_for - 1; tag_ndx > 0; --tag_ndx) {

            tag_release_event_ndx = get_tag_release_event_ndx(release_region_ndx, tag_ndx, n_regions);

            temp_numbers_at_age_m = tagged_natage_m.col(get_tag_release_event_ndx(release_region_ndx, tag_ndx - 1, n_regions)).col(region_ndx);
            temp_numbers_at_age_f = tagged_natage_f.col(get_tag_release_event_ndx(release_region_ndx, tag_ndx - 1, n_regions)).col(region_ndx);
            for(age_ndx = 0; age_ndx < (n_ages - 2); age_ndx++) {
              tagged_natage_m(age_ndx + 1, region_ndx, tag_release_event_ndx) = temp_numbers_at_age_m(age_ndx) * S_m(age_ndx, region_ndx, year_ndx);
              tagged_natage_f(age_ndx + 1, region_ndx, tag_release_event_ndx) = temp_numbers_at_age_f(age_ndx) * S_f(age_ndx, region_ndx, year_ndx);
            }
            tagged_natage_m(n_ages - 1, region_ndx, tag_release_event_ndx) =  temp_numbers_at_age_m(n_ages - 1) * S_m(n_ages - 1, region_ndx, year_ndx) + temp_numbers_at_age_m(n_ages - 2) * S_m(n_ages - 2, region_ndx, year_ndx);
            tagged_natage_f(n_ages - 1, region_ndx, tag_release_event_ndx) =  temp_numbers_at_age_f(n_ages - 1) * S_f(n_ages - 1, region_ndx, year_ndx) + temp_numbers_at_age_f(n_ages - 2) * S_f(n_ages - 2, region_ndx, year_ndx);
          }
          // Once we have aged all tag release partition, clear the first tag releases
          tag_release_event_ndx = get_tag_release_event_ndx(release_region_ndx, 0, n_regions);
          tagged_natage_m.col(tag_release_event_ndx).col(region_ndx).fill(0.0);
          tagged_natage_f.col(tag_release_event_ndx).col(region_ndx).fill(0.0);
        }
      } //if(tag_year_counter > 0)

    } // for(region_ndx = 0; region_ndx < n_regions; ++region_ndx)

    // movement for the tagged partition
    if(tag_year_counter > 0) {
      // Movement  for tagged fish and afterwards apply tag-shedding
      for(tag_ndx = 0; tag_ndx <= n_years_to_retain_tagged_cohorts_for; ++tag_ndx) {
        for(release_region_ndx = 0; release_region_ndx < n_regions; ++release_region_ndx) {
          tag_release_event_ndx = get_tag_release_event_ndx(release_region_ndx, tag_ndx, n_regions);
          if(apply_fixed_movement) {
            tagged_natage_m.col(tag_release_event_ndx) = (tagged_natage_m.col(tag_release_event_ndx).matrix() * fixed_movement_matrix).array();
            tagged_natage_f.col(tag_release_event_ndx) = (tagged_natage_f.col(tag_release_event_ndx).matrix() * fixed_movement_matrix).array();
          } else {
            tagged_natage_m.col(tag_release_event_ndx) = (tagged_natage_m.col(tag_release_event_ndx).matrix() * movement_matrix.matrix()).array();
            tagged_natage_f.col(tag_release_event_ndx) = (tagged_natage_f.col(tag_release_event_ndx).matrix() * movement_matrix.matrix()).array();
          }
          // Apply tag shedding at the end of the year which is just a mortality process
          tagged_natage_m.col(tag_release_event_ndx) *= exp(-annual_tag_shedding_rate);
          tagged_natage_f.col(tag_release_event_ndx) *= exp(-annual_tag_shedding_rate);
        }
      }
    }

    // movement
    if(apply_fixed_movement) {
      natage_m.col(year_ndx + 1) = (natage_m.col(year_ndx + 1).matrix() * fixed_movement_matrix).array();
      natage_f.col(year_ndx + 1) = (natage_f.col(year_ndx + 1).matrix() * fixed_movement_matrix).array();
    } else {
      natage_m.col(year_ndx + 1) = (natage_m.col(year_ndx + 1).matrix() * movement_matrix.matrix()).array();
      natage_f.col(year_ndx + 1) = (natage_f.col(year_ndx + 1).matrix() * movement_matrix.matrix()).array();
    }
    // increment this for the tag-recovery observation
    if(tag_recovery_indicator_by_year(year_ndx) == 1)
      ++tag_recovery_counter;

  } // for(year_ndx = 0; year_ndx < n_years; ++year_ndx) {
  /*
   * Calculate predicted values and evaluate log-likelihoods.
   */
  for(year_ndx = 0; year_ndx < n_years; ++year_ndx) {
    for(region_ndx = 0; region_ndx < n_regions; ++region_ndx) {
      // Check if we do Longline catch at age in this region and year
      if(fixed_catchatage_indicator(region_ndx, year_ndx) == 1) {
        // Get catch at age for males and account for ageing error
        temp_numbers_at_age = catchatage_fixed_m.col(year_ndx).col(region_ndx);
        temp_numbers_at_age_after_ageing_error = (temp_numbers_at_age.matrix().transpose()) * ageing_error_matrix;
        numbers_at_age_and_sex.segment(0,n_ages) = temp_numbers_at_age_after_ageing_error;
        // Now Get catch at age for females and account for ageing error
        temp_numbers_at_age = catchatage_fixed_f.col(year_ndx).col(region_ndx);
        temp_numbers_at_age_after_ageing_error = (temp_numbers_at_age.matrix().transpose()) * ageing_error_matrix;
        numbers_at_age_and_sex.segment(n_ages,n_ages) = temp_numbers_at_age_after_ageing_error;
        // do posfun so no zero proportions
        for(int i = 0; i < numbers_at_age_and_sex.size(); ++i)
          numbers_at_age_and_sex(i) = posfun(numbers_at_age_and_sex(i), eps_for_posfun, pen_posfun);
        // normalise to sum = 1 across both sexes
        numbers_at_age_and_sex /= numbers_at_age_and_sex.sum();
        // Store in predicted container
        pred_fixed_catchatage.col(year_ndx).col(region_ndx) = numbers_at_age_and_sex;
        // Get obs for likelihood calculations
        temp_observed_age_and_sex = obs_fixed_catchatage.col(year_ndx).col(region_ndx);
        // evaluate the likelihood
        if(fixed_catchatage_comp_likelihood == 0) {
          nll(0) -= dmultinom(temp_observed_age_and_sex, numbers_at_age_and_sex, true);
          SIMULATE {
            effective_sample_size = temp_observed_age_and_sex.sum();
            temp_observed_age_and_sex = rmultinom(numbers_at_age_and_sex, effective_sample_size);
            obs_fixed_catchatage.col(year_ndx).col(region_ndx) = temp_observed_age_and_sex;
          }
        } else if (fixed_catchatage_comp_likelihood == 1) {
          N_eff = sum(temp_observed_age_and_sex);
          temp_observed_age_and_sex /= N_eff;
          s1 = 0.0;
          s2 = 0.0;
          for(int ndx = 0; ndx < temp_observed_age_and_sex.size(); ndx++){
            s1 += lgamma(N_eff * temp_observed_age_and_sex(ndx) + 1);
            s2 += lgamma(N_eff * temp_observed_age_and_sex(ndx) + theta_fixed_catchatage * N_eff * numbers_at_age_and_sex(ndx)) - lgamma(theta_fixed_catchatage * N_eff * numbers_at_age_and_sex(ndx));
          }
          nll(0) -= lgamma(N_eff + 1) - s1 + lgamma(theta_fixed_catchatage * N_eff) - lgamma(N_eff + theta_fixed_catchatage * N_eff) + s2;
          SIMULATE {
            temp_observed_age_and_sex = rdirichletmulti(numbers_at_age_and_sex, N_eff, theta_fixed_catchatage);
            obs_fixed_catchatage.col(year_ndx).col(region_ndx) = temp_observed_age_and_sex;
          }
        }
      }
      // Check if we have Trawl fishery catch at age in this region and year
      if(trwl_catchatlgth_indicator(region_ndx, year_ndx) == 1) {
        // male length comp
        temp_numbers_at_lgth = (male_age_length_transition.col(year_ndx).matrix().transpose() * catchatage_trwl_m.col(year_ndx).matrix()).col(0);
        temp_numbers_at_lgth_and_sex.segment(0,n_length_bins) = temp_numbers_at_lgth;
        // female length comp
        temp_numbers_at_lgth = (female_age_length_transition.col(year_ndx).matrix().transpose() * catchatage_trwl_f.col(year_ndx).matrix()).col(0);
        temp_numbers_at_lgth_and_sex.segment(n_length_bins,n_length_bins) = temp_numbers_at_lgth;
        // do posfun so no zero proportions
        for(int i = 0; i < temp_numbers_at_lgth_and_sex.size(); ++i)
          temp_numbers_at_lgth_and_sex(i) = posfun(temp_numbers_at_lgth_and_sex(i), eps_for_posfun, pen_posfun);
        // normalise to sum = 1 across both sexes
        temp_numbers_at_lgth_and_sex /= temp_numbers_at_lgth_and_sex.sum();
        // Store in predicted container
        pred_trwl_catchatlgth.col(year_ndx).col(region_ndx) = temp_numbers_at_lgth_and_sex;
        // Get observations
        temp_observed_lgth_and_sex = obs_trwl_catchatlgth.col(year_ndx).col(region_ndx);
        // evaluate the log-likelihood
        if(trwl_catchatlgth_comp_likelihood == 0) {
          nll(1) -= dmultinom(temp_observed_lgth_and_sex, temp_numbers_at_lgth_and_sex, true);
          SIMULATE {
            effective_sample_size = temp_observed_lgth_and_sex.sum();
            temp_observed_lgth_and_sex = rmultinom(temp_numbers_at_lgth_and_sex, effective_sample_size);
            obs_trwl_catchatlgth.col(year_ndx).col(region_ndx) = temp_observed_lgth_and_sex;
          }
        } else if (trwl_catchatlgth_comp_likelihood == 1) {
          N_eff = sum(temp_observed_lgth_and_sex);
          temp_observed_lgth_and_sex /= N_eff;
          s1 = 0.0;
          s2 = 0.0;
          for(int ndx = 0; ndx < temp_observed_lgth_and_sex.size(); ndx++){
            s1 += lgamma(N_eff * temp_observed_lgth_and_sex(ndx) + 1);
            s2 += lgamma(N_eff * temp_observed_lgth_and_sex(ndx) + theta_trwl_catchatlgth * N_eff * temp_numbers_at_lgth_and_sex(ndx)) - lgamma(theta_trwl_catchatlgth * N_eff * temp_numbers_at_lgth_and_sex(ndx));
          }
          nll(1) -= lgamma(N_eff + 1) - s1 + lgamma(theta_trwl_catchatlgth * N_eff) - lgamma(N_eff + theta_trwl_catchatlgth * N_eff) + s2;
          SIMULATE {
            temp_observed_lgth_and_sex = rdirichletmulti(temp_numbers_at_lgth_and_sex, N_eff, theta_trwl_catchatlgth);
            obs_trwl_catchatlgth.col(year_ndx).col(region_ndx) = temp_observed_lgth_and_sex * N_eff;
          }
        }
      }
      // Check if we have Trawl fishery catch at age in this region and year
      if(fixed_catchatlgth_indicator(region_ndx, year_ndx) == 1) {
        // male length comp
        temp_numbers_at_lgth = (male_age_length_transition.col(year_ndx).matrix().transpose() * catchatage_fixed_m.col(year_ndx).matrix()).col(0);
        temp_numbers_at_lgth_and_sex.segment(0,n_length_bins) = temp_numbers_at_lgth;
        // female length comp
        temp_numbers_at_lgth = (female_age_length_transition.col(year_ndx).matrix().transpose() * catchatage_fixed_f.col(year_ndx).matrix()).col(0);
        temp_numbers_at_lgth_and_sex.segment(n_length_bins,n_length_bins) = temp_numbers_at_lgth;
        // do posfun so no zero proportions
        for(int i = 0; i < temp_numbers_at_lgth_and_sex.size(); ++i)
          temp_numbers_at_lgth_and_sex(i) = posfun(temp_numbers_at_lgth_and_sex(i), eps_for_posfun, pen_posfun);
        // normalise to sum = 1 across both sexes
        temp_numbers_at_lgth_and_sex /= temp_numbers_at_lgth_and_sex.sum();

        // Store in predicted container
        pred_fixed_catchatlgth.col(year_ndx).col(region_ndx) = temp_numbers_at_lgth_and_sex;
        // Get observation for ll evaluation
        temp_observed_lgth_and_sex = obs_fixed_catchatlgth.col(year_ndx).col(region_ndx);
        // evaluate the log-likelihood
        if(fixed_catchatlgth_comp_likelihood == 0) {
          nll(2) -= dmultinom(temp_observed_lgth_and_sex, temp_numbers_at_lgth_and_sex, true);
          SIMULATE {
            effective_sample_size = temp_observed_lgth_and_sex.sum();
            temp_observed_lgth_and_sex = rmultinom(temp_numbers_at_lgth_and_sex, effective_sample_size);
            obs_fixed_catchatlgth.col(year_ndx).col(region_ndx) = temp_observed_lgth_and_sex;
          }
        } else if (fixed_catchatlgth_comp_likelihood == 1) {
          N_eff = sum(temp_observed_lgth_and_sex);
          temp_observed_lgth_and_sex /= N_eff;
          s1 = 0.0;
          s2 = 0.0;
          for(int ndx = 0; ndx < temp_observed_lgth_and_sex.size(); ndx++){
            s1 += lgamma(N_eff * temp_observed_lgth_and_sex(ndx) + 1);
            s2 += lgamma(N_eff * temp_observed_lgth_and_sex(ndx) + theta_fixed_catchatlgth * N_eff * temp_numbers_at_lgth_and_sex(ndx)) - lgamma(theta_fixed_catchatlgth * N_eff * temp_numbers_at_lgth_and_sex(ndx));
          }
          nll(2) -= lgamma(N_eff + 1) - s1 + lgamma(theta_fixed_catchatlgth * N_eff) - lgamma(N_eff + theta_fixed_catchatlgth * N_eff) + s2;
          SIMULATE {
            temp_observed_lgth_and_sex = rdirichletmulti(temp_numbers_at_lgth_and_sex, N_eff, theta_fixed_catchatlgth);
            obs_fixed_catchatlgth.col(year_ndx).col(region_ndx) = temp_observed_lgth_and_sex * N_eff;
          }
        }
      }
      // Check if we have Survey Longline age composition in this region and year
      if(srv_dom_ll_catchatage_indicator(region_ndx, year_ndx) == 1) {
        // Get catch at age for males and account for ageing error
        temp_numbers_at_age = natage_m.col(year_ndx).col(region_ndx).vec() * sel_srv_dom_ll_m.col(srv_dom_ll_sel_by_year_indicator(year_ndx)).vec() * S_m_mid.col(year_ndx).col(region_ndx).vec();
        temp_numbers_at_age_after_ageing_error = (temp_numbers_at_age.matrix().transpose()) * ageing_error_matrix;
        numbers_at_age_and_sex.segment(0,n_ages) = temp_numbers_at_age_after_ageing_error;
        // Now Get catch at age for females and account for ageing error
        temp_numbers_at_age = natage_f.col(year_ndx).col(region_ndx).vec() * sel_srv_dom_ll_f.col(srv_dom_ll_sel_by_year_indicator(year_ndx)).vec() * S_f_mid.col(year_ndx).col(region_ndx).vec();
        temp_numbers_at_age_after_ageing_error = (temp_numbers_at_age.matrix().transpose()) * ageing_error_matrix;
        numbers_at_age_and_sex.segment(n_ages,n_ages) = temp_numbers_at_age_after_ageing_error;
        // do posfun so no zero proportions
        for(int i = 0; i < numbers_at_age_and_sex.size(); ++i)
          numbers_at_age_and_sex(i) = posfun(numbers_at_age_and_sex(i), eps_for_posfun, pen_posfun);
        // normalise to sum = 1 across both sexes
        numbers_at_age_and_sex /= numbers_at_age_and_sex.sum();
        // Store in predicted container
        pred_srv_dom_ll_catchatage.col(year_ndx).col(region_ndx) = numbers_at_age_and_sex;
        // Get observaiton for LL evaluation
        temp_observed_age_and_sex = obs_srv_dom_ll_catchatage.col(year_ndx).col(region_ndx);
        // evaluate the likelihood
        if(srv_dom_ll_catchatage_comp_likelihood == 0) {
          nll(3) -= dmultinom(temp_observed_age_and_sex, numbers_at_age_and_sex, true);
          SIMULATE {
            effective_sample_size = temp_observed_age_and_sex.sum();
            temp_observed_age_and_sex = rmultinom(numbers_at_age_and_sex, effective_sample_size);
            obs_srv_dom_ll_catchatage.col(year_ndx).col(region_ndx) = temp_observed_age_and_sex;
          }
        } else if (srv_dom_ll_catchatage_comp_likelihood == 1) {
          N_eff = sum(temp_observed_age_and_sex);
          temp_observed_age_and_sex /= N_eff;
          s1 = 0.0;
          s2 = 0.0;
          for(int ndx = 0; ndx < temp_observed_age_and_sex.size(); ndx++){
            s1 += lgamma(N_eff * temp_observed_age_and_sex(ndx) + 1);
            s2 += lgamma(N_eff * temp_observed_age_and_sex(ndx) + theta_srv_dom_ll_catchatage * N_eff * numbers_at_age_and_sex(ndx)) - lgamma(theta_srv_dom_ll_catchatage * N_eff * numbers_at_age_and_sex(ndx));
          }
          nll(3) -= lgamma(N_eff + 1) - s1 + lgamma(theta_srv_dom_ll_catchatage * N_eff) - lgamma(N_eff + theta_srv_dom_ll_catchatage * N_eff) + s2;
          SIMULATE {
            temp_observed_age_and_sex = rdirichletmulti(numbers_at_age_and_sex, N_eff, theta_srv_dom_ll_catchatage);
            obs_srv_dom_ll_catchatage.col(year_ndx).col(region_ndx) = temp_observed_age_and_sex;
          }
        }
      }
      // Check if we have Survey Longline biomass obsevation in this region and year
      if(srv_dom_ll_bio_indicator(region_ndx, year_ndx) == 1) {
        if(srv_dom_ll_obs_is_abundance == 1) {
          // numbers calculation
          for(age_ndx = 0; age_ndx < n_ages; age_ndx++)
            pred_srv_dom_ll_bio(region_ndx, year_ndx) += natage_m(age_ndx, region_ndx, year_ndx) * sel_srv_dom_ll_m(age_ndx, srv_dom_ll_sel_by_year_indicator(year_ndx)) * S_m_mid(age_ndx, region_ndx, year_ndx) + natage_f(age_ndx, region_ndx, year_ndx) * sel_srv_dom_ll_f(age_ndx, srv_dom_ll_sel_by_year_indicator(year_ndx)) * S_f_mid(age_ndx, region_ndx, year_ndx) ;
        } else if (srv_dom_ll_obs_is_abundance == 0) {
          // Weight calculation
          for(age_ndx = 0; age_ndx < n_ages; age_ndx++)
            pred_srv_dom_ll_bio(region_ndx, year_ndx) += natage_m(age_ndx, region_ndx, year_ndx) * sel_srv_dom_ll_m(age_ndx, srv_dom_ll_sel_by_year_indicator(year_ndx)) * male_mean_weight_by_age(age_ndx, year_ndx) * S_m_mid(age_ndx, region_ndx, year_ndx) + natage_f(age_ndx, region_ndx, year_ndx) * sel_srv_dom_ll_f(age_ndx, srv_dom_ll_sel_by_year_indicator(year_ndx)) * female_mean_weight_by_age(age_ndx, year_ndx) * S_f_mid(age_ndx, region_ndx, year_ndx) ;
        }
        // apply catchability
        pred_srv_dom_ll_bio(region_ndx, year_ndx) *= srv_dom_ll_q(region_ndx, srv_dom_ll_q_by_year_indicator(year_ndx));

        if(srv_dom_ll_bio_comp_likelihood == 0) {
          nll(4) += square((log(obs_srv_dom_ll_bio(region_ndx, year_ndx) + 0.0001) - log(pred_srv_dom_ll_bio(region_ndx, year_ndx) + 0.0001) ))/ (2.0 * square(obs_srv_dom_ll_se(region_ndx, year_ndx) / obs_srv_dom_ll_bio(region_ndx, year_ndx)));
          // not sure how best to simulate from this likelihood. I think this is right but worth having another look
          SIMULATE {
            obs_srv_dom_ll_bio(region_ndx, year_ndx) = exp(rnorm(log(pred_srv_dom_ll_bio(region_ndx, year_ndx) + 0.0001), obs_srv_dom_ll_se(region_ndx, year_ndx) / obs_srv_dom_ll_bio(region_ndx, year_ndx)));
          }
        } else if(srv_dom_ll_bio_comp_likelihood == 1) {
          nll(4) -= dlnorm(obs_srv_dom_ll_bio(region_ndx, year_ndx), log(pred_srv_dom_ll_bio(region_ndx, year_ndx)) - 0.5 * obs_srv_dom_ll_se(region_ndx, year_ndx) * obs_srv_dom_ll_se(region_ndx, year_ndx), obs_srv_dom_ll_se(region_ndx, year_ndx), true);
          SIMULATE {
            obs_srv_dom_ll_bio(region_ndx, year_ndx) = exp(rnorm(log(pred_srv_dom_ll_bio(region_ndx, year_ndx)) - 0.5 * obs_srv_dom_ll_se(region_ndx, year_ndx) * obs_srv_dom_ll_se(region_ndx, year_ndx), obs_srv_dom_ll_se(region_ndx, year_ndx)));
          }
        }
      }
      // Calculate the Not captured group if tag-likelihood is multinomial release conditioned
      if(tag_likelihood == 2) {
        if(tag_recovery_indicator(year_ndx, region_ndx) == 1) {
          pred_recoveries_multinomial_release.setZero();
          obs_recoveries_multinomial_release = obs_tag_recovery.col(year_ndx).col(region_ndx);
          number_of_tag_releases = obs_recoveries_multinomial_release.sum();
          pred_recoveries_multinomial_release = pred_tag_recovery.col(year_ndx).col(region_ndx).vec();
          // Calculate predicted proportions
          pred_recoveries_multinomial_release /= number_of_tag_releases;
          // calculate the Not recovered group
          pred_recoveries_multinomial_release(pred_recoveries_multinomial_release.size() - 1) = 1 - sum(pred_recoveries_multinomial_release);
          // Save proportions into container
          pred_tag_recovery.col(year_ndx).col(region_ndx) = pred_recoveries_multinomial_release;
          // Evaluate likelihood
          nll(7) -= dmultinom(obs_recoveries_multinomial_release, pred_recoveries_multinomial_release, true);

          SIMULATE {
            obs_recoveries_multinomial_release = rmultinom(pred_recoveries_multinomial_release, number_of_tag_releases);
            obs_tag_recovery.col(year_ndx).col(region_ndx) = obs_recoveries_multinomial_release;
          }
        }
      }
    }
  }
  /*
   * Additional objective function components that are not observations
   */
  // catch assumed to be lognormally distributed
  for(year_ndx = 0; year_ndx < n_years; ++year_ndx) {
    for(region_ndx = 0; region_ndx < n_regions; ++region_ndx) {
      nll(5) -= dnorm(log(fixed_fishery_catch(region_ndx, year_ndx)), log(annual_fixed_catch_pred(region_ndx, year_ndx)) - 0.5 * catch_sd * catch_sd, catch_sd, true);
      nll(6) -= dnorm(log(trwl_fishery_catch(region_ndx, year_ndx)), log(annual_trwl_catch_pred(region_ndx, year_ndx)) - 0.5 * catch_sd * catch_sd, catch_sd, true);
      SIMULATE {
        fixed_fishery_catch(region_ndx, year_ndx) = exp(rnorm(log(annual_fixed_catch_pred(region_ndx, year_ndx)) - 0.5 * catch_sd * catch_sd, catch_sd));
        trwl_fishery_catch(region_ndx, year_ndx) = exp(rnorm(log(annual_trwl_catch_pred(region_ndx, year_ndx)) - 0.5 * catch_sd * catch_sd, catch_sd));
      }
    }
  }

  // Recruitment Penalty
  Type n_rec_devs = 0.0;
  for(region_ndx = 0; region_ndx < trans_rec_dev.dim(0); ++region_ndx) {
    for(year_ndx = 0; year_ndx < trans_rec_dev.dim(1); ++year_ndx) {
      nll(8) += square(trans_rec_dev(region_ndx, year_ndx) - sigma_R_sq / 2.0)/(2.0* sigma_R_sq);
      n_rec_devs += 1.0;
    }
  }
  nll(8) += (n_rec_devs) * ln_sigma_R;

  // Init-dev Penalty
  if(n_init_rec_devs > 0) {
    n_rec_devs = 0.0;
    for(int i = 0; i < ln_init_rec_dev.size(); ++i) {
      nll(9) += square(ln_init_rec_dev(i) - sigma_init_devs_sq / 2.0)/(2.0* sigma_init_devs_sq);
      n_rec_devs += 1.0;
    }
    nll(9) += (n_rec_devs) * ln_sigma_init_devs;
  }
  // pos fun penalty for
  nll(10) = pen_posfun;


  /*
   *  Projection component of the model should never do this during estimation
   *  Strictly a post optimization section of code
   *
   */

  if(do_projection == 1) {
    // first populate projection elements of population containers
    for(proj_year_ndx = n_years; proj_year_ndx < n_projyears; ++proj_year_ndx) {
      for(region_ndx = 0; region_ndx < n_regions; ++region_ndx) {
        for(age_ndx = 0; age_ndx < n_ages; ++age_ndx) {
          Z_f(age_ndx, region_ndx, proj_year_ndx) = M(age_ndx, proj_year_ndx);
          Z_m(age_ndx, region_ndx, proj_year_ndx) = M(age_ndx, proj_year_ndx);
          S_f(age_ndx, region_ndx, proj_year_ndx) = exp(-1.0 * Z_f(age_ndx, region_ndx, proj_year_ndx));
          S_m(age_ndx, region_ndx, proj_year_ndx) = exp(-1.0 * Z_m(age_ndx, region_ndx, proj_year_ndx));
          S_f_mid(age_ndx, region_ndx, proj_year_ndx) = exp(-0.5 * Z_f(age_ndx, region_ndx, proj_year_ndx));
          S_m_mid(age_ndx, region_ndx, proj_year_ndx) = exp(-0.5 * Z_m(age_ndx, region_ndx, proj_year_ndx));
        }
      }
    }
    // Calculate projection period recruitment deviations

    // Deal with Future F or catch

    /*
     * Run projection annual cycle
     */
    for(proj_year_ndx = n_years; proj_year_ndx < n_projyears; ++proj_year_ndx) {
      //std::cerr << "proj year ndx = " << proj_year_ndx <<"\n";
      // in each region we want to calculate recruitment, Ageing and total mortality

      for(region_ndx = 0; region_ndx < n_regions; ++region_ndx) {
        // Recruitment
        //fill in recruitment in current year - a slight ineffieciency as we have already done this for year i.e., proj_year_ndx = 0 above but meh!
        natage_m(0, region_ndx, proj_year_ndx) = (mean_rec(region_ndx) * recruitment_multipliers(region_ndx, proj_year_ndx)) * prop_recruit_male(proj_year_ndx);
        natage_f(0, region_ndx, proj_year_ndx) = (mean_rec(region_ndx) * recruitment_multipliers(region_ndx, proj_year_ndx)) * prop_recruit_female(proj_year_ndx);

        recruitment_yr(proj_year_ndx, region_ndx) = natage_m(0, region_ndx, proj_year_ndx) + natage_f(0, region_ndx, proj_year_ndx);

        // Do something for F or catch


        // Z + Ageing
        m_plus_group = natage_m(n_ages - 1, region_ndx, proj_year_ndx);
        f_plus_group = natage_f(n_ages - 1, region_ndx, proj_year_ndx);
        for(age_ndx = 0; age_ndx < (n_ages - 1); age_ndx++) {
          natage_m(age_ndx + 1, region_ndx, proj_year_ndx + 1) =  natage_m(age_ndx, region_ndx, proj_year_ndx) * S_m(age_ndx, region_ndx, proj_year_ndx);
          natage_f(age_ndx + 1, region_ndx, proj_year_ndx + 1) =  natage_f(age_ndx, region_ndx, proj_year_ndx) * S_f(age_ndx, region_ndx, proj_year_ndx);
          SSB_yr(proj_year_ndx, region_ndx) += natage_f(age_ndx, region_ndx, proj_year_ndx) * pow(S_f(age_ndx, region_ndx, proj_year_ndx), spawning_time_proportion(proj_year_ndx)) * weight_maturity_prod_f(age_ndx, proj_year_ndx);
        }
        // SSB for the plus group
        SSB_yr(proj_year_ndx, region_ndx) += natage_f(n_ages - 1, region_ndx, proj_year_ndx) * pow(S_f(n_ages - 1, region_ndx, proj_year_ndx), spawning_time_proportion(proj_year_ndx)) * weight_maturity_prod_f(n_ages - 1, proj_year_ndx);
        // plus group accumulation
        natage_m(n_ages - 1, region_ndx, proj_year_ndx + 1) +=  m_plus_group * S_m(n_ages - 1, region_ndx, proj_year_ndx);
        natage_f(n_ages - 1, region_ndx, proj_year_ndx + 1) +=  f_plus_group * S_f(n_ages - 1, region_ndx, proj_year_ndx);
        // Calculate Catch at age
        for(age_ndx = 0; age_ndx < n_ages; age_ndx++) {
          catchatage_fixed_m(age_ndx, region_ndx, proj_year_ndx) = F_fixed_m(age_ndx, region_ndx, proj_year_ndx) / Z_m(age_ndx, region_ndx, proj_year_ndx) * natage_m(age_ndx, region_ndx, proj_year_ndx) * (1.0 - S_m(age_ndx, region_ndx, proj_year_ndx));
          catchatage_fixed_f(age_ndx, region_ndx, proj_year_ndx) = F_fixed_f(age_ndx, region_ndx, proj_year_ndx) / Z_f(age_ndx, region_ndx, proj_year_ndx) * natage_f(age_ndx, region_ndx, proj_year_ndx) * (1.0 - S_f(age_ndx, region_ndx, proj_year_ndx));
          catchatage_trwl_m(age_ndx, region_ndx, proj_year_ndx) = F_trwl_m(age_ndx, region_ndx, proj_year_ndx) / Z_m(age_ndx, region_ndx, proj_year_ndx) * natage_m(age_ndx, region_ndx, proj_year_ndx) * (1.0 - S_m(age_ndx, region_ndx, proj_year_ndx));
          catchatage_trwl_f(age_ndx, region_ndx, proj_year_ndx) = F_trwl_f(age_ndx, region_ndx, proj_year_ndx) / Z_f(age_ndx, region_ndx, proj_year_ndx) * natage_f(age_ndx, region_ndx, proj_year_ndx) * (1.0 - S_f(age_ndx, region_ndx, proj_year_ndx));

        }
        // Account for tagged fish in Catch at age

      } // for(region_ndx = 0; region_ndx < n_regions; ++region_ndx)

      // movement
      if(apply_fixed_movement) {
        natage_m.col(proj_year_ndx + 1) = (natage_m.col(proj_year_ndx + 1).matrix() * fixed_movement_matrix).array();
        natage_f.col(proj_year_ndx + 1) = (natage_f.col(proj_year_ndx + 1).matrix() * fixed_movement_matrix).array();
      } else {
        natage_m.col(proj_year_ndx + 1) = (natage_m.col(proj_year_ndx + 1).matrix() * movement_matrix.matrix()).array();
        natage_f.col(proj_year_ndx + 1) = (natage_f.col(proj_year_ndx + 1).matrix() * movement_matrix.matrix()).array();
      }

    } //     for(proj_year_ndx = n_years; proj_year_ndx < (n_years + n_projections_years); ++proj_year_ndx) {
  } //if(do_projection == 1) {






  /*
   * Report section
   */
  REPORT(nll);
  REPORT(mean_rec);
  REPORT(recruitment_multipliers);
  REPORT( recruitment_devs );
  REPORT(init_rec_dev);
  REPORT(Bzero);
  REPORT(init_natage_f);
  REPORT(init_natage_m);
  REPORT(SSB_yr);
  REPORT(natage_f);
  REPORT(natage_m);
  REPORT(init_F_hist);
  REPORT(annual_F_trwl);
  REPORT(annual_F_fixed);
  REPORT(recruitment_yr);
  REPORT(weight_maturity_prod_f);
  REPORT(catchatage_fixed_m);
  REPORT(catchatage_trwl_m);
  REPORT(catchatage_fixed_f);
  REPORT(catchatage_trwl_f);
  REPORT(annual_trwl_catch_pred);
  REPORT(annual_fixed_catch_pred);

  REPORT(movement_matrix);
  REPORT(fixed_movement_matrix);

  REPORT(fixed_sel_pars);
  REPORT(trwl_sel_pars);

  REPORT(sel_fixed_m);
  REPORT(sel_fixed_f);
  REPORT(sel_trwl_f);
  REPORT(sel_trwl_m);
  REPORT(sel_srv_dom_ll_f);
  REPORT(sel_srv_dom_ll_m);


  REPORT(srv_dom_ll_sel_pars);
  REPORT(srv_dom_ll_q);
  // estimated variances
  REPORT(sigma_R);
  REPORT(sigma_init_devs);
  REPORT( catch_sd );


  //
  REPORT(tagged_natage_m); //
  REPORT(tagged_natage_f); //

  // Catchability coeffecients
  REPORT(F_fixed_m);
  REPORT(F_fixed_f);
  REPORT(F_trwl_m);
  REPORT(F_trwl_f);

  REPORT(Z_f);
  REPORT(Z_m);

  REPORT(S_f);
  REPORT(S_m);

  REPORT(tag_reporting_rate);
  REPORT(tag_phi);

  REPORT( prop_recruit_male );
  REPORT( prop_recruit_female );

  // Report model expected/predicted values
  REPORT(pred_fixed_catchatage);
  REPORT(pred_trwl_catchatlgth);
  REPORT(pred_fixed_catchatlgth);
  REPORT(pred_srv_dom_ll_catchatage);
  REPORT(pred_srv_dom_ll_bio);
  REPORT(pred_tag_recovery);

  // Composition parameters
  REPORT( theta_fixed_catchatage);
  REPORT( theta_fixed_catchatlgth);
  REPORT( theta_trwl_catchatlgth);
  REPORT( theta_srv_dom_ll_catchatage);

  // REPORT dimensions so accesor functions and plotting functions only need the $report() object
  REPORT(n_regions);
  REPORT(ages);
  REPORT(years);
  REPORT(length_bins);
  REPORT(n_projections_years);
  REPORT( n_years_to_retain_tagged_cohorts_for );
  // Report indicators - this is mainly to aid plotting functions from report() calls
  REPORT( fixed_catchatage_indicator );
  REPORT( trwl_catchatlgth_indicator );
  REPORT( fixed_catchatlgth_indicator );
  REPORT( srv_dom_ll_catchatage_indicator );
  REPORT( srv_dom_ll_bio_indicator );
  REPORT( tag_recovery_indicator_by_year );
  REPORT( tag_recovery_indicator );
  REPORT( tag_release_event_this_year );
  REPORT( srv_dom_ll_q_transformation );
  REPORT( apply_fixed_movement );

  // likelihood types
  REPORT( tag_likelihood );
  REPORT( fixed_catchatage_comp_likelihood );
  REPORT( trwl_catchatlgth_comp_likelihood );
  REPORT( fixed_catchatlgth_comp_likelihood );
  REPORT( srv_dom_ll_catchatage_comp_likelihood );
  // Report observations
  REPORT( obs_srv_dom_ll_bio );
  REPORT( obs_srv_dom_ll_se );
  REPORT( obs_srv_dom_ll_catchatage );
  REPORT( obs_trwl_catchatlgth );
  REPORT( obs_fixed_catchatage );
  REPORT( obs_fixed_catchatlgth );
  REPORT( fixed_fishery_catch );
  REPORT( trwl_fishery_catch );
  REPORT( obs_tag_recovery );

  // AD reports this will report standard errors for these quantities
  // using TMB::sdreport() method
  ADREPORT(tag_reporting_rate);
  ADREPORT(SSB_yr);
  ADREPORT(SSB_all_areas);
  ADREPORT(movement_matrix);
  ADREPORT(recruitment_yr);
  ADREPORT(annual_F_fixed);
  ADREPORT(annual_F_trwl);
  ADREPORT(init_F_hist);

  ADREPORT(sel_fixed_m);
  ADREPORT(sel_fixed_f);
  ADREPORT(sel_trwl_f);
  ADREPORT(sel_trwl_m);
  ADREPORT(sel_srv_dom_ll_f);
  ADREPORT(sel_srv_dom_ll_m);

  ADREPORT(sigma_R);
  ADREPORT(sigma_init_devs);
  ADREPORT( catch_sd );

  // REMOVE objects after this comment.
  // I created them for reporting interim calculations debugging etc

  return nll.sum();
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif
