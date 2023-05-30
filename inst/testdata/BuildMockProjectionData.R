#'
#' Build MockSablefishModel that is used to validate functionality
#' This R-script will build a simple Sablefish like model and save the data and parameters that can be called
#' by unit-tests. There are no observations in this model.
#'
set.seed(234)
#library(SpatialSablefishAssessment)
years = 2010:2020
min_age = 2
max_age = 31
ages = min_age:max_age
region_key = data.frame(area = c("BS","AI","WGOA","CGOA","EGOA"), TMB_ndx = c(0:4))
length_bins = seq(from = 41, to = 99,by = 2)
M = 0.1
sex_length_lvls = paste0(rep(c("M","F"), each = length(length_bins)), length_bins)
reg_lvls = region_key$area[region_key$TMB_ndx + 1] ## need to be in correct order
year_lvls = years ## need to be in correct order
get_tag_release_ndx = function(region_ndx, release_event_year_ndx, n_regions) {
  return ((release_event_year_ndx - 1) * n_regions + (region_ndx - 1) + 1);
}
#'
#' TMB data
#'
data <- list()
data$ages = ages
data$years = years
data$length_bins = length_bins
data$n_regions = nrow(region_key)
n_regions = data$n_regions
n_ages = length(data$ages)
n_surveys = 1
data$n_surveys = n_surveys
n_length_bins = length(data$length_bins) # the last length bin value is the minimum for a length plus group
data$n_projections_years = 10
data$do_projection = 1
n_projyears = length(data$years) +  data$n_projections_years
n_years = length(data$years)
projyears = min(data$years):(max(data$years) + data$n_projections_years)

data$global_rec_devs = 1
data$rec_devs_sum_to_zero = 0
#' used for sum to zero constraint
Q_sum_to_zero_QR <- function(N) {
  Q_r = vector(length = N * 2);

  for(i in 1:(N - 1)) {
    Q_r[i] = -sqrt((N-i)/(N-i+1.0));
    Q_r[i+N] = 1.0 / sqrt((N-i) * (N-i+1));
  }
  return (Q_r);
}

data$Q_r_for_sum_to_zero = Q_sum_to_zero_QR(n_years);
data$n_init_rec_devs = 0
data$M = matrix(0.104884, nrow = n_ages, ncol = n_projyears)
maturity = c(0.02,0.05,0.09,0.18,0.31,0.49,0.67,0.81,0.9,0.95,0.98,0.99,0.99,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
data$maturity = matrix(maturity, nrow = n_ages, ncol = n_projyears, byrow = F)

## weight at age
weight_at_age_f = c(1.13,1.57,2.02,2.46,2.88,3.27,3.61,3.93,4.2,4.44,4.65,4.83,4.99,5.12,5.24,5.34,5.42,5.49,5.55,5.6,5.64,5.68,5.71,5.73,5.76,5.77,5.79,5.8,5.81,5.85)
weight_at_age_m = c(1.07,1.44,1.78,2.07,2.31,2.5,2.66,2.79,2.89,2.96,3.02,3.07,3.1,3.13,3.15,3.17,3.18,3.19,3.2,3.2,3.21,3.21,3.21,3.21,3.22,3.22,3.22,3.22,3.22,3.22)

# turn to matrix
weight_at_age_f_mat = matrix(weight_at_age_f, byrow = T, nrow = n_projyears, ncol = n_ages)
weight_at_age_m_mat = matrix(weight_at_age_m, byrow = T, nrow = n_projyears, ncol = n_ages)
##
data$female_mean_weight_by_age = t(weight_at_age_f_mat)
data$male_mean_weight_by_age = t(weight_at_age_m_mat)

## age length transition matrix
data$female_age_length_transition = data$male_age_length_transition = array(0, dim = c(n_ages, n_length_bins, n_projyears), dimnames = list(data$ages, data$length_bins, projyears))

m_ALK = matrix(c(0.04,0.01,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.06,0.02,0.01,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.11,0.04,0.01,0.01,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.16,0.08,0.03,0.02,0.01,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.18,0.12,0.07,0.04,0.02,0.01,0.01,0.01,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.17,0.16,0.11,0.07,0.04,0.03,0.02,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.13,0.17,0.15,0.11,0.07,0.05,0.04,0.03,0.02,0.02,0.02,0.02,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.08,0.15,0.17,0.14,0.11,0.09,0.07,0.05,0.05,0.04,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.03,0.03,0.04,0.11,0.16,0.16,0.15,0.12,0.1,0.09,0.08,0.07,0.06,0.06,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.04,0.04,0.04,0.04,0.04,0.04,0.04,0.04,0.04,0.04,0.04,0.02,0.07,0.13,0.16,0.16,0.15,0.14,0.12,0.11,0.1,0.09,0.09,0.08,0.08,0.08,0.08,0.07,0.07,0.07,0.07,0.07,0.07,0.07,0.07,0.07,0.07,0.07,0.07,0.07,0.07,0.01,0.03,0.08,0.13,0.15,0.16,0.16,0.15,0.14,0.13,0.13,0.12,0.12,0.11,0.11,0.11,0.11,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0,0.01,0.05,0.09,0.12,0.14,0.15,0.15,0.15,0.15,0.15,0.14,0.14,0.14,0.14,0.13,0.13,0.13,0.13,0.13,0.13,0.13,0.13,0.12,0.12,0.12,0.12,0.12,0.12,0.12,0,0,0.02,0.05,0.08,0.11,0.13,0.14,0.14,0.15,0.15,0.15,0.15,0.15,0.15,0.14,0.14,0.14,0.14,0.14,0.14,0.14,0.14,0.14,0.14,0.14,0.14,0.14,0.14,0.14,0,0,0.01,0.02,0.05,0.07,0.09,0.1,0.12,0.12,0.13,0.13,0.13,0.14,0.14,0.14,0.14,0.14,0.14,0.14,0.14,0.14,0.14,0.14,0.14,0.13,0.13,0.13,0.13,0.13,0,0,0,0.01,0.02,0.04,0.05,0.07,0.08,0.09,0.1,0.1,0.11,0.11,0.11,0.11,0.12,0.12,0.12,0.12,0.12,0.12,0.12,0.12,0.12,0.12,0.12,0.12,0.12,0.12,0,0,0,0,0.01,0.02,0.03,0.04,0.05,0.06,0.06,0.07,0.07,0.08,0.08,0.08,0.09,0.09,0.09,0.09,0.09,0.09,0.09,0.09,0.09,0.09,0.09,0.09,0.09,0.09,0,0,0,0,0,0.01,0.01,0.02,0.03,0.03,0.04,0.04,0.05,0.05,0.05,0.05,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0,0,0,0,0,0,0,0.01,0.01,0.02,0.02,0.02,0.02,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.04,0.04,0.04,0.04,0.04,0.04,0.04,0.04,0.04,0.04,0,0,0,0,0,0,0,0,0,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0,0,0,0,0,0,0,0,0,0,0,0,0,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
              byrow = F, nrow = n_ages, ncol = n_length_bins)
f_ALK = matrix(c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.01,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.03,0.01,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.07,0.02,0.01,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.13,0.04,0.01,0.01,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.19,0.08,0.03,0.01,0.01,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.2,0.12,0.06,0.03,0.02,0.01,0.01,0.01,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.01,0.17,0.15,0.09,0.05,0.03,0.02,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.11,0.17,0.12,0.07,0.05,0.03,0.02,0.02,0.02,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.05,0.15,0.14,0.1,0.07,0.05,0.04,0.03,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.12,0.15,0.12,0.09,0.07,0.06,0.05,0.04,0.03,0.03,0.03,0.03,0.03,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.03,0.01,0.08,0.13,0.14,0.12,0.09,0.08,0.06,0.05,0.05,0.04,0.04,0.04,0.04,0.04,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0,0.04,0.11,0.13,0.13,0.11,0.1,0.08,0.07,0.07,0.06,0.06,0.05,0.05,0.05,0.05,0.05,0.05,0.04,0.04,0.04,0.04,0.04,0.04,0.04,0.04,0.04,0.04,0.04,0.04,0,0.02,0.07,0.11,0.13,0.12,0.11,0.1,0.09,0.08,0.08,0.07,0.07,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0,0.01,0.04,0.09,0.11,0.12,0.12,0.11,0.1,0.1,0.09,0.08,0.08,0.08,0.08,0.07,0.07,0.07,0.07,0.07,0.07,0.07,0.07,0.07,0.06,0.06,0.06,0.06,0.06,0.06,0,0,0.02,0.06,0.09,0.11,0.11,0.11,0.11,0.1,0.1,0.1,0.09,0.09,0.09,0.08,0.08,0.08,0.08,0.08,0.08,0.08,0.08,0.07,0.07,0.07,0.07,0.07,0.07,0.07,0,0,0.01,0.04,0.07,0.09,0.1,0.11,0.11,0.11,0.1,0.1,0.1,0.1,0.09,0.09,0.09,0.09,0.09,0.09,0.08,0.08,0.08,0.08,0.08,0.08,0.08,0.08,0.08,0.08,0,0,0,0.02,0.04,0.07,0.08,0.09,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.09,0.09,0.09,0.09,0.09,0.09,0.09,0.09,0.09,0.08,0.08,0.08,0.08,0.08,0.08,0,0,0,0.01,0.03,0.05,0.06,0.07,0.08,0.09,0.09,0.09,0.09,0.09,0.09,0.09,0.09,0.09,0.09,0.09,0.09,0.09,0.09,0.09,0.08,0.08,0.08,0.08,0.08,0.08,0,0,0,0,0.01,0.03,0.04,0.06,0.07,0.07,0.08,0.08,0.08,0.08,0.08,0.08,0.08,0.08,0.08,0.08,0.08,0.08,0.08,0.08,0.08,0.08,0.08,0.08,0.08,0.08,0,0,0,0,0.01,0.02,0.03,0.04,0.05,0.06,0.06,0.07,0.07,0.07,0.07,0.07,0.07,0.07,0.07,0.07,0.07,0.07,0.07,0.07,0.07,0.07,0.07,0.07,0.07,0.07,0,0,0,0,0,0.01,0.02,0.02,0.03,0.04,0.05,0.05,0.05,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0,0,0,0,0,0,0.01,0.01,0.02,0.03,0.03,0.04,0.04,0.04,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0,0,0,0,0,0,0,0.01,0.01,0.02,0.02,0.02,0.03,0.03,0.03,0.03,0.04,0.04,0.04,0.04,0.04,0.04,0.04,0.04,0.04,0.04,0.04,0.04,0.04,0.04,0,0,0,0,0,0,0,0,0.01,0.01,0.01,0.02,0.02,0.02,0.02,0.02,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0,0,0,0,0,0,0,0,0,0.01,0.01,0.01,0.01,0.01,0.01,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.03,0,0,0,0,0,0,0,0,0,0,0,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01),
               byrow = F, nrow = n_ages, ncol = n_length_bins)
data$male_age_length_transition[,,] = m_ALK
data$female_age_length_transition[,,] = f_ALK

## Movement - this is used to calculate the simplex for parameters. Cannot have a 1 and 0s
data$movement_matrix = data$fixed_movement_matrix = matrix(0, nrow = n_regions, ncol = n_regions);
diag(data$fixed_movement_matrix) = 1
diag(data$movement_matrix) = 0.9
data$movement_matrix = data$movement_matrix + rlnorm(n = n_regions * n_regions, log(0.01), 0.1)
# renormalise
data$movement_matrix = sweep(data$movement_matrix, 1, STATS = rowSums(data$movement_matrix), "/")

data$spawning_time_proportion = rep(0, n_projyears)
data$sigma_R = 1.2
data$SrType = 3
data$apply_fixed_movement = 1; ##
data$do_recruits_move = 0
#'
#' Fishing inputs
#'
data$F_method = 0
data$F_max = 3
data$F_iterations = 4
##
data$prop_F_hist = 0.0
# drop years outside model years

data$trwl_fishery_catch = matrix(rlnorm(n_regions * n_projyears, log(1), 0.3), nrow = n_regions, ncol = n_years)
data$fixed_fishery_catch = matrix(rlnorm(n_regions * n_projyears, log(1), 0.3), nrow = n_regions, ncol = n_years)

## control variables
data$fixed_sel_type = as.vector(rep(0, 1), mode = "integer")
data$fixed_sel_by_year_indicator = as.vector(rep(0, n_projyears), mode = "integer")

data$trwl_sel_type = as.vector(rep(1, 1), mode = "integer")
data$trwl_sel_by_year_indicator = as.vector(rep(0, n_projyears), mode = "integer")

data$srv_sel_type = matrix(rep(0, 1), ncol = data$n_surveys)
data$srv_sel_by_year_indicator = matrix(rep(0, n_projyears), ncol = data$n_surveys)

#'
#' Tag release stuff
#'
data$apply_Z_on_tagged_fish = 0
data$apply_fishery_tag_reporting = 0;
data$apply_tag_reporting_rate = 0;
tag_release_years = c(2010)
data$tag_release_event_this_year = rep(0, n_years) ## no tag releases
data$tag_release_event_this_year[data$years %in% tag_release_years] = 1
data$male_tagged_cohorts_by_age = array(50, dim = c(n_ages, n_regions, length(tag_release_years)))
data$female_tagged_cohorts_by_age = array(50, dim = c(n_ages, n_regions, length(tag_release_years)))

data$n_years_to_retain_tagged_cohorts_for = 6
data$initial_tag_induced_mortality = rep(0.0, length(tag_release_years))
data$annual_tag_shedding_rate = 0.0

#'
#' Observation data
#'
#data$ageing_error_matrix = sab_inputs$ageing_error
data$ageing_error_matrix = matrix(0, nrow = n_ages, ncol = n_ages)
diag(data$ageing_error_matrix) = 1.0

# reformat year and region into factos so we can keep ordering

### Fixed gear fishery AF
data$fixed_catchatage_indicator = matrix(0, nrow = n_regions, ncol = n_years)
data$obs_fixed_catchatage = array(5, dim = c(n_ages * 2, n_regions, n_years))

data$fixed_catchatage_covar_structure = 0
data$fixed_catchatage_comp_likelihood = 0

### Trawl gear fishery LF
data$trwl_catchatlgth_indicator = matrix(0, nrow = n_regions, ncol = n_years)
data$obs_trwl_catchatlgth = array(5, dim = c(n_length_bins * 2, n_regions, n_years))

data$trwl_catchatlgth_covar_structure = 0
data$trwl_catchatlgth_comp_likelihood = 0

### Fixed gear fishery LF
data$fixed_catchatlgth_indicator = matrix(0, nrow = n_regions, ncol = n_years)
data$obs_fixed_catchatlgth = array(5, dim = c(n_length_bins * 2, n_regions, n_years))
data$fixed_catchatlgth_covar_structure = 0
data$fixed_catchatlgth_comp_likelihood = 0




### Survey LL proportions at age
### Survey LL proportions at age
data$srv_catchatage_indicator = array(0, dim = c(n_regions, n_years, n_surveys))
data$obs_srv_catchatage = array(5, dim = c(n_ages * 2, n_regions, n_years, n_surveys))

data$srv_catchatage_covar_structure = rep(1, n_surveys)
data$srv_catchatage_comp_likelihood = rep(0, n_surveys)

### Survey LL index
data$srv_bio_indicator = array(0, dim = c(n_regions, n_years, n_surveys))
data$obs_srv_bio = array(1, dim = c(n_regions, n_years,n_surveys))
data$obs_srv_se = array(0.2, dim = c(n_regions, n_years, n_surveys))

data$srv_bio_likelihood = rep(1, n_surveys)
data$srv_obs_is_abundance = rep(1, n_surveys)
data$srv_q_by_year_indicator = matrix(0, nrow = n_years, ncol = n_surveys)
data$srv_q_transformation = rep(1, n_surveys) ## logistic
data$q_is_nuisance = rep(0, n_surveys)

tag_recovery_years = 2011:2020
# drop any recovery years before release years
#tag_recovery_years = tag_recovery_years[which(tag_recovery_years %in% (tag_release_years + 1))] ## the plus one is because we don't allow a recovery unless after a year at release
#
data$tag_recovery_indicator_by_year = rep(0, n_years) ## no tag releases
data$tag_recovery_indicator_by_year[data$years %in% tag_recovery_years] = 1
data$tag_recovery_indicator = array(0, dim = c(n_regions * (data$n_years_to_retain_tagged_cohorts_for + 1), n_regions, length(tag_recovery_years)))
data$obs_tag_recovery = array(10, dim = c(n_regions * (data$n_years_to_retain_tagged_cohorts_for + 1), n_regions, length(tag_recovery_years)))

## track each tag cohort for 6 years opst release
for(y_ndx in 1:length(tag_release_years)) {
  for(r_ndx in 1:data$n_regions) {
    for(possible_recovery_years in 1:(data$n_years_to_retain_tagged_cohorts_for)) {
      for(possible_recovery_region in 1:data$n_regions) {
        tmp_recovery_year = tag_release_years[y_ndx] + possible_recovery_years
        if(!tmp_recovery_year %in% tag_recovery_years)
          next;
        recovery_year_ndx = which(tag_recovery_years %in% tmp_recovery_year)
        release_event_ndx = get_tag_release_ndx(r_ndx, possible_recovery_years + 1, data$n_regions)
        data$tag_recovery_indicator[release_event_ndx, possible_recovery_region, recovery_year_ndx] = 1

      }
    }
  }
}

data$tag_likelihood = 0
data$evaluate_tag_likelihood = 1

data$future_recruitment_type = 0
data$year_ndx_for_empirical_resampling = c(0,n_years - 1)
data$future_fishing_type = 0
data$future_fishing_inputs_trwl = array(0.1, dim = c(data$n_regions, data$n_projections_years))
data$future_fishing_inputs_fixed = array(0.1, dim = c(data$n_regions, data$n_projections_years))

#' logit bounds X which is between 0-1 to -inf -> inf based on the logit transformation
#' equivalent to qlogis(X)
#' @param X scalar range [0,1]
#' @export
#' @return Y to be between -inf, inf
logit = function(X) {
  log(X / (1 - X))
}

#' simplex
#' takes a unit vector (or a vector that will be scaled by the mean) of length K and converts to unconstrained K - 1 vector
#' @param xk vector of length K - 1 takes unconstrained values
#' @param sum_to_one whether to rescale xk so it sums to one
#' @return vector of length K - 1 unconstrained values
#' @export
#'
simplex <- function (xk, sum_to_one = TRUE)  {
  zk = vector()
  if (!sum_to_one) {
    xk = xk/sum(xk)
  }
  else {
    if (abs(sum(xk) - 1) > 0.001)
      stop("xk needs to sum = 1, otherwise speify sum_to_one = TRUE")
  }
  K = length(xk)
  zk[1] = xk[1]/(1)
  for (k in 2:(K - 1)) {
    zk[k] = xk[k]/(1 - sum(xk[1:(k - 1)]))
  }
  yk = stats::qlogis(zk) - log(1/(K - 1:(K - 1)))
  return(yk)
}
#'
#' TMB Parameter definition OM values or starting values for EM
#'
parameters <- list()
parameters$ln_mean_rec = rnorm(data$n_regions, log(14), 0.3)

## Fishery selectivity
parameters$ln_fixed_sel_pars = array(0, dim = c(1, 2, 2))
parameters$ln_trwl_sel_pars = array(0, dim = c(1, 2, 2))
parameters$ln_srv_sel_pars = array(0, dim = c(1, 2, 2, 1))

## populate parameters Note some of the male delta values are set to the female values. Line 1800 tem.tpl
parameters$ln_fixed_sel_pars[1,1,1] = 2.111
parameters$ln_fixed_sel_pars[1,2,1] = -0.7118
parameters$ln_fixed_sel_pars[1,1,2] = 1.576
parameters$ln_fixed_sel_pars[1,2,2] = -0.7118

parameters$ln_trwl_sel_pars[1,1,1] = 2.311
parameters$ln_trwl_sel_pars[1,2,1] = 2.2150
parameters$ln_trwl_sel_pars[1,1,2] = 2.011
parameters$ln_trwl_sel_pars[1,2,2] = 2.2150

## NOTE: all delta parameters for all survey selectivities are fixed based on srv_dom_ll_1
parameters$ln_srv_sel_pars[1,1,1, 1] = 2.111
parameters$ln_srv_sel_pars[1,2,1, 1] = -0.711
parameters$ln_srv_sel_pars[1,1,2, 1] = 1.576
parameters$ln_srv_sel_pars[1,2,2, 1] = -0.7118

## movement pars
parameters$transformed_movement_pars = matrix(NA, nrow = n_regions - 1, ncol = n_regions)
for(i in 1:n_regions)
  parameters$transformed_movement_pars[,i] = simplex(data$movement_matrix[i,])

parameters$ln_fixed_F_avg = -2.965016
parameters$ln_fixed_F_devs = array(0, dim = c(n_regions, n_years))

parameters$ln_trwl_F_avg = -2.965016
parameters$ln_trwl_F_devs = array(0, dim = c(n_regions, n_years))

parameters$ln_init_F_avg = parameters$ln_fixed_F_avg
parameters$trans_srv_q = array(logit(0.2), dim = c(data$n_regions, length(unique(data$srv_q_by_year_indicator[,1])), n_surveys))
parameters$trans_rec_dev = array(0, dim = c(1, n_years))
parameters$ln_init_rec_dev = 0
parameters$ln_catch_sd = log(0.02)
parameters$ln_sigma_init_devs = log(0.2)
parameters$ln_sigma_R = log(data$sigma_R)
parameters$trans_SR_pars = log(0.8)
parameters$trans_trwl_catchatlgth_error = log(1)
parameters$trans_fixed_catchatlgth_error = log(1)
parameters$trans_fixed_catchatage_error = log(1)
parameters$trans_srv_catchatage_error = rep(log(1), n_surveys)
parameters$logistic_prop_recruit_male = rep(0, length(data$years))

## reporting rate
parameters$logistic_tag_reporting_rate = array(logit(0.999), dim = c(n_regions, length(tag_recovery_years)))
parameters$ln_tag_phi = log(1)
data$model = "TagIntegratedValidate"
save(data, parameters, region_key, file = file.path("inst", "testdata", "MockProjectionData.RData"))

## rm functions so we don't get a namespace clash
rm(list = c("get_tag_release_ndx", "logit", "simplex","Q_sum_to_zero_QR"))

validate_input_data_and_parameters(data, parameters)

########################
## Check this model data parameter combo doesn't cause issues when making the AD object
########################
data$model = "TagIntegratedValidate"
modA <- TMB::MakeADFun(data = data,
                       parameters = parameters,
                       DLL = "SpatialSablefishAssessment_TMBExports")

data$model = "TagIntegrated"
modA <- TMB::MakeADFun(data = data,
                       parameters = parameters,
                       DLL = "SpatialSablefishAssessment_TMBExports")
