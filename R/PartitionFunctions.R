#
# Partition functions
#
#


#' get_partition
#' @param MLE_report a list that is output from obj$report() usually once an optimsation routine has been done.
#' @param region_key data.frame with colnames area and TMB_ndx for providing real region names to objects
#' @return data frame of parition
#' @export
get_partition = function(MLE_report, region_key = NULL) {
  years = MLE_report$years
  projyears = min(MLE_report$years):(max(MLE_report$years) + MLE_report$n_projections_years)
  projyears = c(projyears, max(projyears) + 1) ## add an extra year for the extra year
  regions = 1:MLE_report$n_regions
  ages = MLE_report$ages

  dimnames(MLE_report$natage_f) = dimnames(MLE_report$natage_m) = list(ages, regions, projyears)
  f_natage = reshape2::melt(MLE_report$natage_f)
  m_natage = reshape2::melt(MLE_report$natage_m)

  colnames(f_natage) = colnames(m_natage) = c("Age", "Region", "Year", "Numbers")
  f_natage$sex = "Female"
  m_natage$sex = "Male"

  full_df = rbind(f_natage, m_natage)
  if(is.null(region_key)) {
    full_df$Region = paste0("Region ", full_df$Region)
  } else {
    full_df$Region = region_key$area[match(full_df$Region, (region_key$TMB_ndx + 1))]
  }
  return(full_df);
}

#'
#' plot_init_nage plot initial numbers at age
#' @param MLE_report a list that is output from obj$report() usually once an optimsation routine has been done.
#' @param data list that is passed to the MakeADfun for the TMB model
#' @param region_key data.frame with colnames area and TMB_ndx for providing real region names to objects
#' @return ggplot2 object that will plot if an observation occurs in a year and region
#' @export
plot_init_nage = function(MLE_report, region_key = NULL) {
  years = MLE_report$years
  regions = 1:MLE_report$n_regions
  ages = MLE_report$ages

  dimnames(MLE_report$init_natage_f) = dimnames(MLE_report$init_natage_m) = list(ages, regions)
  f_init_age = reshape2::melt(MLE_report$init_natage_f)
  m_init_age = reshape2::melt(MLE_report$init_natage_m)

  colnames(f_init_age) = colnames(m_init_age) = c("Age", "Region", "Numbers")
  f_init_age$sex = "Female"
  m_init_age$sex = "Male"

  full_df = rbind(f_init_age, m_init_age)
  if(is.null(region_key)) {
    full_df$Region = paste0("Region ", full_df$Region)
  } else {
    full_df$Region = region_key$area[match(full_df$Region, (region_key$TMB_ndx + 1))]
  }
  gplt = ggplot(full_df) +
    geom_line(aes(x = Age, y = Numbers, col = sex, linetype = sex), linewidth= 1.1) +
    guides( linewidth = "none", linetype = "none") +
    labs(y = "Initial numbers", col = "Sex", linetype = "Sex") +
    facet_wrap(~Region) +
    theme_bw()
  return(gplt)
}


#'
#' plot_partition plot numbers at age over time and space
#' @param MLE_report a list that is output from obj$report() usually once an optimsation routine has been done.
#' @param subset_years if you only want to plot the partition for some years
#' @param region_key data.frame with colnames area and TMB_ndx for providing real region names to objects
#' @return ggplot2 object that will plot if an observation occurs in a year and region
#' @export
plot_partition = function(MLE_report, subset_years = NULL, region_key = NULL) {
  full_df = get_partition(MLE_report, region_key);
  if(!is.null(subset_years)) {
    full_df = full_df %>% filter(Year %in% subset_years)
  }

  gplt = ggplot(full_df) +
    geom_line(aes(x = Age, y = Numbers, col = sex, linetype = sex), linewidth= 1.1) +
    guides( linewidth = "none", linetype = "none") +
    labs(y = "Numbers at age (beginning of the year)", col = "Sex", linetype = "Sex") +
    facet_wrap(Year~Region) +
    theme_bw()
  return(gplt)
}


#'
#' calculate_initial_numbers_at_age
#'
#' @param n_regions integer number of regions
#' @param n_ages integer number of ages
#' @param R0 vector of R0 parameters for each region
#' @param movement_matrix movement matrix containing annual movement rates
#' @param natural_mortality vector of natural mortality rates for each age
#' @return matrix of numbers at age with regions being rows and ages being the cols
#' @export
calculate_initial_numbers_at_age <-function(n_regions, n_ages, R0, movement_matrix, natural_mortality) {
  N_age = matrix(0, nrow = n_regions, ncol = n_ages)
  update_N_age = N_age
  for(i in 1:(n_ages)) {
    # recruitment
    update_N_age[,1] = R0 #* exp(-natural_mortality[1])
    # ageing and mortality
    update_N_age[,2:n_ages] = N_age[,1:(n_ages - 1)] * exp(-natural_mortality[1:(n_ages - 1)])
    # plus group
    update_N_age[,n_ages] = update_N_age[,n_ages] + N_age[,n_ages] * exp(-natural_mortality[n_ages])
    # movement
    N_age = t(movement_matrix) %*% update_N_age
  }
  ## calculate one more annual cycle
  update_N_age[,1] = R0
  # ageing and mortality
  update_N_age[,2:n_ages] = N_age[,1:(n_ages - 1)] * exp(-natural_mortality[1:(n_ages - 1)])
  # plus group
  update_N_age[,n_ages] = update_N_age[,n_ages] + N_age[,n_ages] * exp(-natural_mortality[n_ages])
  # movement
  update_N_age = t(movement_matrix) %*% update_N_age
  ## approximate!
  c = update_N_age[,n_ages] / N_age[,n_ages] - 1
  update_N_age[,n_ages] = N_age[,n_ages] * 1 / (1 - c)

  return(update_N_age);
}


#'
#' calculate_initial_numbers_at_age_age_based_movement
#'
#' @param n_regions integer number of regions
#' @param n_ages integer number of ages
#' @param R0 vector of R0 parameters for each region
#' @param old_movement_matrix movement matrix for older fish
#' @param young_movement_matrix movement matrix for young fish
#' @param age_based_movement_ogive age_based movement for young fish
#' @param natural_mortality vector of natural mortality rates for each age
#' @return matrix of numbers at age with regions being rows and ages being the cols
#' @export
calculate_initial_numbers_at_age_age_based_movement <-function(n_regions, n_ages, R0, old_movement_matrix, young_movement_matrix, age_based_movement_ogive, natural_mortality) {
  N_age = matrix(0, nrow = n_regions, ncol = n_ages)
  update_N_age = N_age
  young_N_age = old_N_age = N_age
  for(i in 1:(n_ages)) {
    # recruitment
    update_N_age[,1] = R0 #* exp(-natural_mortality[1])
    # ageing and mortality
    update_N_age[,2:n_ages] = N_age[,1:(n_ages - 1)] * exp(-natural_mortality[1:(n_ages - 1)])
    # plus group
    update_N_age[,n_ages] = update_N_age[,n_ages] + N_age[,n_ages] * exp(-natural_mortality[n_ages])
    ##
    young_N_age = sweep(update_N_age, 2, age_based_movement_ogive, "*")
    old_N_age = sweep(update_N_age, 2, 1 - age_based_movement_ogive, "*")
    # movement
    N_age = t(young_movement_matrix) %*% young_N_age + t(old_movement_matrix) %*% old_N_age
  }
  ## calculate one more annual cycle
  update_N_age[,1] = R0
  # ageing and mortality
  update_N_age[,2:n_ages] = N_age[,1:(n_ages - 1)] * exp(-natural_mortality[1:(n_ages - 1)])
  # plus group
  update_N_age[,n_ages] = update_N_age[,n_ages] + N_age[,n_ages] * exp(-natural_mortality[n_ages])
  # movement
  young_N_age = sweep(update_N_age, 2, age_based_movement_ogive, "*")
  old_N_age = sweep(update_N_age, 2, 1 - age_based_movement_ogive, "*")
  # movement
  N_age = t(young_movement_matrix) %*% young_N_age + t(old_movement_matrix) %*% old_N_age
  ## approximate!
  c = update_N_age[,n_ages] / N_age[,n_ages] - 1
  update_N_age[,n_ages] = N_age[,n_ages] * 1 / (1 - c)

  return(update_N_age);
}
