#
# Partition functions
#
#

#'
#' get_SSB
#' @param MLE_report a list that is output from obj$report() usually once an optimsation routine has been done.
#' @param region_key data.frame with colnames area and TMB_ndx for providing real region names to objects
#' @return ggplot2 object that will plot if an observation occurs in a year and region
#' @export
get_SSB = function(MLE_report, region_key = NULL) {
  years = MLE_report$years
  regions = 1:MLE_report$n_regions
  ssbs = MLE_report$SSB_yr
  dimnames(ssbs) = list(years, regions)
  molten_ssbs = reshape2::melt(ssbs)
  colnames(molten_ssbs) = c("Year", "Region", "SSB")
  if(is.null(region_key)) {
    molten_ssbs$Region = paste0("Region ", molten_ssbs$Region)
  } else {
    molten_ssbs$Region = region_key$area[match(molten_ssbs$Region, (region_key$TMB_ndx + 1))]
  }
  return(molten_ssbs);
}

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
#' plot_SSB
#' @param MLE_report a list that is output from obj$report() usually once an optimsation routine has been done.
#' @param data list that is passed to the MakeADfun for the TMB model
#' @param region_key data.frame with colnames area and TMB_ndx for providing real region names to objects
#' @return ggplot2 object that will plot if an observation occurs in a year and region
#' @export
plot_SSB = function(MLE_report, region_key = NULL) {
  molten_ssbs = get_SSB(MLE_report, region_key)
  gplt = ggplot(molten_ssbs) +
    geom_line(aes(x = Year, y = SSB, col = Region), linewidth= 1.1) +
    guides(colour = "none", linewidth = "none") +
    labs(y = "Spawning Stock Biomass (SSB)") +
    facet_wrap(~Region) +
    theme_bw()
  return(gplt)
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


