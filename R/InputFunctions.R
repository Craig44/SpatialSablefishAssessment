
#'
#' plot_input_catches
#' @param data list that is passed to the MakeADfun for the TMB model
#' @param region_key data.frame with colnames area and TMB_ndx for providing real region names to objects
#' @return ggplot2 object that will plot if an observation occurs in a year and region
#' @export
plot_input_catches = function(data, region_key = NULL) {
  years = data$years
  regions = 1:data$n_regions
  dimnames(data$fixed_fishery_catch) = dimnames(data$trwl_fishery_catch)  = list(regions, years)
  fixed_catch = reshape2::melt(data$fixed_fishery_catch)
  trwl_catch = reshape2::melt(data$trwl_fishery_catch)
  colnames(fixed_catch) = colnames(trwl_catch) = c("Region", "Year", "Catch")
  fixed_catch$label = "Fixed gear"
  trwl_catch$label = "Trawl"
  full_df = rbind(fixed_catch, trwl_catch)
  if(is.null(region_key)) {
    full_df$Region = paste0("Region ", full_df$Region)
  } else {
    full_df$Region = region_key$area[match(full_df$Region, (region_key$TMB_ndx + 1))]
  }
  gplt = ggplot(full_df) +
    geom_line(aes(x = Year, y = Catch, col = label, linetype = label), linewidth= 1.1) +
    guides( linewidth = "none") +
    labs(y = "Catch", col = "Fishery", linetype = "Fishery") +
    facet_wrap(~Region) +
    theme_bw()
  return(gplt)
}

#'
#' plot_input_observations
#' @param data list that is passed to the MakeADfun for the TMB model
#' @param region_key data.frame with colnames area and TMB_ndx for providing real region names to objects
#' @return ggplot2 object that will plot if an observation occurs in a year and region
#' @export
plot_input_observations = function(data, region_key = NULL) {
  years = data$years
  regions = 1:data$n_regions
  dimnames(data$fixed_catchatage_indicator) = dimnames(data$fixed_catchatlgth_indicator) = dimnames(data$trwl_catchatlgth_indicator) = dimnames(data$srv_dom_ll_catchatage_indicator) = dimnames(data$srv_dom_ll_bio_indicator) = list(regions, years)
  dimnames(data$tag_recovery_indicator_by_release_event_and_recovery_region) = list(1:dim(data$tag_recovery_indicator_by_release_event_and_recovery_region)[1], regions, years[which(data$tag_recovery_indicator == 1)])
  fixed_catchatage = reshape2::melt(data$fixed_catchatage_indicator)
  fixed_catchatlgth  = reshape2::melt(data$fixed_catchatlgth_indicator)
  trwl_catchatlgth = reshape2::melt(data$trwl_catchatlgth_indicator)
  srv_dom_ll_catchatage = reshape2::melt(data$srv_dom_ll_catchatage_indicator)
  srv_dom_ll_bio = reshape2::melt(data$srv_dom_ll_bio_indicator)
  tag_recovery_detailed = NULL
  if(sum(data$tag_recovery_indicator) != 0) {
    tag_recovery_detailed = reshape2::melt(data$tag_recovery_indicator_by_release_event_and_recovery_region)
    colnames(tag_recovery_detailed) = c("Tag release", "Region", "Year", "indicator")
    tag_recovery_detailed$label = "Tag recovery"
    ## collapse tag recoveries across release events
    tag_recovery_detailed = tag_recovery_detailed %>% group_by(Region, Year, label) %>% summarise(indicator = ifelse(sum(indicator)>0, 1, 0))

  }
  colnames(fixed_catchatage) = colnames(fixed_catchatlgth) = colnames(trwl_catchatlgth) = colnames(srv_dom_ll_catchatage) = colnames(srv_dom_ll_bio) = c("Region", "Year", "indicator")
  ## tag releases
  dimnames(data$male_tagged_cohorts_by_age) = dimnames(data$female_tagged_cohorts_by_age) = list(data$ages, regions,  data$years[which(data$tag_release_event_this_year == 1)])
  tag_releases_m = reshape2::melt(data$male_tagged_cohorts_by_age)
  tag_releases_f = reshape2::melt(data$female_tagged_cohorts_by_age)
  tag_releases = rbind(tag_releases_m, tag_releases_f)
  colnames(tag_releases) = c("Age", "Region", "Year", "releases")
  tag_release_df = tag_releases %>% group_by(Region, Year) %>% summarise(indicator = ifelse(sum(releases) > 0, 1, 0))
  tag_release_df$label = "Tag Releases"
  fixed_catchatlgth$label = "Fishery Fixed LF"
  fixed_catchatage$label = "Fishery Fixed AF"
  trwl_catchatlgth$label = "Fishery Trawl LF"
  srv_dom_ll_catchatage$label = "Survey LL AF"
  srv_dom_ll_bio$label = "Survey LL Biomass"
  ## combine
  full_df = rbind(fixed_catchatage, trwl_catchatlgth, srv_dom_ll_catchatage, srv_dom_ll_bio, fixed_catchatlgth, tag_recovery_detailed, tag_release_df)
  if(is.null(region_key)) {
    full_df$Region = paste0("Region ", full_df$Region)
  } else {
    full_df$Region = region_key$area[match(full_df$Region, (region_key$TMB_ndx + 1))]
  }

  full_df$indicator = ifelse(full_df$indicator == 0, NA, 1)
  gplt = ggplot(full_df) +
    geom_point(aes(x = Year, y = label, col = label, size = indicator)) +
    guides(colour = "none", size = "none") +
    labs(y = "") +
    facet_wrap(~Region, ncol = 1) +
    theme_bw()
  return(gplt)
}


#'
#' plot_mean_weight
#' @param data list input for model
#' @return ggplot2
#' @export
plot_mean_weight = function(data) {
  projyears = min(data$years):(max(data$years) + data$n_projections_years)
  dimnames(data$male_mean_weight_by_age) = dimnames(data$female_mean_weight_by_age) = list(data$ages, projyears)

  molten_male_weight_at_age = reshape2::melt(data$male_mean_weight_by_age)
  molten_female_weight_at_age = reshape2::melt(data$female_mean_weight_by_age)

  colnames(molten_female_weight_at_age) = colnames(molten_male_weight_at_age) = c("Age","Year","weight_at_age")
  molten_female_weight_at_age$sex = "Female"
  molten_male_weight_at_age$sex = "Male"
  full_df = rbind(molten_female_weight_at_age, molten_male_weight_at_age)

  gplt = ggplot(full_df, aes(x = Age, y = weight_at_age, col = factor(Year))) +
    geom_line(linewidth = 1.10) +
    labs(x = "Age", y = "Mean weight at age", col = "Year") +
    facet_wrap(~sex) +
    theme_bw()
  return(gplt)
}

#'
#' plot_age_length_matrix
#' @param data list input for model
#' @param subset_years vector of years to subset the age-length plt for
#' @return ggplot2
#' @export
plot_age_length_matrix = function(data, subset_years = NULL) {
  projyears = min(data$years):(max(data$years) + data$n_projections_years)
  dimnames(data$male_age_length_transition) = dimnames(data$female_age_length_transition) = list(data$ages, data$length_bins, projyears)

  molten_male_age_length = reshape2::melt(data$male_age_length_transition)
  molten_female_age_length = reshape2::melt(data$female_age_length_transition)

  colnames(molten_male_age_length) = colnames(molten_female_age_length) = c("Age","Length", "Year","length_at_age")
  molten_female_age_length$sex = "Female"
  molten_male_age_length$sex = "Male"
  full_df = rbind(molten_male_age_length, molten_female_age_length)
  if(!is.null(subset_years))
    full_df = full_df %>% dplyr::filter(Year %in% subset_years)

  gplt = ggplot(full_df, aes(x = Age, y = Length)) +
    geom_tile(aes(fill = length_at_age)) +
    labs(x = "Age", y = "Length", fill = "Proportion") +
    facet_wrap(Year~sex) +
    theme_bw()
  return(gplt)
}
