
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
  tag_recovery_detailed = reshape2::melt(data$tag_recovery_indicator_by_release_event_and_recovery_region)
  colnames(fixed_catchatage) = colnames(fixed_catchatlgth) = colnames(trwl_catchatlgth) = colnames(srv_dom_ll_catchatage) = colnames(srv_dom_ll_bio) = c("Region", "Year", "indicator")
  colnames(tag_recovery_detailed) = c("Tag release", "Region", "Year", "indicator")
  fixed_catchatlgth$label = "Fishery Fixed LF"
  fixed_catchatage$label = "Fishery Fixed AF"
  trwl_catchatlgth$label = "Fishery Trawl LF"
  srv_dom_ll_catchatage$label = "Survey LL AF"
  srv_dom_ll_bio$label = "Survey LL Biomass"
  tag_recovery_detailed$label = "Tag recovery"
  ## collapse tag recoveries across release events
  tag_recovery_detailed = tag_recovery_detailed %>% group_by(Region, Year, label) %>% summarise(indicator = ifelse(sum(indicator)>0, 1, 0))
  ## combine
  full_df = rbind(fixed_catchatage, trwl_catchatlgth, srv_dom_ll_catchatage, srv_dom_ll_bio, fixed_catchatlgth, tag_recovery_detailed)
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
#' plot_tag_release_and_recoveries
#' @param data list that is passed to the MakeADfun for the TMB model
#' @param region_key data.frame with colnames area and TMB_ndx for providing real region names to objects
#' @param release_ndx_to_plot vector of integers to create subset plots
#' @param release_region_to_plt string matching release region to subset the plot for
#' @return ggplot2 object that will plot if an observation occurs in a year and region
#' @export
plot_frequency_of_tag_release_and_recoveries = function(data, region_key = NULL, release_ndx_to_plot = 1:10, release_region_to_plt = NULL) {
  years = data$years
  ages = data$ages
  regions = 1:data$n_regions
  if(!is.null(region_key))
    regions = region_key$area[match(regions, (region_key$TMB_ndx + 1))]

  ## release event years
  release_years = years[which(data$tag_release_event_this_year == 1)]
  ## recovery years
  recovery_years = years[which(data$tag_recovery_indicator == 1)]

  dimnames(data$tag_recovery_indicator_by_release_event_and_recovery_region) = list(1:dim(data$tag_recovery_indicator_by_release_event_and_recovery_region)[1], regions, years[which(data$tag_recovery_indicator == 1)])
  tag_recovery_detailed = reshape2::melt(data$tag_recovery_indicator_by_release_event_and_recovery_region)
  full_df = NULL
  for(y_ndx in 1:length(release_years)) {
    for(r_ndx in 1:length(regions)) {
      n_released_fish = (sum(data$male_tagged_cohorts_by_age[, r_ndx, y_ndx]) + sum(data$female_tagged_cohorts_by_age[, r_ndx, y_ndx]))
      if(n_released_fish > 0) {
        ## we released tagged fish making this a release event
        ## now link to subsequent recoveries
        for(possible_recovery_years in 1:(data$n_years_to_retain_tagged_cohorts_for)) {
          for(possible_recovery_region in 1:length(regions)) {
            tmp_recovery_year = release_years[y_ndx] + possible_recovery_years
            if(!tmp_recovery_year %in% recovery_years)
              next;
            recovery_year_ndx = which(recovery_years %in% tmp_recovery_year)
            release_event_ndx = get_tag_release_ndx(r_ndx, possible_recovery_years + 1, data$n_regions)
            if(data$tag_recovery_indicator_by_release_event_and_recovery_region[release_event_ndx, possible_recovery_region, recovery_year_ndx] == 1) {
              ## recovery here
              tmp_df = data.frame(release_year = release_years[y_ndx], release_region = regions[r_ndx], recovery_year = tmp_recovery_year, recovery_region = regions[possible_recovery_region], n_releases = n_released_fish, n_recoveries = sum(data$obs_tag_recovery[, release_event_ndx, possible_recovery_region, recovery_year_ndx]))
              full_df = rbind(full_df, tmp_df)
            }
          }
        }
      }
    }
  }
  full_df$release_event = paste0(full_df$release_year, "-", full_df$release_region)
  unique_release_events = unique(full_df$release_event )
  if(max(release_ndx_to_plot) > length(unique_release_events)) {
    warning(paste0("you are asking to plot up to ", max(release_ndx_to_plot), " release events, but there are only ",  length(unique_release_events), " release events. changing the max value of release_ndx_to_plot"))
    release_ndx_to_plot = release_ndx_to_plot[-which(release_ndx_to_plot >  length(unique_release_events))]
  }
  subset_df = full_df %>% filter(release_event %in% unique_release_events[release_ndx_to_plot])
  if(!is.null(region_key)) {
    subset_df$recovery_region = factor(subset_df$recovery_region, levels = rev(region_key$area[region_key$TMB_ndx + 1]))
    subset_df$release_region = factor(subset_df$release_region, levels = rev(region_key$area[region_key$TMB_ndx + 1]))
  }
  subset_df$recovery_year = factor(subset_df$recovery_year, levels = recovery_years)
  subset_df$release_event_with_sample_size = paste0(subset_df$release_event, " releases: ", subset_df$n_releases)

  if(!is.null(release_region_to_plt))
    subset_df = subset_df %>% filter(release_region == release_region_to_plt)

  gplt = ggplot(subset_df, aes(x = factor(recovery_year), y = recovery_region)) +
    geom_point(aes(size = n_recoveries)) +
    facet_wrap(~release_event_with_sample_size) +
    labs(x = "Recovery year", y = "Recovery region", size = "Recoveries") +
    theme_bw() +
    scale_size_area() +
    scale_x_discrete(breaks = every_nth(n = 5))

  return(gplt)
}
#'
#' plot_tag_release_AF
#' @param data list that is passed to the MakeADfun for the TMB model
#' @param region_key data.frame with colnames area and TMB_ndx for providing real region names to objects
#' @param release_ndx_to_plot vector of integers to create subset plots
#' @return ggplot2 object that will plot if an observation occurs in a year and region
#' @export
plot_tag_release_AF = function(data, region_key = NULL, release_ndx_to_plot = 1:10) {
  tag_df = get_tag_release_AF(data, region_key)

  tag_df$release_event = paste0(tag_df$release_year, "-", tag_df$release_region)

  unique_release_events = unique(tag_df$release_event )
  if(max(release_ndx_to_plot) > length(unique_release_events)) {
    warning(paste0("you are asking to plot up to ", max(release_ndx_to_plot), " release events, but there are only ",  length(unique_release_events), " release events. changing the max value of release_ndx_to_plot"))
    release_ndx_to_plot = release_ndx_to_plot[-which(release_ndx_to_plot >  length(unique_release_events))]
  }
  subset_df = tag_df %>% filter(release_event %in% unique_release_events[release_ndx_to_plot])
  if(!is.null(region_key)) {
    subset_df$release_region = factor(subset_df$release_region, levels = rev(region_key$area[region_key$TMB_ndx + 1]))
  }

  gplt = ggplot(subset_df, aes(x = age, y = N_age, col = sex, linetype =sex)) +
    geom_line(linewidth = 1.1) +
    facet_wrap(~release_event) +
    labs(x = "Age", y = "AF", col = "Sex", linetype = "Sex") +
    theme_bw()
  return(gplt)
}

