#' get_tag_release_ndx  utility function to get release event ndx given release region and release year
#'
#' @param region_ndx region index 1:n_regions
#' @param release_event_year_ndx release year index 1:n_years_to_retain_tagged_cohorts_for. 1 indicates tagged released this year, 2 is tags released last year etc.
#' @param n_regions number of regions in the model
#' @return a tag release event index
#' @export
get_tag_release_ndx = function(region_ndx, release_event_year_ndx, n_regions) {
  return ((release_event_year_ndx - 1) * n_regions + (region_ndx - 1) + 1);
}


#'
#' get_tag_release_AF
#' @param data list that is passed to the MakeADfun for the TMB model
#' @param region_key data.frame with colnames area and TMB_ndx for providing real region names to objects
#' @return ggplot2 object that will plot if an observation occurs in a year and region
#' @export
get_tag_release_AF = function(data, region_key = NULL) {
  years = data$years
  ages = data$ages
  regions = 1:data$n_regions
  if(!is.null(region_key))
    regions = region_key$area[match(regions, (region_key$TMB_ndx + 1))]

  ## release event years
  release_years = years[which(data$tag_release_event_this_year == 1)]

  full_df = NULL
  for(y_ndx in 1:length(release_years)) {
    for(r_ndx in 1:length(regions)) {
      male_df = data.frame(N_age = data$male_tagged_cohorts_by_age[, r_ndx, y_ndx], sex = "Male", age = ages, release_region = regions[r_ndx], release_year = release_years[y_ndx])
      female_df = data.frame(N_age = data$female_tagged_cohorts_by_age[, r_ndx, y_ndx], sex = "Female", age = ages, release_region = regions[r_ndx], release_year = release_years[y_ndx])
      full_df = rbind(full_df, male_df, female_df)
    }
  }


  return(full_df)
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

#'
#' get_tag_recovery_obs_fitted_values
#' @param MLE_report a list that is output from obj$report() usually once an optimsation routine has been done.
#' @param region_key data.frame with colnames area and TMB_ndx for providing real region names to objects
#' @return data frame in long format
#' @export
get_tag_recovery_obs_fitted_values = function(MLE_report, region_key = NULL) {
  years = MLE_report$years
  ages = MLE_report$ages
  regions = 1:MLE_report$n_regions
  if(!is.null(region_key))
    regions = region_key$area[match(regions, (region_key$TMB_ndx + 1))]

  ## release event years
  release_years = years[which(MLE_report$tag_release_event_this_year == 1)]
  sex_age_lvls = paste0(rep(c("M","F"), each = length(ages)), ages)
  sex_for_report = rep(c("M","F"), each = length(ages))
  age_for_rep = rep(ages, 2)
  ## recovery years
  recovery_years = years[which(MLE_report$tag_recovery_indicator == 1)]

  full_df = NULL
  for(y_ndx in 1:length(recovery_years)) { ## recovery years
    for(r_ndx in 1:length(regions)) { ## recovery regions
      ## we released tagged fish making this a release event

      ## now link to subsequent release events
      for(release_yr_ndx in 1:length(release_years)) {
        for(release_region_ndx in 1:length(regions)) {
          diff_ = recovery_years[y_ndx] - release_years[release_yr_ndx]
          diff_ = min(c(diff_ + 1, data$n_years_to_retain_tagged_cohorts_for + 1))
          release_event_ndx = get_tag_release_ndx(r_ndx, diff_, data$n_regions)
          if(MLE_report$tag_recovery_indicator_by_release_event_and_recovery_region[release_event_ndx, r_ndx, y_ndx] == 1) {
            ## recovery here
            tmp_df = data.frame(sex = sex_for_report, age = age_for_rep, recovery_year = recovery_years[y_ndx], recovery_region = regions[r_ndx], release_region = regions[release_region_ndx], release_year = release_year, observed = MLE_report$obs_tag_recovery[, release_event_ndx, r_ndx, y_ndx], predicted = MLE_report$pred_tag_recovery[, release_event_ndx, r_ndx, y_ndx])
            full_df = rbind(full_df, tmp_df)
          }
        }
      }
    }
  }
  return(full_df)
}

#'
#' plot_tag_recovery_obs
#' @param MLE_report a list that is output from obj$report() usually once an optimsation routine has been done.
#' @param region_key data.frame with colnames area and TMB_ndx for providing real region names to objects
#' @param release_ndx_to_plot vector of integers to create subset plots
#' @return ggplot2
#' @export
plot_tag_recovery_obs= function(MLE_report, region_key = NULL, release_ndx_to_plot = 1:5, sex = "both") {
  if(!sex %in% c("both", "male","female"))
    stop("sex must be 'both', 'male' or 'female'")
  get_tag_df = get_tag_recovery_obs_fitted_values(MLE_report, region_key)
  get_tag_df$release_event = paste0(get_tag_df$release_year, "-", get_tag_df$release_region)
  unique_release_events = unique(get_tag_df$release_event )
  if(max(release_ndx_to_plot) > length(unique_release_events)) {
    warning(paste0("you are asking to plot up to ", max(release_ndx_to_plot), " release events, but there are only ",  length(unique_release_events), " release events. changing the max value of release_ndx_to_plot"))
    release_ndx_to_plot = release_ndx_to_plot[-which(release_ndx_to_plot >  length(unique_release_events))]
  }
  subset_df = get_tag_df %>% filter(release_event %in% unique_release_events[release_ndx_to_plot])
  if(sex != "both") {
    sex_ndx = "M"
    if(sex == "female")
      sex_ndx = "F"
    subset_df = subset_df %>% filter(sex == sex_ndx)
  }

  if(!is.null(region_key)) {
    subset_df$recovery_region = factor(subset_df$recovery_region, levels = (region_key$area[region_key$TMB_ndx + 1]))
    subset_df$release_region = factor(subset_df$release_region, levels = rev(region_key$area[region_key$TMB_ndx + 1]))
  }

  gplt = ggplot(subset_df, aes(x = age, y = predicted, col = factor(recovery_year), linetype = sex)) +
    geom_line(linewidth = 1.1) +
    facet_grid(recovery_region~release_event) +
    labs(col = "Recovery years")
  return(gplt)
}
