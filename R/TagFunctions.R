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
  recovery_years = years[which(data$tag_recovery_indicator_by_year == 1)]

  dimnames(data$tag_recovery_indicator) = list(1:dim(data$tag_recovery_indicator)[1], regions, years[which(data$tag_recovery_indicator_by_year == 1)])
  tag_recovery_detailed = reshape2::melt(data$tag_recovery_indicator)
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
            if(data$tag_recovery_indicator[release_event_ndx, possible_recovery_region, recovery_year_ndx] == 1) {
              ## recovery here
              tmp_df = data.frame(release_year = release_years[y_ndx], release_region = regions[r_ndx], recovery_year = tmp_recovery_year, recovery_region = regions[possible_recovery_region], n_releases = n_released_fish, n_recoveries = sum(data$obs_tag_recovery[release_event_ndx, possible_recovery_region, recovery_year_ndx]))
              full_df = rbind(full_df, tmp_df)
            }
          }
        }
      }
    }
  }
  full_df$release_event = paste0(full_df$release_year, "-", full_df$release_region)
  plot_title = ""
  if(length(release_region_to_plt) == 1) {
    ## assume single release region - so move to title instead of facet panel label
    full_df$release_event_label = full_df$release_year
    plot_title = release_region_to_plt
  } else {
    ## multiple release regions
    full_df$release_event_label = full_df$release_event
  }

  unique_release_events = unique(full_df$release_event )
  if(max(release_ndx_to_plot) > length(unique_release_events)) {
    warning(paste0("you are asking to plot up to ", max(release_ndx_to_plot), " release events, but there are only ",  length(unique_release_events), " release events. changing the max value of release_ndx_to_plot"))
    release_ndx_to_plot = release_ndx_to_plot[-which(release_ndx_to_plot >  length(unique_release_events))]
  }
  subset_df = full_df %>% dplyr::filter(release_event %in% unique_release_events[release_ndx_to_plot])
  if(!is.null(region_key)) {
    subset_df$recovery_region = factor(subset_df$recovery_region, levels = rev(region_key$area[region_key$TMB_ndx + 1]))
    subset_df$release_region = factor(subset_df$release_region, levels = rev(region_key$area[region_key$TMB_ndx + 1]))
  }
  subset_df$recovery_year = factor(subset_df$recovery_year, levels = recovery_years)

  subset_df$release_event_with_sample_size = paste0(subset_df$release_event_label, " n: ", subset_df$n_releases)

  if(!is.null(release_region_to_plt))
    subset_df = subset_df %>% dplyr::filter(release_region == release_region_to_plt)

  subset_df$recovery_color = ifelse(subset_df$n_recoveries == 0, "zero", "positive")
  gplt = ggplot(subset_df, aes(x = factor(recovery_year), y = recovery_region)) +
    geom_point(aes(size = n_recoveries, col = recovery_color)) +
    facet_wrap(~release_event_with_sample_size) +
    guides(col = "none") +
    labs(x = "Recovery year", y = "Recovery region", size = "Recoveries", title = plot_title) +
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
  subset_df = tag_df %>% dplyr::filter(release_event %in% unique_release_events[release_ndx_to_plot])
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
get_tag_recovery_obs_fitted_values = function(MLE_report, region_key = NULL, verbose = FALSE) {
  years = MLE_report$years
  ages = MLE_report$ages
  regions = 1:MLE_report$n_regions
  if(!is.null(region_key))
    regions = region_key$area[match(regions, (region_key$TMB_ndx + 1))]

  ## release event years
  release_years = years[which(MLE_report$tag_recovery_indicator_by_year == 1)]
  sex_age_lvls = paste0(rep(c("M","F"), each = length(ages)), ages)
  sex_for_report = rep(c("M","F"), each = length(ages))
  age_for_rep = rep(ages, 2)

  full_df = NULL
  ## tag data is Poisson or Negative binomial structure
  if(MLE_report$tag_likelihood %in% c(0,1)) {
    ## recovery years
    recovery_years = years[which(MLE_report$tag_recovery_indicator_by_year == 1)]

    #dimnames(MLE_report$tag_recovery_indicator_by_year) = list(1:((MLE_report$n_years_to_retain_tagged_cohorts_for + 1) * MLE_report$n_regions), regions, recovery_years)
    #molten_indicator = reshape2::melt(MLE_report$tag_recovery_indicator_by_year)
    #colnames(molten_indicator) = c("release_event", "recovery_region", "recovery_year", "indicator")
    #molten_indicator$unique_recovery_id = paste0(molten_indicator$release_event,"-", molten_indicator$recovery_region, "-", molten_indicator$recovery_year)
    #sum(molten_indicator$indicator) * length(age_for_rep)
    ## to validate data frame for efficient storing during the loop
    expected_n_records = sum(MLE_report$tag_recovery_indicator)

    for(y_ndx in 1:length(recovery_years)) { ## recovery years
      if(verbose)
        cat("recovery year ", recovery_years[y_ndx],"\n")
      for(r_ndx in 1:length(regions)) { ## recovery regions
        ## now link to release events only include release years prior to recovery year
        possible_release_years = release_years[release_years < recovery_years[y_ndx]]
        for(release_yr_ndx in 1:(MLE_report$n_years_to_retain_tagged_cohorts_for + 1)) {
          for(release_region_ndx in 1:length(regions)) {
            release_event_ndx = get_tag_release_ndx(release_region_ndx, release_yr_ndx, MLE_report$n_regions)
            if(MLE_report$tag_recovery_indicator[release_event_ndx, r_ndx, y_ndx] == 1) {
              this_release_year = as.character(recovery_years[y_ndx] - (release_yr_ndx - 1))
              if(release_yr_ndx == (MLE_report$n_years_to_retain_tagged_cohorts_for + 1))
                this_release_year = "Plus group"
              ## recovery here
              tmp_df = data.frame(release_event = release_event_ndx, recovery_year = recovery_years[y_ndx], recovery_region = regions[r_ndx], release_region = regions[release_region_ndx], release_year = this_release_year, observed = MLE_report$obs_tag_recovery[release_event_ndx, r_ndx, y_ndx], predicted = MLE_report$pred_tag_recovery[release_event_ndx, r_ndx, y_ndx])
              full_df = rbind(full_df, tmp_df)
            }
          }
        }
      }
    }
    full_df$unique_recovery_id = paste0(full_df$release_event,"-", full_df$recovery_region, "-", full_df$recovery_year)
    if(MLE_report$tag_likelihood == 0) {
      full_df$likelihood = "Poisson"
    } else if (MLE_report$tag_likelihood == 0) {
      full_df$likelihood = "Negative Binomial"
    }
  } else {
    expected_n_records = sum(MLE_report$tag_recovery_indicator) * (MLE_report$n_years_to_retain_tagged_cohorts_for * MLE_report$n_regions + 1)
    ## tag data is Multinomial release conditioned
    recovery_regions = rep(regions, MLE_report$n_years_to_retain_tagged_cohorts_for)
    recovery_regions = c(recovery_regions, "Not-recaptured")
    recovery_year_ndx = rep(1:MLE_report$n_years_to_retain_tagged_cohorts_for, each = length(regions)) - 1
    for(y_ndx in 1:length(years)) { ## recovery years
      for(r_ndx in 1:length(regions)) { ## recovery regions
        if(MLE_report$tag_recovery_indicator[y_ndx, r_ndx] == 1) {
          recovery_years = years[y_ndx] + recovery_year_ndx
          recovery_years = c(recovery_years, "Not-recaptured")

          tmp_df = data.frame(release_year = years[y_ndx], release_region = regions[r_ndx], recovery_year = recovery_years, recovery_region = recovery_regions, observed = MLE_report$obs_tag_recovery[, r_ndx, y_ndx], predicted = MLE_report$pred_tag_recovery[, r_ndx, y_ndx])
          full_df = rbind(full_df, tmp_df)

        }
      }
    }
    full_df$likelihood = "Multinomial-release"
  }
  #nrow(full_df)
  #no_derived_values = molten_indicator$unique_recovery_id[molten_indicator$indicator == 1][which(!molten_indicator$unique_recovery_id[molten_indicator$indicator == 1] %in% unique(full_df$unique_recovery_id))]

  if(expected_n_records != nrow(full_df))
    stop(paste0("expected ", expected_n_records , " unique recovery observations. But could only find ", nrow(full_df)))

  full_df$release_event = paste0(full_df$release_year, "-", full_df$release_region)
  return(full_df)
}

#'
#' plot_tag_recovery_obs plot predicted and observed tag-recovery obs
#' @param MLE_report a list that is output from obj$report() usually once an optimsation routine has been done.
#' @param region_key data.frame with colnames area and TMB_ndx for providing real region names to objects
#' @param release_ndx_to_plot vector of integers to create subset plots
#' @return ggplot2
#' @export
plot_tag_recovery_obs = function(MLE_report, region_key = NULL, release_ndx_to_plot = 1:5, sex = "both") {
  if(!sex %in% c("both", "male","female"))
    stop("sex must be 'both', 'male' or 'female'")
  if(!MLE_report$tag_likelihood %in% c(0,1))
    stop("This plot function doesn't work for this tag-likelihood")

  get_tag_df = get_tag_recovery_obs_fitted_values(MLE_report, region_key)
  unique_release_events = unique(get_tag_df$release_event )
  if(max(release_ndx_to_plot) > length(unique_release_events)) {
    warning(paste0("you are asking to plot up to ", max(release_ndx_to_plot), " release events, but there are only ",  length(unique_release_events), " release events. changing the max value of release_ndx_to_plot"))
    release_ndx_to_plot = release_ndx_to_plot[-which(release_ndx_to_plot >  length(unique_release_events))]
  }
  subset_df = get_tag_df %>% dplyr::filter(release_event %in% unique_release_events[release_ndx_to_plot])
  if(sex != "both") {
    sex_ndx = "M"
    if(sex == "female")
      sex_ndx = "F"
    subset_df = subset_df %>% dplyr::filter(sex == sex_ndx)
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


#'
#' plot_tag_recovery_fits plot predicted and observed tag-recovery obs
#' @param MLE_report a list that is output from obj$report() usually once an optimsation routine has been done.
#' @param region_key data.frame with colnames area and TMB_ndx for providing real region names to objects
#' @param plt_type string for the plot type
#' \itemize{
#'   \item `aggregate`: aggregated observed and predicted
#'   \item `stand_resid`: boxplots of standardised residuals for release region and recovery region i.e., boxplots are all the release-events
#'   \item `rqr_resid`: boxplots of randomised quantile residuals for release region and recovery region i.e., boxplots are all the release-events
#'   \item `raw_resid`: points and lines of raw residuals by recovery year (aggregated over release events) facets are release region by recovery region i.e., can by multiple release years for a recovery year
#' }
#' @return ggplot2 facets are columns recovery regions rows are release regions
#' @export
plot_tag_recovery_fits <- function(MLE_report, region_key = NULL, plt_type = "aggregate") {
  if(!plt_type %in% c("aggregate", "stand_resid", "rqr_resid", "raw_resid"))
    stop('plt_type needs to be one of "aggregate", "stand_resid", "rqr_resid", "raw_resid"')
  if(!MLE_report$tag_likelihood %in% c(0,1))
    stop("This plot function doesn't work for this tag-likelihood")

  tag_fits_by_age_sex = get_tag_recovery_obs_fitted_values(MLE_report, region_key)
  regions = paste0("Region ", 1:MLE_report$n_regions)
  if(!is.null(region_key))
    regions = region_key$area[region_key$TMB_ndx + 1]

  if(!is.null(region_key)) {
    tag_fits_by_age_sex$recovery_region = factor(tag_fits_by_age_sex$recovery_region, levels = rev(regions), ordered = T)
    tag_fits_by_age_sex$release_region = factor(tag_fits_by_age_sex$release_region, levels = rev(regions), ordered = T)
  }
  gplt = NULL
  ## aggregate to release region and recovery region
  if(plt_type == "aggregate") {
    tag_aggregated = tag_fits_by_age_sex %>% group_by(release_region, recovery_region) %>% summarise(observed = sum(observed), predicted = sum(predicted))
    tag_aggregated$resid = tag_aggregated$observed - tag_aggregated$predicted
    tag_aggregated$resid_sign = ifelse(tag_aggregated$resid < 0, "negative", "positive")

    gplt = ggplot(tag_aggregated, aes(x = release_region, y = recovery_region, fill = resid)) +
      geom_tile() +
      scale_fill_gradient(low = "red", high = "blue") +
      geom_text(aes(x = recovery_region, y = release_region, label = round(resid,2)), color = "black", size = 4) +
      labs(x = "Recovery", y = "Release", fill = "Aggregated\nresiduals")

  } else {
    tag_fits = tag_fits_by_age_sex %>% group_by(release_event, unique_recovery_id) %>% summarise(observed = sum(observed), predicted = sum(predicted), recovery_year = unique(recovery_year), recovery_region = unique(recovery_region), release_region = unique(release_region), release_year = unique(release_year))
    ## residuals
    tag_fits$resid = tag_fits$observed - tag_fits$predicted
    tag_fits$resid_sign = ifelse(tag_fits$resid < 0, "negative", "positive")
    if(MLE_report$tag_likelihood == 0) {
      ## Poisson
      tag_fits$rqr = pnorm(ppois(tag_fits$observed, tag_fits$predicted))
      tag_fits$stand_resid = tag_fits$resid / tag_fits$predicted
    } else if(MLE_report$tag_likelihood == 1) {
      ## Negative Binomial
      tag_fits$rqr = pnorm(pnbinom(tag_fits$observed, size = MLE_report$tag_phi, mu = tag_fits$predicted))
      tag_var = tag_fits$predicted + (tag_fits$predicted^2/MLE_report$tag_phi)
      tag_fits$stand_resid = tag_fits$resid / sqrt(tag_var)
    }

    if(plt_type == "stand_resid") {
      gplt = ggplot(data = tag_fits) +
        geom_boxplot(aes(y = stand_resid)) +
        ylim(-10,10) +
        geom_hline(yintercept = 0, col = "red", linewidth = 1.1, linetype = "dashed") +
        theme_bw()+
        facet_grid(recovery_region~release_region)
    }
    if(plt_type == "stand_resid") {
      gplt = ggplot(data = tag_fits) +
        geom_boxplot(aes(y = rqr)) +
        theme_bw()+
        facet_grid(recovery_region~release_region)
    }
    if(plt_type == "raw_resid") {
      gplt = ggplot(data = tag_fits) +
        geom_line(aes(x = recovery_year, y = resid)) +
        geom_point(aes(x = recovery_year, y = resid)) +
        theme_bw()+
        labs(x = "recovery year",y = "Raw residual") +
        facet_grid(recovery_region~release_region)
    }
  }

  return(gplt)
}

#'
#' get_tag_reporting_rate get a data frame of reporting rate parameters
#' @param MLE_report a list that is output from obj$report() usually once an optimsation routine has been done.
#' @param region_key data.frame with colnames area and TMB_ndx for providing real region names to objects
#' @export
#' @return data.frame of tag-reporting rates by tag-recovery year recovery region
get_tag_reporting_rate <- function(MLE_report, region_key = NULL) {
  region_label = paste0("Region", 1:MLE_report$n_regions)
  if(!is.null(region_key))
    region_label = region_key$area[region_key$TMB_ndx + 1]
  recovery_years = MLE_report$years[which(MLE_report$tag_recovery_indicator_by_year == 1)]
  dimnames(MLE_report$tag_reporting_rate) = list(region_label, recovery_years)
  molten_report_rates = reshape2::melt(MLE_report$tag_reporting_rate)
  colnames(molten_report_rates) = c("Region", "Year", "ReportingRate")
  ## no recovery events so no reporting rates
  if(length(recovery_years) == 0)
    return(NULL)
  return(molten_report_rates)
}

#' plot_tag_reporting_rate create a ggplot of reporring rate by recovery year and region
#' @param MLE_report a list that is output from obj$report() usually once an optimsation routine has been done.
#' @param region_key data.frame with colnames area and TMB_ndx for providing real region names to objects
#' @export
#' @return ggplot of reporting rates
plot_tag_reporting_rate <- function(MLE_report, region_key = NULL) {
  report_rate_df = get_tag_reporting_rate(MLE_report, region_key)

  gplt = ggplot(data = report_rate_df, aes(x = Year, y = ReportingRate, col = Region, linetype = Region)) +
    geom_line(linewidth = 1.1) +
    labs(x = "Tag recovery year", y = "Tag reporting rate", col = "Region", linetype = "Region") +
    ylim(0, NA) +
    theme_classic() +
    theme(
          axis.text = element_text(size = 14),
          axis.title = element_text(size = 14),
          strip.text = element_text(size=14),
          legend.text = element_text(size=14))
  return(gplt)
}
