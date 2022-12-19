#
#
#'
#' plot_AF
#' @param MLE_report a list that is output from obj$report() usually once an optimsation routine has been done.
#' @param label character labeling the observation you want to plot. See below for options
#' \itemize{
#'   \item trwl
#'   \item fixed
#'   \item srv_dom_ll
#' }
#' @param subset_years vector of years to plot it for
#' @param sex character that allows users to specify if the want sex specific plots
#' @param region_key data.frame with colnames area and TMB_ndx for providing real region names to objects
#' @return ggplot2 object that will plot if an observation occurs in a year and region
#' @export
plot_AF = function(MLE_report, label = "fixed", subset_years = NULL, sex = "both", region_key = NULL) {
  if(!label %in% c("fixed","srv_dom_ll"))
    stop("label not one of the expected values.")
  if(!sex %in% c("both", "male", "female"))
    stop('label not one of the expected values. Expected one of the following "both", "male", "female"')

  years = MLE_report$years
  regions = 1:MLE_report$n_regions
  ages = MLE_report$ages
  ## get objects
  obs_indicator = get(paste0(label,"_catchatage_indicator"), MLE_report)
  obs_df = get(paste0("obs_",label,"_catchatage"), MLE_report)
  pred_df = get(paste0("pred_",label,"_catchatage"), MLE_report)
  dimnames(obs_df) = dimnames(pred_df) = list(c(paste0("M_",ages), paste0("F_",ages)), regions, years)
  dimnames(obs_indicator) = list(regions, years)
  NA_ndx = which(obs_indicator == 0, arr.ind = T)
  if(nrow(NA_ndx) > 0) {
    for(i in 1:nrow(NA_ndx)) {
      obs_df[,NA_ndx[i,1], NA_ndx[i,2]] = NA
      pred_df[,NA_ndx[i,1], NA_ndx[i,2]] = NA
    }
  }
  molten_obs = reshape2::melt(obs_df)
  molten_pred = reshape2::melt(pred_df)
  colnames(molten_obs) = c("S_Age", "Region", "Year", "Observed")
  molten_obs$Predicted = molten_pred$value
  if(is.null(region_key)) {
    molten_obs$Region = paste0("Region ", molten_obs$Region)
  } else {
    molten_obs$Region = region_key$area[match(molten_obs$Region, (region_key$TMB_ndx + 1))]
  }

  molten_obs$Age = as.numeric(substring(molten_obs$S_Age, first = 3))
  molten_obs$Sex = ifelse(substring(molten_obs$S_Age, first = 0, last = 1) == "M", "Male", "Female")
  full_df = molten_obs
  ## multiple predicted proportions by effective sample size
  full_df= full_df %>% group_by(Year, Region) %>% mutate(Predicted = Predicted * sum(Observed))

  if(!is.null(subset_years)) {
    full_df = full_df %>% filter(Year %in% subset_years)
  }
  if(sex == "male")
    full_df = full_df %>% filter(Sex == "Male")
  if(sex == "female")
    full_df = full_df %>% filter(Sex == "Female")

  ## plot
  gplt = ggplot(full_df, aes(x = Age)) +
    geom_point(aes(y = Observed, col = "Observed", shape = Sex, group = Sex)) +
    geom_line(aes(y = Predicted, col = "Predicted", linetype = Sex, group = Sex), linewidth= 1.1) +
    guides( linewidth = "none") +
    labs(y = "AF", col = "", linetype = "") +
    facet_grid(Year ~ Region) +
    theme_bw()
  return(gplt)
}
#'
#' plot_mean_age
#' @param plot_mean_age a list that is output from obj$report() usually once an optimsation routine has been done.
#' @param label character labeling the observation you want to plot. See below for options
#' \itemize{
#'   \item all
#'   \item trwl
#'   \item fixed
#'   \item srv_dom_ll
#' }
#' @param subset_years vector of years to plot it for
#' @param sex character that allows users to specify if the want sex specific plots
#' @param region_key data.frame with colnames area and TMB_ndx for providing real region names to objects
#' @return ggplot2 object that will plot if an observation occurs in a year and region
#' @export
plot_mean_age = function(MLE_report, label = "fixed", subset_years = NULL, sex = "both", region_key = NULL) {
  if(!label %in% c("all","fixed","srv_dom_ll"))
    stop("label not one of the expected values.")
  if(!sex %in% c("both", "male", "female"))
    stop('label not one of the expected values. Expected one of the following "both", "male", "female"')

  years = MLE_report$years
  regions = 1:MLE_report$n_regions
  ages = MLE_report$ages
  ## get objects
  if(label != "all") {
    obs_indicator = get(paste0(label,"_catchatage_indicator"), MLE_report)
    obs_df = get(paste0("obs_",label,"_catchatage"), MLE_report)
    pred_df = get(paste0("pred_",label,"_catchatage"), MLE_report)
    dimnames(obs_df) = dimnames(pred_df) = list(c(paste0("M_",ages), paste0("F_",ages)), regions, years)
    dimnames(obs_indicator) = list(regions, years)
    NA_ndx = which(obs_indicator == 0, arr.ind = T)
    if(nrow(NA_ndx) > 0) {
      for(i in 1:nrow(NA_ndx)) {
        obs_df[,NA_ndx[i,1], NA_ndx[i,2]] = NA
        pred_df[,NA_ndx[i,1], NA_ndx[i,2]] = NA
      }
    }
    molten_obs = reshape2::melt(obs_df)
    molten_pred = reshape2::melt(pred_df)
    colnames(molten_obs) = c("S_Age", "Region", "Year", "Observed")
    molten_obs$Predicted = molten_pred$value
    if(is.null(region_key)) {
      molten_obs$Region = paste0("Region ", molten_obs$Region)
    } else {
      molten_obs$Region = region_key$area[match(molten_obs$Region, (region_key$TMB_ndx + 1))]
    }

    molten_obs$Age = as.numeric(substring(molten_obs$S_Age, first = 3))
    molten_obs$Sex = ifelse(substring(molten_obs$S_Age, first = 0, last = 1) == "M", "Male", "Female")
    molten_obs$label = label

    full_df = molten_obs
  } else {
    full_df = NULL;
    obs_labs = c("trwl","fixed","srv_dom_ll")
    for(i in 1:length(obs_labs)) {
      obs_indicator = get(paste0(obs_labs[i],"_catchatage_indicator"), MLE_report)
      obs_df = get(paste0("obs_",obs_labs[i],"_catchatage"), MLE_report)
      pred_df = get(paste0("pred_",obs_labs[i],"_catchatage"), MLE_report)
      dimnames(obs_df) = dimnames(pred_df) = list(c(paste0("M_",ages), paste0("F_",ages)), regions, years)
      dimnames(obs_indicator) = list(regions, years)
      NA_ndx = which(obs_indicator == 0, arr.ind = T)
      if(nrow(NA_ndx) > 0) {
        for(i in 1:nrow(NA_ndx)) {
          obs_df[,NA_ndx[i,1], NA_ndx[i,2]] = NA
          pred_df[,NA_ndx[i,1], NA_ndx[i,2]] = NA
        }
      }
      molten_obs = reshape2::melt(obs_df)
      molten_pred = reshape2::melt(pred_df)
      colnames(molten_obs) = c("S_Age", "Region", "Year", "Observed")
      molten_obs$Predicted = molten_pred$value
      if(is.null(region_key)) {
        molten_obs$Region = paste0("Region ", molten_obs$Region)
      } else {
        molten_obs$Region = region_key$area[match(molten_obs$Region, (region_key$TMB_ndx + 1))]
      }
      molten_obs$Age = as.numeric(substring(molten_obs$S_Age, first = 3))
      molten_obs$Sex = ifelse(substring(molten_obs$S_Age, first = 0, last = 1) == "M", "Male", "Female")
      molten_obs$label = obs_labs[i]
      full_df = rbind(full_df, molten_obs)
    }
  }

  ## multiple predicted proportions by effective sample size
  full_df= full_df %>% group_by(Year, Region, label) %>% mutate(N_eff = sum(Observed), Observed_mean_prop = Observed / N_eff)

  full_df= full_df %>% group_by(Year, Region, label, Sex) %>% summarise(Ey = sum(Age * Predicted), Oy = sum(Age * Observed_mean_prop), E_squared_y = sum(Age^2 * Predicted), N_eff = mean(N_eff))
  full_df$Ry = full_df$Oy - full_df$Ey
  full_df$SEy = sqrt((full_df$E_squared_y - full_df$Ey^2) / full_df$N_eff)
  full_df$'Std.res' <- (full_df$Oy - full_df$Ey)/full_df$SEy
  ## I think this is the final Francis weighting value TODO: to check
  Nmult <- 1 / var(full_df$'Std.res',na.rm=TRUE)
  # Find the adjusted confidence intervals
  full_df$ObsloAdj <- full_df$Oy - 2 * full_df$SEy / sqrt(Nmult)
  full_df$ObshiAdj <- full_df$Oy + 2 * full_df$SEy / sqrt(Nmult)

  if(!is.null(subset_years)) {
    full_df = full_df %>% filter(Year %in% subset_years)
  }
  if(sex == "male")
    full_df = full_df %>% filter(Sex == "Male")
  if(sex == "female")
    full_df = full_df %>% filter(Sex == "Female")

  ## plot
  gplt = ggplot(full_df, aes(x = Year)) +
    geom_point(aes(y = Oy, col = "Observed", shape = Sex, group = Sex)) +
    geom_line(aes(y = Ey, col = "Predicted", linetype = Sex, group = Sex), linewidth= 1.1) +
    geom_errorbar(aes(ymin=ObsloAdj, ymax=ObshiAdj, col = "Observed"), width=.2, position=position_dodge(.9)) +
    guides( linewidth = "none") +
    labs(y = "Mean age", col = "", linetype = "") +
    facet_grid(label ~ Region) +
    theme_bw()
  return(gplt)
}

#'
#' plot_LF
#' @param MLE_report a list that is output from obj$report() usually once an optimsation routine has been done.
#' @param label character labeling the observation you want to plot. See below for options
#' \itemize{
#'   \item trwl
#'   \item fixed
#' }
#' @param subset_years vector of years to plot it for
#' @param sex character that allows users to specify if the want sex specific plots
#' @param region_key data.frame with colnames area and TMB_ndx for providing real region names to objects
#' @return ggplot2 object that will plot if an observation occurs in a year and region
#' @export
plot_LF = function(MLE_report, label = "fixed", subset_years = NULL, sex = "both", region_key = NULL) {
  if(!label %in% c("fixed","trwl"))
    stop("label not one of the expected values.")
  if(!sex %in% c("both", "male", "female"))
    stop('label not one of the expected values. Expected one of the following "both", "male", "female"')

  years = MLE_report$years
  regions = 1:MLE_report$n_regions
  length_bins = MLE_report$length_bins
  ## get objects
  obs_indicator = get(paste0(label,"_catchatlgth_indicator"), MLE_report)
  obs_df = get(paste0("obs_",label,"_catchatlgth"), MLE_report)
  pred_df = get(paste0("pred_",label,"_catchatlgth"), MLE_report)
  dimnames(obs_df) = dimnames(pred_df) = list(c(paste0("M_",length_bins), paste0("F_",length_bins)), regions, years)
  dimnames(obs_indicator) = list(regions, years)
  NA_ndx = which(obs_indicator == 0, arr.ind = T)
  if(nrow(NA_ndx) > 0) {
    for(i in 1:nrow(NA_ndx)) {
      obs_df[,NA_ndx[i,1], NA_ndx[i,2]] = NA
      pred_df[,NA_ndx[i,1], NA_ndx[i,2]] = NA
    }
  }
  molten_obs = reshape2::melt(obs_df)
  molten_pred = reshape2::melt(pred_df)
  colnames(molten_obs) = c("S_Length", "Region", "Year", "Observed")
  molten_obs$Predicted = molten_pred$value
  if(is.null(region_key)) {
    molten_obs$Region = paste0("Region ", molten_obs$Region)
  } else {
    molten_obs$Region = region_key$area[match(molten_obs$Region, (region_key$TMB_ndx + 1))]
  }

  molten_obs$Length = as.numeric(substring(molten_obs$S_Length, first = 3))
  molten_obs$Sex = ifelse(substring(molten_obs$S_Length, first = 0, last = 1) == "M", "Male", "Female")
  full_df = molten_obs
  ## multiple predicted proportions by effective sample size
  full_df= full_df %>% group_by(Year, Region) %>% mutate(Predicted = Predicted * sum(Observed))

  if(!is.null(subset_years)) {
    full_df = full_df %>% filter(Year %in% subset_years)
  }
  if(sex == "male")
    full_df = full_df %>% filter(Sex == "Male")
  if(sex == "female")
    full_df = full_df %>% filter(Sex == "Female")

  ## plot
  gplt = ggplot(full_df, aes(x = Length)) +
    geom_point(aes(y = Observed, col = "Observed", shape = Sex, group = Sex)) +
    geom_line(aes(y = Predicted, col = "Predicted", linetype = Sex, group = Sex), linewidth= 1.1) +
    guides( linewidth = "none") +
    labs(y = "LF", col = "", linetype = "") +
    facet_grid(Year ~ Region) +
    theme_bw()
  return(gplt)
}
#'
#' plot_mean_length
#' @param plot_mean_length a list that is output from obj$report() usually once an optimsation routine has been done.
#' @param label character labeling the observation you want to plot. See below for options
#' \itemize{
#'   \item all
#'   \item trwl
#'   \item fixed
#' }
#' @param subset_years vector of years to plot it for
#' @param sex character that allows users to specify if the want sex specific plots
#' @param region_key data.frame with colnames area and TMB_ndx for providing real region names to objects
#' @return ggplot2 object that will plot if an observation occurs in a year and region
#' @export
plot_mean_length = function(MLE_report, label = "fixed", subset_years = NULL, sex = "both", region_key = NULL) {
  if(!label %in% c("all","fixed","trwl"))
    stop("label not one of the expected values.")
  if(!sex %in% c("both", "male", "female"))
    stop('label not one of the expected values. Expected one of the following "both", "male", "female"')

  years = MLE_report$years
  regions = 1:MLE_report$n_regions
  length_bins = MLE_report$length_bins
  ## get objects
  if(label != "all") {
    obs_indicator = get(paste0(label,"_catchatlgth_indicator"), MLE_report)
    obs_df = get(paste0("obs_",label,"_catchatlgth"), MLE_report)
    pred_df = get(paste0("pred_",label,"_catchatlgth"), MLE_report)
    dimnames(obs_df) = dimnames(pred_df) = list(c(paste0("M_",length_bins), paste0("F_",length_bins)), regions, years)
    dimnames(obs_indicator) = list(regions, years)
    NA_ndx = which(obs_indicator == 0, arr.ind = T)
    if(nrow(NA_ndx) > 0) {
      for(i in 1:nrow(NA_ndx)) {
        obs_df[,NA_ndx[i,1], NA_ndx[i,2]] = NA
        pred_df[,NA_ndx[i,1], NA_ndx[i,2]] = NA
      }
    }
    molten_obs = reshape2::melt(obs_df)
    molten_pred = reshape2::melt(pred_df)
    colnames(molten_obs) = c("S_Length", "Region", "Year", "Observed")
    molten_obs$Predicted = molten_pred$value
    if(is.null(region_key)) {
      molten_obs$Region = paste0("Region ", molten_obs$Region)
    } else {
      molten_obs$Region = region_key$area[match(molten_obs$Region, (region_key$TMB_ndx + 1))]
    }

    molten_obs$Length = as.numeric(substring(molten_obs$S_Length, first = 3))
    molten_obs$Sex = ifelse(substring(molten_obs$S_Length, first = 0, last = 1) == "M", "Male", "Female")
    molten_obs$label = label

    full_df = molten_obs
  } else {
    full_df = NULL;
    obs_labs = c("trwl","fixed")
    for(i in 1:length(obs_labs)) {
      obs_indicator = get(paste0(obs_labs[i],"_catchatlgth_indicator"), MLE_report)
      obs_df = get(paste0("obs_",obs_labs[i],"_catchatlgth"), MLE_report)
      pred_df = get(paste0("pred_",obs_labs[i],"_catchatlgth"), MLE_report)
      dimnames(obs_df) = dimnames(pred_df) = list(c(paste0("M_",length_bins), paste0("F_",length_bins)), regions, years)
      dimnames(obs_indicator) = list(regions, years)
      NA_ndx = which(obs_indicator == 0, arr.ind = T)
      if(nrow(NA_ndx) > 0) {
        for(i in 1:nrow(NA_ndx)) {
          obs_df[,NA_ndx[i,1], NA_ndx[i,2]] = NA
          pred_df[,NA_ndx[i,1], NA_ndx[i,2]] = NA
        }
      }
      molten_obs = reshape2::melt(obs_df)
      molten_pred = reshape2::melt(pred_df)
      colnames(molten_obs) = c("S_Length", "Region", "Year", "Observed")
      molten_obs$Predicted = molten_pred$value
      if(is.null(region_key)) {
        molten_obs$Region = paste0("Region ", molten_obs$Region)
      } else {
        molten_obs$Region = region_key$area[match(molten_obs$Region, (region_key$TMB_ndx + 1))]
      }
      molten_obs$Length = as.numeric(substring(molten_obs$S_Length, first = 3))
      molten_obs$Sex = ifelse(substring(molten_obs$S_Length, first = 0, last = 1) == "M", "Male", "Female")
      molten_obs$label = obs_labs[i]
      full_df = rbind(full_df, molten_obs)
    }
  }

  ## multiple predicted proportions by effective sample size
  full_df= full_df %>% group_by(Year, Region, label) %>% mutate(N_eff = sum(Observed), Observed_mean_prop = Observed / N_eff)

  full_df= full_df %>% group_by(Year, Region, label, Sex) %>% summarise(Ey = sum(Length * Predicted), Oy = sum(Length * Observed_mean_prop), E_squared_y = sum(Length^2 * Predicted), N_eff = mean(N_eff))
  full_df$Ry = full_df$Oy - full_df$Ey
  full_df$SEy = sqrt((full_df$E_squared_y - full_df$Ey^2) / full_df$N_eff)
  full_df$'Std.res' <- (full_df$Oy - full_df$Ey)/full_df$SEy
  ## I think this is the final Francis weighting value TODO: to check
  Nmult <- 1 / var(full_df$'Std.res',na.rm=TRUE)
  # Find the adjusted confidence intervals
  full_df$ObsloAdj <- full_df$Oy - 2 * full_df$SEy / sqrt(Nmult)
  full_df$ObshiAdj <- full_df$Oy + 2 * full_df$SEy / sqrt(Nmult)

  if(!is.null(subset_years)) {
    full_df = full_df %>% filter(Year %in% subset_years)
  }
  if(sex == "male")
    full_df = full_df %>% filter(Sex == "Male")
  if(sex == "female")
    full_df = full_df %>% filter(Sex == "Female")

  ## plot
  gplt = ggplot(full_df, aes(x = Year)) +
    geom_point(aes(y = Oy, col = "Observed", shape = Sex, group = Sex)) +
    geom_line(aes(y = Ey, col = "Predicted", linetype = Sex, group = Sex), linewidth= 1.1) +
    geom_errorbar(aes(ymin=ObsloAdj, ymax=ObshiAdj, col = "Observed"), width=.2, position=position_dodge(.9)) +
    guides( linewidth = "none") +
    labs(y = "Mean length", col = "", linetype = "") +
    facet_grid(label ~ Region) +
    theme_bw()
  return(gplt)
}

#'
#' plot_catch_fit
#' @param MLE_report a list that is output from obj$report() usually once an optimsation routine has been done.
#' @param region_key data.frame with colnames area and TMB_ndx for providing real region names to objects
#' @return ggplot2 object that will plot if an observation occurs in a year and region
#' @export

plot_catch_fit = function(MLE_report, region_key = NULL) {
  full_df = get_catches(MLE_report, region_key)
  gplt = ggplot() +
    geom_point(data = full_df %>% filter(type == "Observed"), aes(x = Year, y = Catch, col = type)) +
    geom_line(data = full_df %>% filter(type == "Predicted"), aes(x = Year, y = Catch, col = type), linewidth= 1.1, linetype = "dashed") +
    guides( linewidth = "none", linetype = "none") +
    labs(y = "Catch", col = "", linetype = "") +
    facet_wrap(label~Region) +
    theme_bw()
  return(gplt)
}
#'
#' plot_index_fit
#' @param MLE_report a list that is output from obj$report() usually once an optimsation routine has been done.
#' @param region_key data.frame with colnames area and TMB_ndx for providing real region names to objects
#' @return ggplot2 object that will plot if an observation occurs in a year and region
#' @export
plot_index_fit = function(MLE_report, region_key = NULL) {

  full_df = get_index(MLE_report, region_key)

  gplt = ggplot(full_df, aes(x = Year)) +
    geom_point(aes(y = Observed, col = "Observed")) +
    geom_line(aes(y = Predicted, col = "Predicted"), linewidth= 1.1, linetype = "dashed") +
    geom_errorbar(aes(ymin=L_CI, ymax=U_CI, col = "Observed"), width=.2, position=position_dodge(.9)) +
    guides( linewidth = "none", linetype = "none") +
    labs(y = "Index", col = "", linetype = "") +
    facet_wrap(~Region, ncol = 2, scales = "free_y") +
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

      ## now link to subsequent recoveries
      for(release_year_ndx in 1:(MLE_report$n_years_to_retain_tagged_cohorts_for)) {
        release_year = recovery_years[y_ndx] - (release_year_ndx - 1)
        for(release_region_ndx in 1:length(regions)) {

          release_event_ndx = get_tag_release_ndx(release_region_ndx, release_year_ndx, MLE_report$n_regions)

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
#' @return data frame in long format
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

#'
#' get_index
#' @param MLE_report a list that is output from obj$report() usually once an optimsation routine has been done.
#' @param region_key data.frame with colnames area and TMB_ndx for providing real region names to objects
#' @return ggplot2 object that will plot if an observation occurs in a year and region
#' @export
get_index = function(MLE_report, region_key = NULL) {
  years = MLE_report$years
  regions = 1:MLE_report$n_regions
  dimnames(MLE_report$obs_srv_dom_ll_bio) = dimnames(MLE_report$pred_srv_dom_ll_bio) =   dimnames(MLE_report$obs_srv_dom_ll_se) = list(regions, years)
  MLE_report$obs_srv_dom_ll_bio[MLE_report$srv_dom_ll_bio_indicator == 0] = NA
  MLE_report$obs_srv_dom_ll_se[MLE_report$srv_dom_ll_bio_indicator == 0] = NA
  MLE_report$pred_srv_dom_ll_bio[MLE_report$srv_dom_ll_bio_indicator == 0] = NA

  index_obs = reshape2::melt(MLE_report$obs_srv_dom_ll_bio)
  index_se = reshape2::melt(MLE_report$obs_srv_dom_ll_se)
  index_fit = reshape2::melt(MLE_report$pred_srv_dom_ll_bio)
  colnames(index_obs) = c("Region", "Year", "Observed")
  colnames(index_fit) = c("Region", "Year", "Predicted")
  colnames(index_se) = c("Region", "Year", "SE")
  CIs = lognormal_CI(index_obs$Observed, sigma = index_se$SE, CI = 0.95)
  index_obs$Predicted = index_fit$Predicted
  index_obs$SE = index_se$SE
  full_df = index_obs
  full_df$U_CI = CIs$upper
  full_df$L_CI = CIs$lower

  if(is.null(region_key)) {
    full_df$Region = paste0("Region ", full_df$Region)
  } else {
    full_df$Region = region_key$area[match(full_df$Region, (region_key$TMB_ndx + 1))]
  }
  full_df = full_df %>% filter(!is.na(Observed))
  return(full_df)
}
