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
    full_df = full_df %>% dplyr::filter(Year %in% subset_years)
  }
  if(sex == "male")
    full_df = full_df %>% dplyr::filter(Sex == "Male")
  if(sex == "female")
    full_df = full_df %>% dplyr::filter(Sex == "Female")

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
    full_df = full_df %>% dplyr::filter(Year %in% subset_years)
  }
  if(sex == "male")
    full_df = full_df %>% dplyr::filter(Sex == "Male")
  if(sex == "female")
    full_df = full_df %>% dplyr::filter(Sex == "Female")

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
    full_df = full_df %>% dplyr::filter(Year %in% subset_years)
  }
  if(sex == "male")
    full_df = full_df %>% dplyr::filter(Sex == "Male")
  if(sex == "female")
    full_df = full_df %>% dplyr::filter(Sex == "Female")

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
    full_df = full_df %>% dplyr::filter(Year %in% subset_years)
  }
  if(sex == "male")
    full_df = full_df %>% dplyr::filter(Sex == "Male")
  if(sex == "female")
    full_df = full_df %>% dplyr::filter(Sex == "Female")

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
    geom_point(data = full_df %>% dplyr::filter(type == "Observed"), aes(x = Year, y = Catch, col = type)) +
    geom_line(data = full_df %>% dplyr::filter(type == "Predicted"), aes(x = Year, y = Catch, col = type), linewidth= 1.1, linetype = "dashed") +
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
  full_df = full_df %>% dplyr::filter(!is.na(Observed))
  return(full_df)
}
