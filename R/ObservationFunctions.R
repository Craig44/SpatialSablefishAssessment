#' get_negloglike get a data frame of negative log likelihoods
#' @param MLE_report a list that is output from obj$report() usually once an optimsation routine has been done.
#' @return data frame with age-frequency info
#' @export
get_negloglike <- function(MLE_report) {
  nll_df = data.frame(negloglike = round(MLE_report$nll,4), observations = c("Fixed AF", "Trawl LF", "Fixed LF","Survey AF","Survey abund","Fixed catch","Trawl catch","Tag recovery", "Recruitment", "Initialisation devs", "posfun penalty"))
  return(nll_df)
}

#' get_AF accessor function to get age-frequency data
#' @param MLE_report a list that is output from obj$report() usually once an optimsation routine has been done.
#' @param label character labeling the observation you want to plot. See below for options
#' \itemize{
#'   \item `fixed`
#'   \item `srv_dom_ll`
#' }
#' @param subset_years vector of years to plot it for
#' @param sex character that allows users to specify if the want sex specific plots
#' \itemize{
#'   \item `both`
#'   \item `male`
#'   \item `female`
#' }
#' @param region_key data.frame with colnames area and TMB_ndx for providing real region names to objects
#' @return data frame with age-frequency info
#' @export

get_AF <- function(MLE_report, label = "fixed", subset_years = NULL, sex = "both", region_key = NULL) {
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
  full_df$label = label

  ## remove rows that have observed NA
  full_df = full_df %>% dplyr::filter(!is.na(Observed))

  return(full_df)
}

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
  full_df = get_AF(MLE_report = MLE_report, label = label, subset_years = subset_years, sex = sex, region_key = region_key)
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
    full_df = get_AF(MLE_report = MLE_report, label = label, subset_years = subset_years, sex = sex, region_key = region_key)
  } else {
    full_df = NULL;
    obs_labs = c("srv_dom_ll","fixed")
    for(i in 1:length(obs_labs)) {
      tmp_df = get_AF(MLE_report = MLE_report, label = obs_labs[i], subset_years = subset_years, sex = sex, region_key = region_key)
      full_df = rbind(full_df, tmp_df)
    }
  }
  ## drop NA's
  full_df = full_df %>% filter(!is.na(Observed))
  ## multiple predicted proportions by effective sample size
  full_df= full_df %>% group_by(Year, Region, label) %>% mutate(N_eff = sum(Observed), Observed_prop = Observed / N_eff, Predicted_prop = Predicted / sum(Predicted))

  full_df= full_df %>% group_by(Year, Region, label, Sex) %>% summarise(Ey = sum(Age * Predicted_prop), Oy = sum(Age * Observed_prop), E_squared_y = sum(Age^2 * Predicted_prop), N_eff = mean(N_eff))
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
    geom_point(aes(y = Oy, col = "Observed", shape = Sex, group = Sex), size = 1.6) +
    geom_line(aes(y = Ey, col = "Predicted", linetype = Sex, group = Sex), linewidth= 1.2) +
    geom_point(aes(y = Ey, col = "Predicted", shape = Sex, group = Sex), size = 1) +
    geom_errorbar(aes(ymin=ObsloAdj, ymax=ObshiAdj, col = "Observed"), width=.2, position=position_dodge(.9)) +
    guides( linewidth = "none") +
    labs(y = "Mean age", col = "", linetype = "") +
    facet_grid(label ~ Region) +
    theme_bw()
  return(gplt)
}
#'
#' get_LF
#' @param MLE_report a list that is output from obj$report() usually once an optimsation routine has been done.
#' @param label character labeling the observation you want to plot. See below for options
#' \itemize{
#'   \item trwl
#'   \item fixed
#' }
#' @param subset_years vector of years to plot it for
#' @param sex character that allows users to specify if the want sex specific plots
#' @param region_key data.frame with colnames area and TMB_ndx for providing real region names to objects
#' @return long data frame with LF infor
#' @export
get_LF = function(MLE_report, label = "fixed", subset_years = NULL, sex = "both", region_key = NULL) {
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
  full_df$label = label
  ## remove rows that have observed NA
  full_df = full_df %>% dplyr::filter(!is.na(Observed))

  return(full_df)
}



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
  full_df = get_LF(MLE_report = MLE_report, label = label, subset_years = subset_years, sex = sex, region_key = region_key)

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
    full_df = get_LF(MLE_report = MLE_report, label = label, subset_years = subset_years, sex = sex, region_key = region_key)
  } else {
    full_df = NULL;
    obs_labs = c("trwl","fixed")
    for(i in 1:length(obs_labs)) {
      tmp_df = get_LF(MLE_report = MLE_report, label = obs_labs[i], subset_years = subset_years, sex = sex, region_key = region_key)
      full_df = rbind(full_df, tmp_df)
    }
  }

  ## multiple predicted proportions by effective sample size
  full_df= full_df %>% group_by(Year, Region, label) %>% mutate(N_eff = sum(Observed), Observed_prop = Observed / N_eff, Predicted_prop = Predicted / sum(Predicted))

  full_df= full_df %>% group_by(Year, Region, label, Sex) %>% summarise(Ey = sum(Length * Predicted_prop), Oy = sum(Length * Observed_prop), E_squared_y = sum(Length^2 * Predicted_prop), N_eff = mean(N_eff))
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
    geom_point(aes(y = Oy, col = "Observed", shape = Sex, group = Sex), size = 1.6) +
    geom_line(aes(y = Ey, col = "Predicted", linetype = Sex, group = Sex), linewidth= 1.1) +
    geom_point(aes(y = Ey, col = "Predicted", shape = Sex, group = Sex), size = 1.2) +
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
  ## convert the SE of an estimator to a standard deviation that is the
  ## right scale for the lognormal distribution
  ## first calculate CV = sigma/mean then pass this to the log_sigma function
  index_se$SE = log_sigma(index_se$SE / index_se$Observed)

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


#'
#'
#' Francis_reweighting
#' @details Calculates stage two weights for each region and observation using method TA1.8 from \insertCite{francis2011data}{SpatialSablefishAssessment}
#' @param MLE_report a list that is output from obj$report() usually once an optimsation routine has been done.
#' @param region_key data.frame with colnames area and TMB_ndx for providing real region names to objects
#' @return named list length_multipliers
#' \itemize{
#'   \item length_multipliers stage 2 weights for length observations
#'   \item age_multipliers  stage 2 weights for age observations
#'   \item mean_age_df mean age info used to calculate weights
#'   \item mean_len_df mean length info used to calculate weights
#' }
#' @export
#' @references
#' \insertAllCited{}

Francis_reweighting <- function(MLE_report, region_key = NULL) {
  years = data$years
  regions = 1:data$n_regions
  age_obs_label = c("fixed","srv_dom_ll")
  len_obs_label = c("fixed","trwl")
  age_comp_df = len_comp_df = NULL
  for(age_ndx in 1:length(age_obs_label))
    age_comp_df = rbind(age_comp_df, get_AF(MLE_report = MLE_report, region_key = region_key, label = age_obs_label[age_ndx], sex = "both"))
  for(len_ndx in 1:length(len_obs_label))
    len_comp_df = rbind(len_comp_df, get_LF(MLE_report = MLE_report, region_key = region_key, label = len_obs_label[len_ndx], sex = "both"))
  ## drop NA observed
  len_comp_df = len_comp_df %>% filter(!is.na(Observed))
  age_comp_df = age_comp_df %>% filter(!is.na(Observed))
  ## get effective sample size
  len_comp_df = len_comp_df %>% group_by(Year, Region, label) %>% mutate(Nassumed = sum(Observed), O_prop = Observed / Nassumed, P_prop = Predicted / Nassumed,
                                                                         O_length = O_prop * Length, P_length = P_prop * Length,  P_length_sq = P_prop * Length^2)
  age_comp_df = age_comp_df %>% group_by(Year, Region, label) %>% mutate(Nassumed = sum(Observed), O_prop = Observed / Nassumed, P_prop = Predicted / Nassumed,
                                                                         O_age = O_prop * Age, P_age = P_prop * Age,  P_age_sq = P_prop * Age^2)
  ## summarise for each obs and year
  mean_len_df = len_comp_df %>% group_by(Year, Region, label) %>% summarise(O_mean_length = sum(O_length),P_mean_length = sum(P_length), P_mean_length_sq = sum(P_length_sq), Nassumed = mean(Nassumed)) %>%
    mutate(stand_mean_length = sqrt(P_mean_length_sq - P_mean_length^2), resid_mean_length = O_mean_length - P_mean_length)
  mean_age_df = age_comp_df %>% group_by(Year, Region, label) %>% summarise(O_mean_age = sum(O_age),P_mean_age = sum(P_age), P_mean_age_sq = sum(P_age_sq), Nassumed = mean(Nassumed)) %>%
    mutate(stand_mean_age = sqrt(P_mean_age_sq - P_mean_age^2),resid_mean_age = O_mean_age - P_mean_age)
  ## get the multiplier over all years for each observation
  length_multipliers = mean_len_df %>% group_by(Region, label) %>% summarise(multiplier = 1/var(resid_mean_length * sqrt(Nassumed)/stand_mean_length, na.rm = T))
  age_multipliers = mean_age_df %>% group_by(Region, label) %>% summarise(multiplier = 1/var(resid_mean_age * sqrt(Nassumed)/stand_mean_age, na.rm = T))
  return(list(length_multipliers = length_multipliers, age_multipliers = age_multipliers, mean_age_df = mean_age_df, mean_len_df = mean_len_df))
}


#'
#'
#' simulate_observations
#' @details Simulate observations conditional on MLE estimates
#' @param obj A TMB object that has been build using `TMB::MakeADFun`
#' @param n_sims an integer specifying how many simulated data sets you want
#' @param region_key data.frame with colnames area and TMB_ndx for providing real region names to objects
#' @return named list containing simualted observations for all key obsevations
#' @export

simulate_observations <- function(obj, n_sims = 200, region_key = NULL) {
  fixed_effect_pars = get_tmb_fixed_effects(obj)
  all_pars = obj$env$last.par.best
  sd_rep = sdreport(obj, getJointPrecision = T)
  sim_pars =  MASS::mvrnorm(n = n_sims, mu = fixed_effect_pars, Sigma = sd_rep$cov.fixed)
  sim_srv_bio = sim_srv_AF = sim_fixed_AF = sim_fixed_LF = sim_trwl_LF = sim_tag_recovery = NULL
  for(sim_iter in 1:n_sims) {
    if(sim_iter %% 50 == 0)
      cat("simulation iteration: ", sim_iter, "\n")
    ## simualte
    this_sim = obj$simulate(par = sim_pars[sim_iter,], complete = T)
    ## store sim obs
    # survey biomass
    index_df = get_index(this_sim, region_key = region_key)
    index_df$sim = sim_iter
    sim_srv_bio = rbind(sim_srv_bio, index_df)
    # survey AF
    srv_AF = get_AF(MLE_report = this_sim, label = "srv_dom_ll", region_key = region_key)
    srv_AF$sim = sim_iter
    sim_srv_AF = rbind(sim_srv_AF, srv_AF)
    # Fixed AF
    fixed_AF = get_AF(MLE_report = this_sim, label = "fixed", region_key = region_key)
    fixed_AF$sim = sim_iter
    sim_fixed_AF = rbind(sim_fixed_AF, fixed_AF)
    # Fixed LF
    fixed_LF = get_LF(MLE_report = this_sim, label = "fixed", region_key = region_key)
    fixed_LF$sim = sim_iter
    sim_fixed_LF = rbind(sim_fixed_LF, fixed_LF)
    # Trawl LF
    trwl_LF = get_LF(MLE_report = this_sim, label = "trwl", region_key = region_key)
    trwl_LF$sim = sim_iter
    sim_trwl_LF = rbind(sim_trwl_LF, trwl_LF)
    # Tag data
    tag_data = get_tag_recovery_obs_fitted_values(MLE_report = this_sim, region_key = region_key)
    tag_data$sim = sim_iter
    sim_tag_recovery = rbind(sim_tag_recovery, tag_data)
  }
  return(list(sim_tag_recovery = sim_tag_recovery, sim_trwl_LF = sim_trwl_LF, sim_fixed_LF = sim_fixed_LF, sim_fixed_AF = sim_fixed_AF, sim_srv_AF = sim_srv_AF, sim_srv_bio = sim_srv_bio))
}
