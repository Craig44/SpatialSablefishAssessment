#' get_negloglike get a data frame of negative log likelihoods
#' @param MLE_report a list that is output from obj$report() usually once an optimsation routine has been done.
#' @return data frame with age-frequency info
#' @export
get_negloglike <- function(MLE_report) {
  nll_df= NULL
  if(MLE_report$model_type == 0) {
    nll_label = c("fishery-ll-age_comp",
    "fishery-trwl-male_length_comp",
    "fishery-trwl-female_length_comp",
    "srv-Domestic-ll_index",
    "srv-Japanese-ll_index",
    "fishery-ll_index",
    "srv-Domestic-ll-age_comp",
    "srv-Domestic-ll-male_length_comp",
    "srv-Domestic-ll-female_length_comp",
    "srv-Japanese-ll-age_comp",
    "srv-Japanese-ll-male_length_comp",
    "srv-Japanese-ll-female_length_comp",
    "srv-GOA-trwl-age_comp",
    "srv-GOA-trwl-male_length_comp",
    "srv-GOA-trwl-female_length_comp",
    "fishery-ll-male_length_comp",
    "fishery-ll-female_length_comp",
    "srv-GOA-trwl-index",
    "Historic-Japanese-Fishery-index",
    "Historic-Japanese-Fishery-length_comp",
    "fishery-ll-catch_Sum_of_squares",
    "fishery-trwlcatch_Sum_of_squares",
    "Recruitment-penalty",
    "Longline_F_penalty",
    "Trawl_F_penalty",
    "Q_priors",
    "M_priors",
    "M_penalty"
    )

    nll_df = data.frame(negloglike = round(MLE_report$nll,4), weighted_negloglike = round(MLE_report$nll_weighted, 4), weights = round(MLE_report$nll_weighted / MLE_report$nll, 4),
                        observations = nll_label, distribution = NA)

  } else {
    AF_fixed_like = ifelse(MLE_report$fixed_catchatage_comp_likelihood == 0, "Multinomial", "Dirichlet-Multinomial")
    LF_fixed_like = ifelse(MLE_report$fixed_catchatlgth_comp_likelihood == 0, "Multinomial", "Dirichlet-Multinomial")
    AF_srv_like = paste(ifelse(MLE_report$srv_catchatage_comp_likelihood == 0, "Multinomial", "Dirichlet-Multinomial"), collapse = ", ")
    LF_trwl_like = ifelse(MLE_report$trwl_catchatlgth_comp_likelihood == 0, "Multinomial", "Dirichlet-Multinomial")
    tag_likelihood = switch(MLE_report$tag_likelihood + 1,
                            "Poissson",
                            "Negative Binomial",
                            "Multinomial")

    nll_df = data.frame(negloglike = round(MLE_report$nll,4),
                        observations = c("Fixed AF", "Trawl LF", "Fixed LF","Survey AF","Survey abund","Fixed catch","Trawl catch","Tag recovery", "Recruitment", "Initialisation devs", "posfun penalty", "F-penalty"),
                        distribution = c(AF_fixed_like, LF_trwl_like, LF_fixed_like, AF_srv_like, "lognormal", "lognormal", "lognormal", tag_likelihood, "lognormal", "lognormal", "", "SSE"))
  }
  return(nll_df)
}
#' get_n_datasets get a data frame of the number of data sets
#' @param MLE_report a list that is output from obj$report() usually once an optimsation routine has been done.
#' @return data frame with age-frequency info
#' @export
get_n_datasets <- function(MLE_report) {
  if(MLE_report$model_type == 0) {
    cat("Skipping this function has not been modified for the model type: Assessment")
    return(NULL)
  }
  years = data$years
  regions = 1:data$n_regions
  surveys = 1:data$n_surveys

  dimnames(data$fixed_catchatage_indicator) = dimnames(data$fixed_catchatlgth_indicator) = dimnames(data$trwl_catchatlgth_indicator) = dimnames(data$srv_bio_indicator) = list(regions, years)
  dimnames(data$srv_catchatage_indicator) = list(regions, years, surveys)
  dimnames(data$tag_recovery_indicator) = list(1:dim(data$tag_recovery_indicator)[1], regions, years[which(data$tag_recovery_indicator_by_year == 1)])
  fixed_catchatage = reshape2::melt(data$fixed_catchatage_indicator)
  fixed_catchatlgth  = reshape2::melt(data$fixed_catchatlgth_indicator)
  trwl_catchatlgth = reshape2::melt(data$trwl_catchatlgth_indicator)
  srv_catchatage = reshape2::melt(data$srv_catchatage_indicator)
  srv_bio = reshape2::melt(data$srv_bio_indicator)
  tag_recovery_detailed = NULL
  if(sum(data$tag_recovery_indicator) != 0) {
    tag_recovery_detailed = reshape2::melt(data$tag_recovery_indicator)
    colnames(tag_recovery_detailed) = c("Tag release", "Region", "Year", "indicator")
    tag_recovery_detailed$label = "Tag recovery"
    ## collapse tag recoveries across release events
    tag_recovery_detailed = tag_recovery_detailed %>% group_by(Region, Year, label) %>% summarise(indicator = ifelse(sum(indicator)>0, 1, 0))

  }
  colnames(fixed_catchatage) = colnames(fixed_catchatlgth) = colnames(trwl_catchatlgth) =  colnames(srv_bio) = c("Region", "Year", "indicator")
  colnames(srv_catchatage) = c("Region", "Year", "Survey","indicator")
  ## tag releases
  tag_release_df = NULL
  if((sum(data$male_tagged_cohorts_by_age) + sum(data$female_tagged_cohorts_by_age)) > 0) {
    dimnames(data$male_tagged_cohorts_by_age) = dimnames(data$female_tagged_cohorts_by_age) = list(data$ages, regions,  data$years[which(data$tag_release_event_this_year == 1)])
    tag_releases_m = reshape2::melt(data$male_tagged_cohorts_by_age)
    tag_releases_f = reshape2::melt(data$female_tagged_cohorts_by_age)
    tag_releases = rbind(tag_releases_m, tag_releases_f)
    colnames(tag_releases) = c("Age", "Region", "Year", "releases")
    tag_release_df = tag_releases %>% group_by(Region, Year) %>% summarise(indicator = ifelse(sum(releases) > 0, 1, 0))
    tag_release_df$label = "Tag Releases"
  }
  fixed_catchatlgth$label = "Fishery Fixed LF"
  fixed_catchatage$label = "Fishery Fixed AF"
  trwl_catchatlgth$label = "Fishery Trawl LF"
  srv_catchatage$label = "Survey LL AF"
  srv_bio$label = "Survey LL Biomass"
  ## combine
  full_df = rbind(fixed_catchatage, trwl_catchatlgth, srv_catchatage, srv_bio, fixed_catchatlgth, tag_recovery_detailed, tag_release_df)

  if(is.null(region_key)) {
    full_df$Region = paste0("Region ", full_df$Region)
  } else {
    full_df$Region = region_key$area[match(full_df$Region, (region_key$TMB_ndx + 1))]
  }

  full_df$indicator = ifelse(full_df$indicator == 0, NA, 1)

  AF_fixed_like = ifelse(MLE_report$fixed_catchatage_comp_likelihood == 0, "Multinomial", "Dirichlet-Multinomial")
  LF_fixed_like = ifelse(MLE_report$fixed_catchatlgth_comp_likelihood == 0, "Multinomial", "Dirichlet-Multinomial")
  AF_srv_like = ifelse(MLE_report$srv_catchatage_comp_likelihood == 0, "Multinomial", "Dirichlet-Multinomial")
  LF_trwl_like = ifelse(MLE_report$trwl_catchatlgth_comp_likelihood == 0, "Multinomial", "Dirichlet-Multinomial")
  tag_likelihood = switch(MLE_report$tag_likelihood + 1,
                          "Poissson",
                          "Negative Binomial",
                          "Multinomial")

  nll_df = data.frame(negloglike = round(MLE_report$nll,4),
                      observations = c("Fixed AF", "Trawl LF", "Fixed LF","Survey AF","Survey abund","Fixed catch","Trawl catch","Tag recovery", "Recruitment", "Initialisation devs", "posfun penalty"),
                      distribution = c(AF_fixed_like, LF_trwl_like, LF_fixed_like, AF_srv_like, "lognormal", "lognormal", "lognormal", tag_likelihood, "lognormal", "lognormal", ""))
  return(nll_df)
}
#' get_AF accessor function to get age-frequency data
#' @param MLE_report a list that is output from obj$report() usually once an optimsation routine has been done.
#' @param observation character labeling the observation you want to plot. See below for options. Ignored if `MLE_report$model_type == 0` and will return all age comps.
#' \itemize{
#'   \item `fixed`
#'   \item `srv`
#' }
#' @param subset_years vector of years to get observation for it. Ignored if `MLE_report$model_type == 0`.
#' @param sex character that allows users to specify if the want sex specific plots. Ignored if `MLE_report$model_type == 0`.
#' \itemize{
#'   \item `both`
#'   \item `male`
#'   \item `female`
#' }
#' @param region_key data.frame with colnames area and TMB_ndx for providing real region names to objects
#' @param survey_labels character vector for each of the n_surveys
#' @return data frame with age-frequency info
#' @export

get_AF <- function(MLE_report, observation = "fixed", subset_years = NULL, sex = "both", region_key = NULL, survey_labels = NULL) {
  full_df = NULL
  if(MLE_report$model_type == 0) {
    #
    age_obs = c("obs_srv_nmfs_trwl_age", "obs_srv_jap_ll_age", "obs_srv_dom_ll_age", "obs_ll_catchatage")
    age_indicators = c("srv_nmfs_trwl_age_indicator", "srv_jap_ll_age_indicator", "srv_dom_ll_age_indicator", "ll_catchatage_indicator")
    age_pred = c("pred_srv_nmfs_trwl_age", "pred_srv_jap_ll_age", "pred_srv_dom_ll_age", "pred_ll_catchatage")
    age_obs_label = c("NMFS survey", "Japanese LL survey", "Domestic LL survey", "Fixed gear fishery")
    for(i in 1:length(age_indicators)) {
      obs_indicator = get(age_indicators[i], MLE_report)
      if(sum(obs_indicator) == 0)
        next;
      obs_df = get(age_obs[i], MLE_report)
      pred_df = get(age_pred[i], MLE_report)
      label = age_obs_label[i]
      obs_years = MLE_report$years[obs_indicator == 1]
      dimnames(obs_df) = dimnames(pred_df) = list(MLE_report$ages,obs_years)
      molten_obs = reshape2::melt(obs_df)
      molten_pred = reshape2::melt(pred_df)
      colnames(molten_obs) = c("Age", "Year", "Observed")
      molten_obs$Predicted = molten_pred$value
      if(is.null(region_key)) {
        molten_obs$Region = paste0("Region ", 1)
      } else {
        molten_obs$Region = region_key$area[1]
      }
      molten_obs$observation = label
      full_df = rbind(full_df, molten_obs)
    }
    ## multiple predicted proportions by effective sample size
    full_df= full_df %>% group_by(Year, Region, observation) %>% mutate(Predicted = Predicted * sum(Observed))
    full_df$Sex = "Combined"
  } else {

    if(!observation %in% c("fixed","srv"))
      stop("observation not one of the expected values.")
    if(!sex %in% c("both", "male", "female"))
      stop('sex not one of the expected values. Expected one of the following "both", "male", "female"')
    surveys = paste0("Survey ", 1:MLE_report$n_surveys)
    if(!is.null(survey_labels))
      surveys = survey_labels
    years = MLE_report$years
    regions = 1:MLE_report$n_regions
    ages = MLE_report$ages
    ## get objects
    obs_indicator = get(paste0(observation,"_catchatage_indicator"), MLE_report)
    obs_df = get(paste0("obs_",observation,"_catchatage"), MLE_report)
    pred_df = get(paste0("pred_",observation,"_catchatage"), MLE_report)
    if(observation == "srv") {
      dimnames(obs_df) = dimnames(pred_df) = list(c(paste0("M_",ages), paste0("F_",ages)), regions, years, surveys)
      dimnames(obs_indicator) = list(regions, years, surveys)
      NA_ndx = which(obs_indicator == 0, arr.ind = T)
      ## TODO: this will need to be tightened up to deal with multiple surveys
      if(nrow(NA_ndx) > 0) {
        for(i in 1:nrow(NA_ndx)) {
          ##cat("i = ", i, " ", sum(obs_df[,NA_ndx[i,1], NA_ndx[i,2],NA_ndx[i,3]]), "\n");
          obs_df[,NA_ndx[i,1], NA_ndx[i,2],NA_ndx[i,3]] = NA
          pred_df[,NA_ndx[i,1], NA_ndx[i,2],NA_ndx[i,3]] = NA
        }
      }
      molten_obs = reshape2::melt(obs_df)
      molten_pred = reshape2::melt(pred_df)
      colnames(molten_obs) = c("S_Age", "Region", "Year", "Survey", "Observed")
    } else {
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
      molten_obs$Survey = factor(1) ## add a dummy variable so the code can be streamlined with suryve AF
    }
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
    full_df= full_df %>% group_by(Year, Region, Survey) %>% mutate(Predicted = Predicted * sum(Observed))

    if(sex == "male")
      full_df = full_df %>% dplyr::filter(Sex == "Male")
    if(sex == "female")
      full_df = full_df %>% dplyr::filter(Sex == "Female")
    full_df$observation = observation

    ## remove rows that have observed NA
    full_df = full_df %>% dplyr::filter(!is.na(Observed))
  }

  if(!is.null(subset_years))
    full_df = full_df %>% dplyr::filter(Year %in% subset_years)

  return(full_df)
}

#
#
#'
#' plot_AF
#' @param MLE_report a list that is output from obj$report() usually once an optimsation routine has been done.
#' @param observation character labeling the observation you want to plot. See below for options
#' \itemize{
#'   \item `trwl`
#'   \item `fixed`
#'   \item `srv`
#' }
#' If `MLE_report$model_type == 0` i.a., assessment model
#' \itemize{
#'   \item `nmfs`
#'   \item `srv_dom_ll`
#'   \item `srv_jap_ll`
#'   \item `fixed_gear`
#' }
#' @param subset_years vector of years to plot it for
#' @param sex character that allows users to specify if the want sex specific plots
#' \itemize{
#'   \item `male`
#'   \item `female`
#'   \item `both`
#' }
#' @param region_key data.frame with colnames area and TMB_ndx for providing real region names to objects
#' @param survey_labels character vector for each of the n_surveys
#' @return ggplot2 object that will plot if an observation occurs in a year and region
#' @export
plot_AF = function(MLE_report, observation = "srv", subset_years = NULL, sex = "both", region_key = NULL, survey_labels = NULL) {
  gplt = NULL
  if(MLE_report$model_type == 0) {
    full_df = get_AF(MLE_report = MLE_report, subset_years = subset_years, region_key = region_key)
    if(!observation %in% c("nmfs", "srv_dom_ll", "srv_jap_ll", "fixed_gear"))
      stop(paste0('observation, needs to be one of the following "nmfs", "srv_dom_ll", "srv_jap_ll", "fixed_gear"'))
    this_obs = NULL
    if(observation == "nmfs") {
      this_obs = full_df %>% filter(observation == "NMFS survey")
    } else if (observation == "srv_dom_ll") {
      this_obs = full_df %>% filter(observation == "Domestic LL survey")
    } else if (observation == "srv_jap_ll") {
      this_obs = full_df %>% filter(observation == "Japanese LL survey")
    } else if (observation == "fixed_gear") {
      this_obs = full_df %>% filter(observation == "Fixed gear fishery")
    }
    gplt = ggplot(this_obs, aes(x = Age)) +
      geom_point(aes(y = Observed, col = "Observed")) +
      geom_line(aes(y = Predicted, col = "Predicted"), linewidth= 1.1) +
      guides( linewidth = "none") +
      labs(y = "AF", col = "", linetype = "") +
      facet_wrap(~Year) +
      theme_bw()
  } else {
    full_df = get_AF(MLE_report = MLE_report, observation = observation, subset_years = subset_years, sex = sex, region_key = region_key, survey_labels = survey_labels)

    ## plot
    gplt = ggplot(full_df, aes(x = Age)) +
      geom_point(aes(y = Observed, col = "Observed", shape = Sex, group = Sex)) +
      geom_line(aes(y = Predicted, col = "Predicted", linetype = Sex, group = Sex), linewidth= 1.1) +
      guides( linewidth = "none") +
      labs(y = "AF", col = "", linetype = "") +
      facet_grid(Year ~ Region) +
      theme_bw()
  }
  return(gplt)
}
#'
#' get_mean_age
#' @param MLE_report a list that is output from obj$report() usually once an optimsation routine has been done.
#' @param observation character labeling the observation you want to plot. See below for options
#' \itemize{
#'   \item all
#'   \item trwl
#'   \item fixed
#'   \item srv
#' }
#' @param subset_years vector of years to plot it for
#' @param sex character that allows users to specify if the want sex specific plots
#' \itemize{
#'   \item `male`
#'   \item `female`
#'   \item `both`
#' }
#' @param region_key data.frame with colnames area and TMB_ndx for providing real region names to objects
#' @param survey_labels character vector for each of the n_surveys
#' @return a data frame of mean ages
#' @export
get_mean_age = function(MLE_report, observation = "fixed", subset_years = NULL, sex = "both", region_key = NULL, survey_labels = NULL) {
  full_df = NULL
  years = MLE_report$years
  regions = 1:MLE_report$n_regions
  ages = MLE_report$ages
  if(MLE_report$model_type == 0) {
    if(!observation %in% c("all"))
      stop("observation needs to be 'all'")
    full_df = get_AF(MLE_report = MLE_report, subset_years = subset_years, region_key = region_key)
    full_df$Survey = 1;

  } else {
    if(!observation %in% c("all","fixed","srv"))
      stop("observation not one of the expected values.")
    if(!sex %in% c("both", "male", "female"))
      stop('sex not one of the expected values. Expected one of the following "both", "male", "female"')

    ## get objects
    if(observation != "all") {
      full_df = get_AF(MLE_report = MLE_report, observation = observation, subset_years = subset_years, sex = sex, region_key = region_key, survey_labels = survey_labels)
    } else {
      full_df = NULL;
      obs_labs = c("srv","fixed")
      for(i in 1:length(obs_labs)) {
        tmp_df = get_AF(MLE_report = MLE_report, observation = obs_labs[i], subset_years = subset_years, sex = sex, region_key = region_key)
        full_df = rbind(full_df, tmp_df)
      }
    }
  }

  ## drop NA's
  full_df = full_df %>% filter(!is.na(Observed))
  ## multiple predicted proportions by effective sample size
  full_df= full_df %>% group_by(Year, Region, observation, Survey) %>% mutate(N_eff = sum(Observed), Observed_prop = Observed / N_eff, Predicted_prop = Predicted / sum(Predicted))

  full_df= full_df %>% group_by(Year, Region, observation, Sex, Survey) %>% summarise(Ey = sum(Age * Predicted_prop), Oy = sum(Age * Observed_prop), E_squared_y = sum(Age^2 * Predicted_prop), N_eff = mean(N_eff))
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
  return(full_df)
}

#'
#' plot_mean_age
#' @param MLE_report a list that is output from obj$report() usually once an optimsation routine has been done.
#' @param observation character labeling the observation you want to plot. See below for options
#' \itemize{
#'   \item all
#'   \item trwl
#'   \item fixed
#'   \item srv
#' }
#' @param subset_years vector of years to plot it for
#' @param sex character that allows users to specify if the want sex specific plots
#' \itemize{
#'   \item `male`
#'   \item `female`
#'   \item `both`
#' }
#' @param region_key data.frame with colnames area and TMB_ndx for providing real region names to objects
#' @param survey_labels character vector for each of the n_surveys
#' @return ggplot2 object that will plot if an observation occurs in a year and region
#' @export
plot_mean_age = function(MLE_report, observation = "fixed", subset_years = NULL, sex = "both", region_key = NULL, survey_labels = NULL) {
  if(!observation %in% c("all","fixed","srv"))
    stop("observation not one of the expected values.")
  if(!sex %in% c("both", "male", "female"))
    stop('sex not one of the expected values. Expected one of the following "both", "male", "female"')

  full_df = get_mean_age(MLE_report = MLE_report, observation = observation, subset_years = subset_years, sex = sex, region_key = region_key, survey_labels = survey_labels)
  gplt = NULL
  if(observation == "srv") {
    ## plot
    gplt = ggplot(full_df, aes(x = Year)) +
      geom_point(aes(y = Oy, col = "Observed", shape = Sex, group = Sex), size = 1.6) +
      geom_line(aes(y = Ey, col = "Predicted", linetype = Sex, group = Sex), linewidth= 1.2) +
      geom_point(aes(y = Ey, col = "Predicted", shape = Sex, group = Sex), size = 1) +
      geom_errorbar(aes(ymin=ObsloAdj, ymax=ObshiAdj, col = "Observed"), width=.2, position=position_dodge(.9)) +
      guides( linewidth = "none") +
      labs(y = "Mean age", col = "", linetype = "") +
      facet_grid(Survey ~ Region) +
      theme_bw()
  } else {
    ## plot
    gplt = ggplot(full_df, aes(x = Year)) +
      geom_point(aes(y = Oy, col = "Observed", shape = Sex, group = Sex), size = 1.6) +
      geom_line(aes(y = Ey, col = "Predicted", linetype = Sex, group = Sex), linewidth= 1.2) +
      geom_point(aes(y = Ey, col = "Predicted", shape = Sex, group = Sex), size = 1) +
      geom_errorbar(aes(ymin=ObsloAdj, ymax=ObshiAdj, col = "Observed"), width=.2, position=position_dodge(.9)) +
      guides( linewidth = "none") +
      labs(y = "Mean age", col = "", linetype = "") +
      facet_grid(observation ~ Region) +
      theme_bw()
  }

  return(gplt)
}
#'
#' get_LF
#' @param MLE_report a list that is output from obj$report() usually once an optimsation routine has been done.
#' @param observation character labeling the observation you want to plot. See below for options
#' \itemize{
#'   \item trwl
#'   \item fixed
#' }
#' @param subset_years vector of years to plot it for
#' @param sex character that allows users to specify if the want sex specific plots
#' \itemize{
#'   \item `male`
#'   \item `female`
#'   \item `both`
#' }
#' @param region_key data.frame with colnames area and TMB_ndx for providing real region names to objects
#' @return long data frame with LF infor
#' @export
get_LF = function(MLE_report, observation = "fixed", subset_years = NULL, sex = "both", region_key = NULL) {
  if(!sex %in% c("both", "male", "female"))
    stop('sex not one of the expected values. Expected one of the following "both", "male", "female"')
  length_bins = MLE_report$length_bins

  full_df = NULL
  if(MLE_report$model_type == 0) {
    #
    lgth_obs = c("obs_srv_jap_fishery_ll_lgth", "obs_ll_catchatlgth", "obs_trwl_catchatlgth", "obs_srv_nmfs_trwl_lgth", "obs_srv_dom_ll_lgth", "obs_srv_jap_ll_lgth")
    lgth_indicators = c("srv_jap_fishery_ll_lgth_indicator", "ll_catchatlgth_indicator", "trwl_catchatlgth_indicator", "srv_nmfs_trwl_lgth_indicator" , "srv_dom_ll_lgth_indicator" ,  "srv_jap_ll_lgth_indicator")
    lgth_pred = c("pred_srv_jap_fishery_ll_lgth", "pred_ll_catchatlgth", "pred_trwl_catchatlgth", "pred_srv_nmfs_trwl_lgth", "pred_srv_dom_ll_lgth", "pred_srv_jap_ll_lgth")
    lgth_obs_label = c("Japanese LL fishery", "Fixed gear fishery", "Trawl gear fishery", "NMFS survey","Domestic LL survey", "Japanese LL survey")
    for(i in 1:length(lgth_indicators)) {
      if(sex != "both") {
        sex_label = substring(sex, first =0, last = 1)
        obs_indicator = get(lgth_indicators[i], MLE_report)
        if(sum(obs_indicator) > 0) { ## no data in for this observation
          if(lgth_obs_label[i] == "Japanese LL fishery") {
            obs_df = get(lgth_obs[i], MLE_report)
            pred_df = get(lgth_pred[i], MLE_report)
          } else {
            obs_df = get(paste0(lgth_obs[i], "_", sex_label), MLE_report)
            pred_df = get(paste0(lgth_pred[i], "_", sex_label), MLE_report)
          }
          label = lgth_obs_label[i]
          obs_years = MLE_report$years[obs_indicator == 1]
          dimnames(obs_df) = dimnames(pred_df) = list(MLE_report$length_bins,obs_years)
          molten_obs = reshape2::melt(obs_df)
          molten_pred = reshape2::melt(pred_df)
          colnames(molten_obs) = c("Length", "Year", "Observed")
          molten_obs$Predicted = molten_pred$value
          if(is.null(region_key)) {
            molten_obs$Region = paste0("Region ", 1)
          } else {
            molten_obs$Region = region_key$area[1]
          }
          molten_obs$observation = label
          molten_obs$Sex = sex
          if(lgth_obs_label[i] == "Japanese LL fishery")
            molten_obs$Sex = "Combined"

          full_df = rbind(full_df, molten_obs)
        }
      } else {
        for(sex_ndx in c("male", "female")) {
          sex_label = substring(sex_ndx, first =0, last = 1)
          obs_indicator = get(lgth_indicators[i], MLE_report)
          if(sum(obs_indicator) > 0) { ## no data in for this observation
            if(lgth_obs_label[i] == "Japanese LL fishery") {
              obs_df = get(lgth_obs[i], MLE_report)
              pred_df = get(lgth_pred[i], MLE_report)
            } else {
              obs_df = get(paste0(lgth_obs[i], "_", sex_label), MLE_report)
              pred_df = get(paste0(lgth_pred[i], "_", sex_label), MLE_report)
            }
            label = lgth_obs_label[i]
            obs_years = MLE_report$years[obs_indicator == 1]
            dimnames(obs_df) = dimnames(pred_df) = list(MLE_report$length_bins, obs_years)
            molten_obs = reshape2::melt(obs_df)
            molten_pred = reshape2::melt(pred_df)
            colnames(molten_obs) = c("Length", "Year", "Observed")
            molten_obs$Predicted = molten_pred$value
            if(is.null(region_key)) {
              molten_obs$Region = paste0("Region ", 1)
            } else {
              molten_obs$Region = region_key$area[1]
            }
            molten_obs$observation = label
            molten_obs$Sex = sex_ndx
            if(lgth_obs_label[i] == "Japanese LL fishery")
              molten_obs$Sex = "Combined"
            full_df = rbind(full_df, molten_obs)
          }
        }
      }
    }
    ## multiple predicted proportions by effective sample size
    full_df= full_df %>% group_by(Year, Region, observation, Sex) %>% mutate(Predicted = Predicted * sum(Observed))
    if(!is.null(subset_years))
      full_df = full_df %>% dplyr::filter(Year %in% subset_years)
  } else {
    if(!observation %in% c("fixed","trwl"))
      stop("observation not one of the expected values.")

    years = MLE_report$years
    regions = 1:MLE_report$n_regions
    ## get objects
    obs_indicator = get(paste0(observation,"_catchatlgth_indicator"), MLE_report)
    obs_df = get(paste0("obs_",observation,"_catchatlgth"), MLE_report)
    pred_df = get(paste0("pred_",observation,"_catchatlgth"), MLE_report)
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
    full_df$observation = observation
    ## remove rows that have observed NA
    full_df = full_df %>% dplyr::filter(!is.na(Observed))
  }
  return(full_df)
}



#' plot_LF
#' @param MLE_report a list that is output from obj$report() usually once an optimsation routine has been done.
#' @param observation character labeling the observation you want to plot. See below for options
#' \itemize{
#'   \item trwl
#'   \item fixed
#' }
#' @param subset_years vector of years to plot it for
#' @param sex character that allows users to specify if the want sex specific plots
#' @param region_key data.frame with colnames area and TMB_ndx for providing real region names to objects
#' @return ggplot2 object that will plot if an observation occurs in a year and region
#' @export
plot_LF = function(MLE_report, observation = "fixed", subset_years = NULL, sex = "both", region_key = NULL) {
  full_df = get_LF(MLE_report = MLE_report, observation = observation, subset_years = subset_years, sex = sex, region_key = region_key)

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
#' get_mean_length
#' @param MLE_report a list that is output from obj$report() usually once an optimsation routine has been done.
#' @param observation character labeling the observation you want to plot. See below for options
#' \itemize{
#'   \item all
#'   \item trwl
#'   \item fixed
#' }
#' @param subset_years vector of years to plot it for
#' @param sex character that allows users to specify if the want sex specific plots
#' @param region_key data.frame with colnames area and TMB_ndx for providing real region names to objects
#' @return data frame of mean lengths
#' @export
get_mean_length = function(MLE_report, observation = "fixed", subset_years = NULL, sex = "both", region_key = NULL) {
  if(!observation %in% c("all","fixed","trwl"))
    stop("observation not one of the expected values.")
  if(!sex %in% c("both", "male", "female"))
    stop('sex not one of the expected values. Expected one of the following "both", "male", "female"')

  years = MLE_report$years
  regions = 1:MLE_report$n_regions
  length_bins = MLE_report$length_bins
  ## get objects
  if(observation != "all") {
    full_df = get_LF(MLE_report = MLE_report, observation = observation, subset_years = subset_years, sex = sex, region_key = region_key)
  } else {
    full_df = NULL;
    obs_labs = c("trwl","fixed")
    for(i in 1:length(obs_labs)) {
      tmp_df = get_LF(MLE_report = MLE_report, observation = obs_labs[i], subset_years = subset_years, sex = sex, region_key = region_key)
      full_df = rbind(full_df, tmp_df)
    }
  }

  ## multiple predicted proportions by effective sample size
  if(MLE_report$model_type == 0) {
    full_df= full_df %>% group_by(Year, Region, observation, Sex) %>% mutate(N_eff = sum(Observed), Observed_prop = Observed / N_eff, Predicted_prop = Predicted / sum(Predicted))
    full_df= full_df %>% group_by(Year, Region, observation, Sex) %>% summarise(Ey = sum(Length * Predicted_prop), Oy = sum(Length * Observed_prop), E_squared_y = sum(Length^2 * Predicted_prop), N_eff = mean(N_eff))
  } else  {
    full_df= full_df %>% group_by(Year, Region, observation) %>% mutate(N_eff = sum(Observed), Observed_prop = Observed / N_eff, Predicted_prop = Predicted / sum(Predicted))
    full_df= full_df %>% group_by(Year, Region, observation, Sex) %>% summarise(Ey = sum(Length * Predicted_prop), Oy = sum(Length * Observed_prop), E_squared_y = sum(Length^2 * Predicted_prop), N_eff = mean(N_eff))
  }
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
  return(full_df)
}

#'
#' plot_mean_length
#' @param MLE_report a list that is output from obj$report() usually once an optimsation routine has been done.
#' @param observation character labeling the observation you want to plot. See below for options
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
plot_mean_length = function(MLE_report, observation = "fixed", subset_years = NULL, sex = "both", region_key = NULL) {
  if(!observation %in% c("all","fixed","trwl"))
    stop("observation not one of the expected values.")
  if(!sex %in% c("both", "male", "female"))
    stop('sex not one of the expected values. Expected one of the following "both", "male", "female"')
  full_df = get_mean_length(MLE_report, observation, subset_years, sex, region_key)

  ## plot
  gplt = ggplot(full_df, aes(x = Year)) +
    geom_point(aes(y = Oy, col = "Observed", shape = Sex, group = Sex), size = 1.6) +
    geom_line(aes(y = Ey, col = "Predicted", linetype = Sex, group = Sex), linewidth= 1.1) +
    geom_point(aes(y = Ey, col = "Predicted", shape = Sex, group = Sex), size = 1.2) +
    geom_errorbar(aes(ymin=ObsloAdj, ymax=ObshiAdj, col = "Observed"), width=.2, position=position_dodge(.9)) +
    guides( linewidth = "none") +
    labs(y = "Mean length", col = "", linetype = "") +
    facet_grid(observation ~ Region) +
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
    geom_point(data = full_df %>% dplyr::filter(type == "Observed"), aes(x = Year, y = Catch, col = "Observed"), size = 2.1) +
    geom_line(data = full_df %>% dplyr::filter(type == "Predicted"), aes(x = Year, y = Catch, col = "Predicted"), linewidth= 1.1, linetype = "dashed") +
    guides( linewidth = "none", linetype = "none") +
    labs(y = "Catch", col = "", linetype = "") +
    facet_grid(Region~Fishery) +
    theme_bw() +
    theme(legend.position = "bottom",
          axis.text = element_text(size = 14),
          axis.title = element_text(size = 14),
          strip.text = element_text(size=14),
          plot.title = element_text(size = 20, face = "bold"),
          legend.text = element_text(size=14))
  return(gplt)
}
#'
#' plot_index_fit
#' @param MLE_report a list that is output from obj$report() usually once an optimsation routine has been done.
#' @param region_key data.frame with colnames area and TMB_ndx for providing real region names to objects
#' @param resid
#' \itemize{
#'   \item F - conventional observed with predicted with standard errors
#'   \item T - Pearson residuals with smoother
#' }
#' @return ggplot2 object that will plot if an observation occurs in a year and region
#' @export
plot_index_fit = function(MLE_report, region_key = NULL, resid = F) {
  full_df = get_index(MLE_report, region_key)
  gplt = NULL
  if(MLE_report$model_type == 0) {
    if(resid) {
      gplt = ggplot(full_df, aes(x = Year)) +
        geom_point(aes(y = Pearsons_residuals, col = "Pearson residuals")) +
        geom_smooth(aes(y = Pearsons_residuals, col = "Smoother", fill = "Smoother"), alpha = 0.6) +
        geom_hline(yintercept = c(-2,2), col = "gray60", linetype = "dashed") +
        guides( linewidth = "none", linetype = "none") +
        labs(y = "Index", col = "", linetype = "", fill = "") +
        facet_wrap(~observation, ncol = 2, scales = "free_y") +
        theme_bw()
    } else {
      gplt = ggplot(full_df, aes(x = Year)) +
        geom_point(aes(y = Observed, col = "Observed")) +
        geom_line(aes(y = Predicted, col = "Predicted"), linewidth= 1.1, linetype = "dashed") +
        geom_errorbar(aes(ymin=L_CI, ymax=U_CI, col = "Observed"), width=.2, position=position_dodge(.9)) +
        guides( linewidth = "none", linetype = "none") +
        labs(y = "Index", col = "", linetype = "") +
        facet_wrap(~observation, ncol = 2, scales = "free_y") +
        theme_bw()
    }

  } else {
    gplt = ggplot(full_df, aes(x = Year)) +
      geom_point(aes(y = Observed, col = "Observed")) +
      geom_line(aes(y = Predicted, col = "Predicted"), linewidth= 1.1, linetype = "dashed") +
      geom_errorbar(aes(ymin=L_CI, ymax=U_CI, col = "Observed"), width=.2, position=position_dodge(.9)) +
      guides( linewidth = "none", linetype = "none") +
      labs(y = "Index", col = "", linetype = "") +
      facet_wrap(Survey~Region, ncol = 2, scales = "free_y") +
      theme_bw()
  }

  return(gplt)
}


#'
#' get_index
#' @param MLE_report a list that is output from obj$report() usually once an optimsation routine has been done.
#' @param region_key data.frame with colnames area and TMB_ndx for providing real region names to objects
#' @param survey_labels character vector for each of the n_surveys
#' @return ggplot2 object that will plot if an observation occurs in a year and region
#' @export
get_index = function(MLE_report, region_key = NULL, survey_labels = NULL) {
  years = MLE_report$years
  regions = 1:MLE_report$n_regions
  full_df = NULL
  if(MLE_report$model_type == 0) {
    if(sum(MLE_report$srv_jap_fishery_ll_bio_indicator) > 0) {
      jap_fishery_cpue = data.frame(Year = MLE_report$years[MLE_report$srv_jap_fishery_ll_bio_indicator == 1],SE = MLE_report$se_jap_fishery_ll_bio, Observed = MLE_report$obs_jap_fishery_ll_bio, Predicted = MLE_report$pred_jap_fishery_ll_bio, observation = "Japanese CPUE")
      ## convert the SE of an estimator to a standard deviation that is the
      ## right scale for the lognormal distribution
      ## first calculate CV = sigma/mean then pass this to the log_sigma function
      if(MLE_report$srv_jap_fishery_ll_bio_likelihood == 0)
        jap_fishery_cpue$SE = log_sigma(jap_fishery_cpue$SE / jap_fishery_cpue$Observed)
      full_df = rbind(full_df, jap_fishery_cpue)
    }
    if(sum(MLE_report$ll_cpue_indicator) > 0) {
      fixed_cpue = data.frame(Year = MLE_report$years[MLE_report$ll_cpue_indicator == 1],SE = MLE_report$se_ll_cpue, Observed = MLE_report$obs_ll_cpue, Predicted = MLE_report$pred_ll_cpue, observation = "Fixed gear CPUE")
      if(MLE_report$ll_cpue_likelihood == 0)
        fixed_cpue$SE = log_sigma(fixed_cpue$SE / fixed_cpue$Observed)
      full_df = rbind(full_df, fixed_cpue)
    }
    if(sum(MLE_report$srv_jap_ll_bio_indicator) > 0) {
      jap_srv_ll = data.frame(Year = MLE_report$years[MLE_report$srv_jap_ll_bio_indicator == 1],SE = MLE_report$se_jap_ll_bio, Observed = MLE_report$obs_jap_ll_bio, Predicted = MLE_report$pred_jap_ll_bio, observation = "Japanese LL survey")
      if(MLE_report$srv_jap_ll_bio_likelihood == 0)
        jap_srv_ll$SE = log_sigma(jap_srv_ll$SE / jap_srv_ll$Observed)
      full_df = rbind(full_df, jap_srv_ll)
    }
    if(sum(MLE_report$srv_nmfs_trwl_bio_indicator) > 0) {
      srv_nmfs_trwl = data.frame(Year = MLE_report$years[MLE_report$srv_nmfs_trwl_bio_indicator == 1],SE = MLE_report$se_nmfs_trwl_bio, Observed = MLE_report$obs_nmfs_trwl_bio, Predicted = MLE_report$pred_nmfs_trwl_bio, observation = "NMFS trawl survey")
      if(MLE_report$srv_nmfs_trwl_bio_likelihood == 0)
        srv_nmfs_trwl$SE = log_sigma(srv_nmfs_trwl$SE / srv_nmfs_trwl$Observed)
      full_df = rbind(full_df, srv_nmfs_trwl)
    }
    if(sum(MLE_report$srv_dom_ll_bio_indicator) > 0) {
      srv_dom_ll = data.frame(Year = MLE_report$years[MLE_report$srv_dom_ll_bio_indicator == 1],SE = MLE_report$se_dom_ll_bio, Observed = MLE_report$obs_dom_ll_bio, Predicted = MLE_report$pred_dom_ll_bio, observation = "Domestic LL survey")
      if(MLE_report$srv_dom_ll_bio_likelihood == 0)
        srv_dom_ll$SE = log_sigma(srv_dom_ll$SE / srv_dom_ll$Observed)
      full_df = rbind(full_df, srv_dom_ll)
    }
    ## combine all datasets
    full_df$Region = 1
  } else {
    surveys = 1:MLE_report$n_surveys
    dimnames(MLE_report$obs_srv_bio) = dimnames(MLE_report$pred_srv_bio) =   dimnames(MLE_report$obs_srv_se) = list(regions, years, surveys)
    MLE_report$obs_srv_bio[MLE_report$srv_bio_indicator == 0] = NA
    MLE_report$obs_srv_se[MLE_report$srv_bio_indicator == 0] = NA
    MLE_report$pred_srv_bio[MLE_report$srv_bio_indicator == 0] = NA

    index_obs = reshape2::melt(MLE_report$obs_srv_bio)
    index_se = reshape2::melt(MLE_report$obs_srv_se)
    index_fit = reshape2::melt(MLE_report$pred_srv_bio)
    colnames(index_obs) = c("Region", "Year", "Survey", "Observed")
    colnames(index_fit) = c("Region", "Year", "Survey", "Predicted")
    colnames(index_se) = c("Region", "Year", "Survey", "SE")
    index_obs$Predicted = index_fit$Predicted
    ## convert the SE of an estimator to a standard deviation that is the
    ## right scale for the lognormal distribution
    ## first calculate CV = sigma/mean then pass this to the log_sigma function
    for(s_ndx in 1:MLE_report$n_surveys) {
      if(MLE_report$srv_bio_likelihood[s_ndx] == 0) {
        index_se$SE[which(index_se$Survey == s_ndx & index_se$SE != 0)] = log_sigma(index_se$SE[which(index_se$Survey == s_ndx & index_se$SE != 0)] / index_obs$Observed[which(index_se$Survey == s_ndx & index_se$SE != 0)])
      }
    }
    if(!is.null(survey_labels)) {
      index_obs$Survey = survey_labels[index_obs$Survey]
    } else {
      surveys = paste0("Survey ", 1:MLE_report$n_surveys)
      index_obs$Survey = surveys[index_obs$Survey]
    }

    index_obs$SE = index_se$SE
    full_df = index_obs
  }
  CIs = lognormal_CI(full_df$Observed, sigma = full_df$SE, CI = 0.95)
  full_df$U_CI = CIs$upper
  full_df$L_CI = CIs$lower

  full_df$Pearsons_residuals = (log(full_df$Observed/full_df$Predicted) + 0.5*full_df$SE^2)/full_df$SE

  if(is.null(region_key)) {
    full_df$Region = paste0("Region ", full_df$Region)
  } else {
    full_df$Region = region_key$area[match(full_df$Region, (region_key$TMB_ndx + 1))]
  }
  full_df = full_df %>% dplyr::filter(!is.na(Observed))
  ## change name from label to
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
  age_comp_df = get_AF(MLE_report = MLE_report, region_key = region_key, observation = "all", sex = "both")
  len_comp_df = get_LF(MLE_report = MLE_report, region_key = region_key, observation = "all", sex = "both")
  ## drop NA observed
  len_comp_df = len_comp_df %>% filter(!is.na(Observed))
  age_comp_df = age_comp_df %>% filter(!is.na(Observed))
  ## get effective sample size
  len_comp_df = len_comp_df %>% group_by(Year, Region, observation) %>% mutate(Nassumed = sum(Observed), O_prop = Observed / Nassumed, P_prop = Predicted / Nassumed,
                                                                         O_length = O_prop * Length, P_length = P_prop * Length,  P_length_sq = P_prop * Length^2)
  age_comp_df = age_comp_df %>% group_by(Year, Region, observation) %>% mutate(Nassumed = sum(Observed), O_prop = Observed / Nassumed, P_prop = Predicted / Nassumed,
                                                                         O_age = O_prop * Age, P_age = P_prop * Age,  P_age_sq = P_prop * Age^2)
  ## summarise for each obs and year
  mean_len_df = len_comp_df %>% group_by(Year, Region, observation) %>% summarise(O_mean_length = sum(O_length),P_mean_length = sum(P_length), P_mean_length_sq = sum(P_length_sq), Nassumed = mean(Nassumed)) %>%
    mutate(stand_mean_length = sqrt(P_mean_length_sq - P_mean_length^2), resid_mean_length = O_mean_length - P_mean_length)
  mean_age_df = age_comp_df %>% group_by(Year, Region, observation) %>% summarise(O_mean_age = sum(O_age),P_mean_age = sum(P_age), P_mean_age_sq = sum(P_age_sq), Nassumed = mean(Nassumed)) %>%
    mutate(stand_mean_age = sqrt(P_mean_age_sq - P_mean_age^2),resid_mean_age = O_mean_age - P_mean_age)
  ## get the multiplier over all years for each observation
  length_multipliers = mean_len_df %>% group_by(Region, observation) %>% summarise(multiplier = 1/var(resid_mean_length * sqrt(Nassumed)/stand_mean_length, na.rm = T))
  age_multipliers = mean_age_df %>% group_by(Region, observation) %>% summarise(multiplier = 1/var(resid_mean_age * sqrt(Nassumed)/stand_mean_age, na.rm = T))
  ## some regions may only have a one sample which will cause NA's
  ## impute a weight with the mean weight
  length_multipliers$multiplier[is.na(length_multipliers$multiplier)] = mean(length_multipliers$multiplier, na.rm = T)
  age_multipliers$multiplier[is.na(age_multipliers$multiplier)] = mean(age_multipliers$multiplier, na.rm = T)

  return(list(length_multipliers = length_multipliers, age_multipliers = age_multipliers, mean_age_df = mean_age_df, mean_len_df = mean_len_df))
}


#'
#'
#' simulate_observations
#' @details Simulate observations conditional on MLE estimates
#' @param obj A TMB object that has been build using `TMB::MakeADFun`
#' @param n_sims an integer specifying how many simulated data sets you want
#' @param sd_report if you have already run `TMB::sdreport` then you can pass this function it. Which can save a time
#' @param include_param_uncertainty boolean whether to simulate parameters from a multivariate normal distribution for each simulation
#' @param region_key data.frame with colnames area and TMB_ndx for providing real region names to objects
#' @param survey_labels character vector for each of the n_surveys
#' @return named list containing simualted observations for all key obsevations
#' @export

simulate_observations <- function(obj, n_sims = 200, sd_report = NULL, include_param_uncertainty = F, region_key = NULL, survey_labels = NULL) {
  fixed_effect_pars = get_tmb_fixed_effects(obj)
  all_pars = obj$env$last.par.best
  if(include_param_uncertainty) {
    if(is.null(sd_report))
      sd_report = sdreport(obj, getJointPrecision = T)
    sim_pars =  MASS::mvrnorm(n = n_sims, mu = fixed_effect_pars, Sigma = sd_report$cov.fixed)
  }
  mle_rep = obj$report(all_pars)
  ## get the real observed data and save them in the following datasets
  obs_index_df = get_index(MLE_report = mle_rep, region_key = region_key, survey_labels = survey_labels)
  obs_srv_AF_df = get_AF(MLE_report = mle_rep, observation = "srv", sex = "both", region_key = region_key, survey_labels = survey_labels)
  obs_fixed_AF_df = get_AF(MLE_report = mle_rep, observation = "fixed", region_key = region_key)
  obs_fixed_LF_df = get_LF(MLE_report = mle_rep, observation = "fixed", region_key = region_key)
  obs_trwl_LF_df = get_LF(MLE_report = mle_rep, observation = "trwl", region_key = region_key)
  if(sum(mle_rep$tag_recovery_indicator_by_year) > 0)
    obs_tag_data_df = get_tag_recovery_obs_fitted_values(MLE_report = mle_rep, region_key = region_key)


  sim_srv_bio = sim_srv_AF = sim_fixed_AF = sim_fixed_LF = sim_trwl_LF = sim_tag_recovery = NULL
  for(sim_iter in 1:n_sims) {
    if(sim_iter %% 50 == 0)
      cat("simulation iteration: ", sim_iter, "\n")
    this_sim = NULL
    ## simualte
    if(include_param_uncertainty) {
      this_sim = obj$simulate(par = sim_pars[sim_iter,], complete = T)
    } else {
      this_sim = obj$simulate(par = all_pars, complete = T)
    }
    ## store sim obs
    # survey biomass
    index_df = get_index(this_sim, region_key = region_key, survey_labels = survey_labels) %>% rename(Simulated = Observed)
    index_df$Observed = obs_index_df$Observed
    index_df$sim = sim_iter
    sim_srv_bio = rbind(sim_srv_bio, index_df)
    # survey AF
    srv_AF = get_AF(MLE_report = this_sim, observation = "srv", region_key = region_key, survey_labels = survey_labels) %>% rename(Simulated = Observed)
    srv_AF$Observed = obs_srv_AF_df$Observed
    srv_AF$sim = sim_iter
    sim_srv_AF = rbind(sim_srv_AF, srv_AF)
    # Fixed AF
    fixed_AF = get_AF(MLE_report = this_sim, observation = "fixed", region_key = region_key) %>% rename(Simulated = Observed)
    fixed_AF$Observed = obs_fixed_AF_df$Observed
    fixed_AF$sim = sim_iter
    sim_fixed_AF = rbind(sim_fixed_AF, fixed_AF)
    # Fixed LF
    fixed_LF = get_LF(MLE_report = this_sim, observation = "fixed", region_key = region_key) %>% rename(Simulated = Observed)
    fixed_LF$Observed = obs_fixed_LF_df$Observed
    fixed_LF$sim = sim_iter
    sim_fixed_LF = rbind(sim_fixed_LF, fixed_LF)
    # Trawl LF
    trwl_LF = get_LF(MLE_report = this_sim, observation = "trwl", region_key = region_key) %>% rename(Simulated = Observed)
    trwl_LF$Observed = obs_trwl_LF_df$Observed
    trwl_LF$sim = sim_iter
    sim_trwl_LF = rbind(sim_trwl_LF, trwl_LF)
    # Tag data
    if(sum(mle_rep$tag_recovery_indicator_by_year) > 0) {
      tag_data = get_tag_recovery_obs_fitted_values(MLE_report = this_sim, region_key = region_key) %>% rename(Simulated = observed, Predicted = predicted)
      tag_data$Observed  = obs_tag_data_df$observed
      tag_data$sim = sim_iter
      sim_tag_recovery = rbind(sim_tag_recovery, tag_data)
    }
  }
  return(list(sim_tag_recovery = sim_tag_recovery, sim_trwl_LF = sim_trwl_LF, sim_fixed_LF = sim_fixed_LF, sim_fixed_AF = sim_fixed_AF, sim_srv_AF = sim_srv_AF, sim_srv_bio = sim_srv_bio))
}


#' calculate_simulated_residuals
#' @details take simulated observations created from `simulate_observations` and creating simulated scaled residuals using the DHARMa package \insertCite{dharma}{SpatialSablefishAssessment}
#' @param type string defining type
#' @param sim_ob a data frame which should be one of the elements that are created from `simulate_observations`
#' @param type string defining type
#' \itemize{
#'   \item `abundance`: abundance
#'   \item `AF`: Age frequency
#'   \item `LF`: Length frequency
#'   \item `tag-recovery`: tag recovery observations
#' }
#' @export
#' @return data.frame of scaled simulated residuals for each release and recovery event
#' @references
#' \insertAllCited{}
calculate_simulated_residuals <- function(sim_ob, type = "abundance") {
  if(!type %in% c("abundance", "AF", "LF", "tag-recovery"))
    stop('type needs to be one of the following c("abundance", "AF", "LF", "tag-recovery")')
  full_simualted_resids = NULL
  if(type == "abundance") {
    regions = unique(sim_ob$Region)
    surveys = unique(sim_ob$Survey)
    for(r in 1:length(regions)) {
      for(s in 1:length(surveys)) {
        this_sim_vals = sim_ob %>% filter(Region == regions[r], Survey == surveys[s]) %>%
          pivot_wider(id_cols = Year, names_from = sim, values_from = Simulated) %>%
          dplyr::select(!Year)
        this_sim_obs = sim_ob %>% filter(Region == regions[r], sim == 1, Survey == surveys[s]) %>% dplyr::select(Observed)
        this_pred = sim_ob %>% filter(Region == regions[r], sim == 1, Survey == surveys[s]) %>% dplyr::select(Predicted)
        Years = sim_ob %>% filter(Region == regions[r], sim == 1, Survey == surveys[s]) %>% dplyr::select(Year)

        this_dharma = suppressMessages(createDHARMa(simulatedResponse = as.matrix(this_sim_vals), observedResponse = this_sim_obs$Observed, fittedPredictedResponse  = this_pred$Predicted))
        this_scaled_df = data.frame(scaled_resids = this_dharma$scaledResiduals, qnorm_transformed_scaled_resids = qnorm(this_dharma$scaledResiduals), observed = this_dharma$observedResponse, mean_sim_vals = this_dharma$fittedPredictedResponse, Region = regions[r], Year =Years, Survey = surveys[s])
        full_simualted_resids = rbind(full_simualted_resids, this_scaled_df)
      }
    }
    full_simualted_resids$type = "Abundance"
  } else if (type == "LF") {
    regions = unique(sim_ob$Region)
    for(r in 1:length(regions)) {
      this_df = sim_ob %>% filter(Region == regions[r])
      years_this_region = unique(this_df$Year)
      for(y in 1:length(years_this_region)) {
        this_sim_vals =  this_df %>% filter(Year == years_this_region[y]) %>%
          pivot_wider(id_cols = S_Length, names_from = sim, values_from = Simulated) %>%
          ungroup() %>%
          dplyr::select(!S_Length)
        this_sim_obs = this_df %>% filter(Year == years_this_region[y], sim == 1) %>% ungroup() %>% dplyr::select(Observed)
        S_Lengths = this_df %>% ungroup() %>% filter(Year == years_this_region[y], sim == 1) %>% dplyr::select(S_Length)

        this_dharma = suppressMessages(createDHARMa(simulatedResponse = as.matrix(this_sim_vals), observedResponse = this_sim_obs$Observed, integerResponse = T))
        this_scaled_df = data.frame(scaled_resids = this_dharma$scaledResiduals, qnorm_transformed_scaled_resids = qnorm(this_dharma$scaledResiduals), observed = this_dharma$observedResponse, mean_sim_vals = this_dharma$fittedPredictedResponse, Region = regions[r], Year = years_this_region[y], S_Lengths = S_Lengths)
        full_simualted_resids = rbind(full_simualted_resids, this_scaled_df)
      }
    }
    full_simualted_resids$sex = substring(full_simualted_resids$S_Length, first = 1, last = 1)
    full_simualted_resids$length = as.numeric(substring(full_simualted_resids$S_Length, first = 3))
    full_simualted_resids$type = "LF"
  } else if (type == "AF") {
    regions = unique(sim_ob$Region)
    surveys = unique(sim_ob$Survey)
    for(s in 1:length(surveys)) {
      for(r in 1:length(regions)) {
        this_df = sim_ob %>% filter(Region == regions[r], Survey == surveys[s])
        years_this_region = unique(this_df$Year)

        for(y in 1:length(years_this_region)) {
          this_sim_vals =  this_df %>% filter(Year == years_this_region[y]) %>%
            pivot_wider(id_cols = S_Age, names_from = sim, values_from = Simulated) %>%
            ungroup() %>%
            dplyr::select(!S_Age)
          this_sim_obs = this_df %>% filter(Year == years_this_region[y], sim == 1) %>% ungroup() %>% dplyr::select(Observed)
          S_Ages = this_df %>% ungroup() %>% filter(Year == years_this_region[y], sim == 1) %>% dplyr::select(S_Age)

          this_dharma = suppressMessages(createDHARMa(simulatedResponse = as.matrix(this_sim_vals), observedResponse = this_sim_obs$Observed, integerResponse = T))
          this_scaled_df = data.frame(scaled_resids = this_dharma$scaledResiduals, qnorm_transformed_scaled_resids = qnorm(this_dharma$scaledResiduals), observed = this_dharma$observedResponse, mean_sim_vals = this_dharma$fittedPredictedResponse, Region = regions[r], Survey = surveys[s], Year = years_this_region[y], S_Ages = S_Ages)
          full_simualted_resids = rbind(full_simualted_resids, this_scaled_df)
        }
      }
    }
    full_simualted_resids$sex = substring(full_simualted_resids$S_Age, first = 1, last = 1)
    full_simualted_resids$age = as.numeric(substring(full_simualted_resids$S_Age, first = 3))
    full_simualted_resids$type = "AF"

  } else if(type == "tag-recovery") {
    rel_event = unique(sim_ob$release_event)
    sim_ob$recovery_event = paste0(sim_ob$recovery_year, "-", sim_ob$recovery_region)
    for(r in 1:length(rel_event)) {
      this_sim_vals =  sim_ob %>% filter(release_event == rel_event[r]) %>%
        pivot_wider(id_cols = recovery_event, names_from = sim, values_from = Simulated) %>%
        dplyr::select(!recovery_event)
      this_sim_obs = sim_ob %>% filter(release_event == rel_event[r], sim == 1) %>% ungroup() %>% dplyr::select(Observed)
      rec_event = sim_ob %>% ungroup() %>% filter(release_event == rel_event[r], sim == 1) %>% dplyr::select(recovery_event)

      if(length(this_sim_obs$Observed) < 3) {
        cat("skipping simulated obs for release event ", rel_event[r], ". The DHarama package wont calculate residauls when number of observations < 3. Found ", length(this_sim_obs$Observed), "\n")
      } else {
        this_dharma = suppressMessages(createDHARMa(simulatedResponse = as.matrix(this_sim_vals), observedResponse = this_sim_obs$Observed, integerResponse = T))
        this_scaled_df = data.frame(scaled_resids = this_dharma$scaledResiduals, qnorm_transformed_scaled_resids = qnorm(this_dharma$scaledResiduals), observed = this_dharma$observedResponse, mean_sim_vals = this_dharma$fittedPredictedResponse, recovery_event = rec_event, release_event = rel_event[r])
        full_simualted_resids = rbind(full_simualted_resids, this_scaled_df)
      }
    }
    full_simualted_resids$recovery_region = Reduce(c, lapply(strsplit(full_simualted_resids$recovery_event, split = "-"), FUN = function(x){x[2]}))
    full_simualted_resids$recovery_year = Reduce(c, lapply(strsplit(full_simualted_resids$recovery_event, split = "-"), FUN = function(x){x[1]}))
    full_simualted_resids$release_region = Reduce(c, lapply(strsplit(full_simualted_resids$release_event, split = "-"), FUN = function(x){x[2]}))
    full_simualted_resids$release_year = Reduce(c, lapply(strsplit(full_simualted_resids$release_event, split = "-"), FUN = function(x){x[1]}))
    full_simualted_resids$type = "Tag"
  }
  return(full_simualted_resids);
}


#' calculate_pearson_residuals
#' @details Take a model fit and calculate Pearson resiudals for all observations.
#' @param MLE_report a list that is output from obj$report() usually once an optimsation routine has been done.
#' @param region_key data.frame with colnames area and TMB_ndx for providing real region names to objects
#' @export
#' @return data.frame of pearson residuals for each observation
calculate_pearson_residuals <- function(MLE_report, region_key = NULL) {
  return(NULL)
}

#' summarise_AF_quant_resids
#'
#' @param AF_sim_resids an AF object created from `calculate_simulated_residuals`
#' @param sex character specifying which sex to plot
#' \itemize{
#'   \item `M`: Male
#'   \item `F`: Female
#' }
#' @param obs_label a character specifying the title possible observation label
#' @export
#' @return a joint plot
summarise_AF_quant_resids <- function(AF_sim_resids, sex = "M", obs_label = "") {
  title_label = paste0(ifelse(sex == "M", "Male", "Female"), " ", obs_label)
  AF_sim_resids$year_class  = AF_sim_resids$Year - AF_sim_resids$age
  yr_plt = ggplot(AF_sim_resids %>% filter(sex == sex), aes(x = Year, group = Year, y = (qnorm_transformed_scaled_resids))) +
    geom_boxplot() +
    labs(x = "Year", y = "Simualted quantile residuals") +
    ggtitle(title_label) +
    geom_hline(yintercept = 0, col = "red", linetype = "dashed") +
    facet_wrap(~Region, ncol = 1) +
    scale_x_discrete(breaks = every_nth(n = 10)) +
    theme_bw()
  yr_class_plt = ggplot(AF_sim_resids %>% filter(sex == sex), aes(x = year_class, group = year_class, y = (qnorm_transformed_scaled_resids))) +
    geom_boxplot() +
    labs(x = "Year class", y = "") +
    ggtitle("") +
    geom_hline(yintercept = 0, col = "red", linetype = "dashed") +
    facet_wrap(~Region, ncol = 1) +
    scale_x_discrete(breaks = every_nth(n = 10)) +
    theme_bw()
  age_plt = ggplot(AF_sim_resids %>% filter(sex == sex), aes(x = age, group = age, y = (qnorm_transformed_scaled_resids))) +
    geom_boxplot() +
    labs(x = "Age", y = "") +
    ggtitle("") +
    facet_wrap(~Region, ncol = 1) +
    geom_hline(yintercept = 0, col = "red", linetype = "dashed") +
    #scale_x_discrete(breaks = every_nth(n = 10)) +
    theme_bw()

  joint_tag_plt = ggarrange(yr_plt, yr_class_plt, age_plt,
                            ncol = 3)

  return(joint_tag_plt)
}


#' summarise_LF_quant_resids
#'
#' @param LF_sim_resids an LF object created from `calculate_simulated_residuals`
#' @param sex character specifying which sex to plot
#' \itemize{
#'   \item `M`: Male
#'   \item `F`: Female
#' }
#' @param obs_label a character specifying the title possible observation label
#' @export
#' @return a joint plot
summarise_LF_quant_resids <- function(LF_sim_resids, sex = "M", obs_label = "") {
  title_label = paste0(ifelse(sex == "M", "Male", "Female"), " ", obs_label)
  yr_plt = ggplot(LF_sim_resids %>% filter(sex == sex), aes(x = Year, group = Year, y = (qnorm_transformed_scaled_resids))) +
    geom_boxplot() +
    labs(x = "Year", y = "Simualted quantile residuals") +
    ggtitle(title_label) +
    geom_hline(yintercept = 0, col = "red", linetype = "dashed") +
    facet_wrap(~Region, ncol = 1) +
    scale_x_discrete(breaks = every_nth(n = 10)) +
    theme_bw()
  len_plt = ggplot(LF_sim_resids %>% filter(sex == sex), aes(x = length, group = length, y = (qnorm_transformed_scaled_resids))) +
    geom_boxplot() +
    labs(x = "Length bin", y = "") +
    ggtitle("") +
    facet_wrap(~Region, ncol = 1) +
    geom_hline(yintercept = 0, col = "red", linetype = "dashed") +
    scale_x_discrete(breaks = every_nth(n = 10)) +
    theme_bw()

  joint_tag_plt = ggarrange(yr_plt, len_plt,
                            ncol = 3)

  return(joint_tag_plt)
}

#' summarise_tag_quant_resids
#'
#' @param tag_sim_resids an LF object created from `calculate_simulated_residuals`
#' @param recovery_year boolean whether to plot by recovery year = T or release year = F
#' @export
#' @return a joint plot
summarise_tag_quant_resids <- function(tag_sim_resids, recovery_year = T) {
  if(recovery_year) {
    rel_reg_rec_yr_plt = ggplot(tag_sim_resids, aes(x = recovery_year, group = recovery_year, y = (qnorm_transformed_scaled_resids))) +
      geom_boxplot() +
      labs(x = "Recovery Year", y = "") +
      ggtitle("Release Region") +
      geom_hline(yintercept = 0, col = "red", linetype = "dashed") +
      facet_wrap(~release_region, ncol = 1) +
      scale_x_discrete(breaks = every_nth(n = 10)) +
      theme_bw()
    rec_reg_rec_yr_plt = ggplot(tag_sim_resids, aes(x = recovery_year, group = recovery_year, y = (qnorm_transformed_scaled_resids))) +
      geom_boxplot() +
      labs(x = "Recovery Year", y = "") +
      ggtitle("Recovery Region") +
      facet_wrap(~recovery_region, ncol = 1) +
      geom_hline(yintercept = 0, col = "red", linetype = "dashed") +
      scale_x_discrete(breaks = every_nth(n = 10)) +
      theme_bw()
    joint_tag_plt = ggarrange(rel_reg_rec_yr_plt, rec_reg_rec_yr_plt,
                              ncol = 2)
  } else {
    rel_reg_rel_yr_plt = ggplot(tag_sim_resids, aes(x = release_year, group = release_year, y = (qnorm_transformed_scaled_resids))) +
      geom_boxplot() +
      labs(x = "Release Year", y = "Simualted quantile residuals") +
      ggtitle("Release Region") +
      geom_hline(yintercept = 0, col = "red", linetype = "dashed") +
      facet_wrap(~release_region, ncol = 1) +
      scale_x_discrete(breaks = every_nth(n = 10)) +
      theme_bw()
    rec_reg_rel_yr_plt = ggplot(tag_sim_resids, aes(x = release_year, group = release_year, y = (qnorm_transformed_scaled_resids))) +
      geom_boxplot() +
      labs(x = "Release Year", y = "") +
      ggtitle("Recovery Region") +
      facet_wrap(~recovery_region, ncol = 1) +
      geom_hline(yintercept = 0, col = "red", linetype = "dashed") +
      scale_x_discrete(breaks = every_nth(n = 10)) +
      theme_bw()
    joint_tag_plt = ggarrange(rel_reg_rel_yr_plt, rec_reg_rel_yr_plt,ncol = 2)
  }
  return(joint_tag_plt)
}

#' get_comp_sample_size
#' @details get compositional input sample sizes and effective sample sizes
#' @param MLE_report a list that is output from obj$report() usually once an optimsation routine has been done.
#' @param data list of input data
#' @param region_key data.frame with colnames area and TMB_ndx for providing real region names to objects
#' @param survey_labels character vector for each of the n_surveys
#' @return data.frame
#' @export
get_comp_sample_size <- function(MLE_report, data, region_key = NULL, survey_labels = NULL) {
  years = MLE_report$years
  regions = 1:MLE_report$n_regions
  if(is.null(region_key)) {
    regions = paste0("Region ", regions)
  } else {
    regions = region_key$area[match(regions, (region_key$TMB_ndx + 1))]
  }
  ages = MLE_report$ages
  length_bins = MLE_report$length_bins
  sample_sizes = NULL
  if(MLE_report$model_type == 0) {
    LFs = get_LF(MLE_report = MLE_report, observation = "all", subset_years = NULL, region_key = region_key, sex = "both")
    this_sample_size = LFs %>% group_by(Year, Region, observation, Sex) %>% summarise(N_eff = sum(Observed), N_input = sum(Observed))
    this_sample_size$observation = paste0(this_sample_size$observation, " LF")
    sample_sizes = rbind(sample_sizes, this_sample_size)
    AFs = get_AF(MLE_report = MLE_report, observation = "all", subset_years = NULL, region_key = region_key, sex = "both")
    this_sample_size = LFs %>% group_by(Year, Region, observation, Sex) %>% summarise(N_eff = sum(Observed), N_input = sum(Observed))
    this_sample_size$observation = paste0(this_sample_size$observation, " AF")
    sample_sizes = rbind(sample_sizes, this_sample_size)

  } else {
    ## check likelihoods
    ## fixed catch at age
    dimnames(data$obs_fixed_catchatage)= list(c(paste0("M_",ages), paste0("F_",ages)), regions, years)
    molten_fixed_catchatage = reshape2::melt(data$obs_fixed_catchatage)
    colnames(molten_fixed_catchatage) = c("s_age", "Region", "Year", "value")
    this_sample_size = molten_fixed_catchatage %>% group_by(Region, Year) %>% summarise(N_input = sum(value), N_eff = N_input)
    if(data$fixed_catchatage_comp_likelihood == 1)
      this_sample_size$N_eff = (1 + this_sample_size$N_input * MLE_report$theta_fixed_catchatage ) / (1 + MLE_report$theta_fixed_catchatage)
    this_sample_size$observation = "Fixed Catch-at-age"
    sample_sizes = rbind(sample_sizes, this_sample_size)

    ## fixed catch at length
    dimnames(data$obs_fixed_catchatlgth)= list(c(paste0("M_",length_bins), paste0("F_",length_bins)), regions, years)
    molten_fixed_catchatlgth = reshape2::melt(data$obs_fixed_catchatlgth)
    colnames(molten_fixed_catchatlgth) = c("s_lgth", "Region", "Year", "value")
    this_sample_size = molten_fixed_catchatlgth %>% group_by(Region, Year) %>% summarise(N_input = sum(value), N_eff = N_input)
    if(data$fixed_catchatlgth_comp_likelihood == 1)
      this_sample_size$N_eff = (1 + this_sample_size$N_input * MLE_report$theta_fixed_catchatlgth ) / (1 + MLE_report$theta_fixed_catchatlgth)
    this_sample_size$observation = "Fixed Catch-at-length"
    sample_sizes = rbind(sample_sizes, this_sample_size)

    ## Trawl catch at length
    dimnames(data$obs_trwl_catchatlgth)= list(c(paste0("M_",length_bins), paste0("F_",length_bins)), regions, years)
    molten_trwl_catchatlgth = reshape2::melt(data$obs_trwl_catchatlgth)
    colnames(molten_trwl_catchatlgth) = c("s_lgth", "Region", "Year", "value")
    this_sample_size = molten_trwl_catchatlgth %>% group_by(Region, Year) %>% summarise(N_input = sum(value), N_eff = N_input)
    if(data$trwl_catchatlgth_comp_likelihood == 1)
      this_sample_size$N_eff = (1 + this_sample_size$N_input * MLE_report$theta_trwl_catchatlgth) / (1 + MLE_report$theta_trwl_catchatlgth)
    this_sample_size$observation = "Trawl Catch-at-length"
    sample_sizes = rbind(sample_sizes, this_sample_size)
    ## Survey catch at age
    surveys = paste0("Survey ", 1:MLE_report$n_surveys)
    if(!is.null(survey_labels))
      surveys = survey_labels
    dimnames(data$obs_srv_catchatage)= list(c(paste0("M_",ages), paste0("F_",ages)), regions, years, surveys)
    molten_srv_catchatage = reshape2::melt(data$obs_srv_catchatage)
    colnames(molten_srv_catchatage) = c("s_age", "Region", "Year", "Survey" ,"value")
    this_sample_size = molten_srv_catchatage %>% group_by(Region, Year, Survey) %>% summarise(N_input = sum(value), N_eff = N_input)
    if(any(data$srv_catchatage_comp_likelihood == 1))
      this_sample_size$N_eff = (1 + this_sample_size$N_input * MLE_report$theta_srv_catchatage ) / (1 + MLE_report$theta_srv_catchatage)
    this_sample_size$observation = paste0("Survey Catch-at-age-", this_sample_size$Survey)
    sample_sizes = rbind(sample_sizes, this_sample_size)
  }
  return(sample_sizes)
}

#' plot_comp_sample_size
#' @details plot compositional input sample sizes or effective sample sizes
#' @param MLE_report a list that is output from obj$report() usually once an optimsation routine has been done.
#' @param region_key data.frame with colnames area and TMB_ndx for providing real region names to objects
#' @param data list of input data
#' @param survey_labels character vector for each of the n_surveys
#' @param plot_n_eff boolean if T plot effective sample size else plot input sample size
#' @return data.frame
#' @export
plot_comp_sample_size <- function(MLE_report, data, region_key = NULL, survey_labels = NULL, plot_n_eff = T) {
  sample_size_df = get_comp_sample_size(MLE_report, data, region_key, survey_labels)
  sample_size_df = sample_size_df %>% filter(N_input > 0)
  #sample_size_df$N_input =  dplyr::na_if(sample_size_df$N_input, 0)
  #sample_size_df$N_eff =  dplyr::na_if(sample_size_df$N_eff, 0)
  gplt = NULL
  if(MLE_report$model_type == 0) {
    if(plot_n_eff) {
      gplt <- ggplot(sample_size_df, aes(x = Year, y = N_eff, col = observation)) +
        geom_line(linewidth = 0.8, linetype= "dashed") +
        geom_point(size = 1.6) +
        theme_bw() +
        labs(x = "Year", y = "Effective sample size") +
        theme(legend.position = "none") +
        facet_wrap(~observation)
    } else {
      gplt <- ggplot(sample_size_df, aes(x = Year, y = N_input, col = observation)) +
        geom_line(linewidth = 0.8, linetype= "dashed") +
        geom_point(size = 1.6) +
        theme_bw() +
        labs(x = "Year", y = "Input sample size") +
        theme(legend.position = "none") +
        facet_wrap(~observation)
    }
  } else {
    if(plot_n_eff) {
      gplt <- ggplot(sample_size_df, aes(x = Year, y = N_eff, col = observation)) +
        geom_line(linewidth = 1.1) +
        geom_point(size = 1.1) +
        theme_bw() +
        labs(x = "Year", y = "Effective sample size") +
        facet_wrap(~Region)
    } else {
      gplt <- ggplot(sample_size_df, aes(x = Year, y = N_input, col = observation)) +
        geom_line(linewidth = 1.1) +
        geom_point(size = 1.1) +
        theme_bw() +
        labs(x = "Year", y = "Input sample size") +
        facet_wrap(~Region)
    }
  }

  return(gplt)
}
