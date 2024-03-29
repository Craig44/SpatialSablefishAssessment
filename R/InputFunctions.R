
#'
#' plot_input_catches
#' @param data list that is passed to the MakeADfun for the TMB model
#' @param region_key data.frame with colnames area and TMB_ndx for providing real region names to objects
#' @return ggplot2 object that will plot if an observation occurs in a year and region
#' @export
plot_input_catches = function(data, region_key = NULL) {
  years = data$years
  if(data$model == "Assessment") {
    regions = 1:1
    data$fixed_fishery_catch = matrix(data$ll_fishery_catch, nrow = 1)
    data$trwl_fishery_catch = matrix(data$trwl_fishery_catch, nrow = 1)
  } else {
    regions = 1:data$n_regions
  }
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
    theme_bw() +
    theme(axis.text = element_text(size = 13),
          axis.title = element_text(size = 13))
  return(gplt)
}
#'
#' plot_age_error_matrix
#' @param data list that is passed to the MakeADfun for the TMB model
#' @param fill_text add value in ggplot
#' @return ggplot2 object that will plot the ageing error matrix
#' @export
plot_age_error_matrix <- function(data, fill_text = F) {
  ageing_error_mat = data$ageing_error_matrix
  dimnames(ageing_error_mat) = list(data$ages, data$ages)
  if(any(abs(rowSums(ageing_error_mat) - 1) > 0.01))
    stop("ageing_error_mat expected to have rows that sum close to 1. This is not the case")
  molten_age_error = reshape2::melt(ageing_error_mat)
  colnames(molten_age_error) = c("Observed age", "Possible ages", "Proportion")
  gplt = ggplot(molten_age_error, aes(y = `Observed age`, x = `Possible ages`, fill = Proportion)) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = "red") +
    #geom_text(aes(y = `Observed age`, x = `Possible ages`, label = round(Proportion,2)), color = "black", size = 4) +
    theme_bw()
  if(fill_text)
    gplt = gplt + geom_text(aes(y = `Observed age`, x = `Possible ages`, label = round(Proportion,2)), color = "black", size = 2)
  return(gplt)
}


#'
#' get_input_observations
#' @param data list that is passed to the MakeADfun for the TMB model
#' @param region_key data.frame with colnames area and TMB_ndx for providing real region names to objects
#' @param survey_labels character vector for each of the n_surveys
#' @return data.frame when observation occurs in a year and region
#' @importFrom dplyr bind_rows
#' @export
get_input_observations = function(data, region_key = NULL, survey_labels = NULL) {
  years = data$years
  full_df = NULL
  if(data$model == "Assessment") {
    obs_indicator_df = data.frame(year = data$years,
               fixed_gear_age = data$ll_catchatage_indicator,
               fixed_gear_lgth = data$ll_catchatlgth_indicator,
               fixed_gear_bio = data$ll_cpue_indicator,
               trawl_gear_lgth = data$trwl_catchatlgth_indicator,
               longline_survey_age = data$srv_dom_ll_age_indicator,
               longline_survey_lgth = data$srv_dom_ll_lgth_indicator,
               longline_survey_bio = data$srv_dom_ll_bio_indicator,
               nmfs_survey_age = data$srv_nmfs_trwl_age_indicator,
               nmfs_survey_lgth = data$srv_nmfs_trwl_lgth_indicator,
               nmfs_survey_bio = data$srv_nmfs_trwl_bio_indicator,
               japanese_survey_age = data$srv_jap_ll_age_indicator,
               japanese_survey_lgth = data$srv_jap_ll_lgth_indicator,
               japanese_survey_bio = data$srv_jap_ll_bio_indicator,
               japanese_gear_bio = data$srv_jap_fishery_ll_bio_indicator,
               japanese_gear_lgth = data$srv_jap_fishery_ll_lgth_indicator
      )
    full_df = obs_indicator_df %>% pivot_longer(!year)

    colnames(full_df) = c("Year", "label", "indicator")
    full_df$indicator = ifelse(full_df$indicator == 0, NA, 1)

    region_lab = "Alaska"
    if(!is.null(region_key))
      region_lab = region_key$area[region_key$TMB_ndx + 1]
    full_df$Region = region_lab

    fixed_catch = trawl_catch = expand.grid(years, unique(full_df$Region))
    colnames(fixed_catch) = colnames(trawl_catch) = c("Year", "Region")
    fixed_catch$indicator = 1
    trawl_catch$indicator = 1
    fixed_catch$indicator = 1
    trawl_catch$indicator = 1
    trawl_catch$label = "trawl_gear_catch"
    fixed_catch$label = "fixed_gear_catch"
    full_df = bind_rows(full_df, trawl_catch, fixed_catch)

    ## Create labels
    full_df$type = Reduce(c, lapply(strsplit(full_df$label, split = "_"), function(x) {x[3]}))
    full_df$obs_type = Reduce(c, lapply(strsplit(full_df$label, split = "_"), function(x) {x[2]}))
    full_df$source = Reduce(c, lapply(strsplit(full_df$label, split = "_"), function(x) {x[1]}))
    ##
    full_df = full_df %>% mutate(type =
                         case_when(type == "age"  ~ "Age",
                                   type == "bio"  ~ "Abundance",
                                   type == "catch"  ~ "Catch",
                                   type == "lgth"  ~ "Length"),
                         obs_type = case_when(obs_type == "gear" ~ "Fishery",
                                            obs_type == "survey" ~ "Survey"),
                         source = case_when(source == "fixed" ~ "Fixed",
                                            source == "japanese" ~ "Japanese LL",
                                            source == "longline" ~ "Domestic LL",
                                            source == "nmfs" ~ "NMFS",
                                            source == "trawl" ~ "Trawl"
                                            )
                         )


  } else {
    regions = 1:data$n_regions
    dimnames(data$fixed_catchatage_indicator) = dimnames(data$fixed_catchatlgth_indicator) = dimnames(data$trwl_catchatlgth_indicator)= list(regions, years)
    surveys = paste0("Survey ", 1:data$n_surveys)
    if(!is.null(survey_labels))
      surveys = survey_labels
    dimnames(data$srv_catchatage_indicator) = list(regions, years, paste0(surveys, " AF"))
    dimnames(data$srv_bio_indicator) = list(regions, years, paste0(surveys, " bio"))
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
      tag_recovery_detailed$type = "tag"
    }
    colnames(fixed_catchatage) = colnames(fixed_catchatlgth) = colnames(trwl_catchatlgth) = c("Region", "Year", "indicator")
    colnames(srv_catchatage) = colnames(srv_bio) = c("Region", "Year", "label", "indicator")
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
      tag_release_df$type = "tag"
    }
    fixed_catchatlgth$label = "Fishery Fixed LF"
    fixed_catchatlgth$type = "length"

    fixed_catchatage$label = "Fishery Fixed AF"
    fixed_catchatage$type = "age"

    trwl_catchatlgth$label = "Fishery Trawl LF"
    trwl_catchatlgth$type = "length"

    srv_catchatage$type = "age"
    srv_bio$type = "abundance"
    ## combine
    full_df = rbind(fixed_catchatage, trwl_catchatlgth, srv_catchatage, srv_bio, fixed_catchatlgth, tag_recovery_detailed, tag_release_df)

    if(is.null(region_key)) {
      full_df$Region = paste0("Region ", full_df$Region)
    } else {
      full_df$Region = region_key$area[match(full_df$Region, (region_key$TMB_ndx + 1))]
    }

    full_df$indicator = ifelse(full_df$indicator == 0, NA, 1)

    fixed_catch = trawl_catch = expand.grid(years, unique(full_df$Region))
    colnames(fixed_catch) = colnames(trawl_catch) = c("Year", "Region")
    fixed_catch$indicator = 1
    trawl_catch$indicator = 1
    fixed_catch$indicator = 1
    trawl_catch$indicator = 1
    trawl_catch$label = "Trawl catch"
    fixed_catch$label = "Fixed catch"
    trawl_catch$type = "catch"
    fixed_catch$type = "catch"

    full_df = bind_rows(full_df, trawl_catch, fixed_catch)
    full_df$type = factor(full_df$type, levels = rev(c("catch", "age", "length", "abundance", "tag")), ordered = T)
  }
  return(full_df)
}

#'
#' plot_input_observations
#' @param data list that is passed to the MakeADfun for the TMB model
#' @param region_key data.frame with colnames area and TMB_ndx for providing real region names to objects
#' @param survey_labels character vector for each of the n_surveys
#' @return ggplot2 object that will plot if an observation occurs in a year and region
#' @export
plot_input_observations = function(data, region_key = NULL, survey_labels = NULL) {

  full_df = get_input_observations(data, region_key, survey_labels)
  gplt = NULL
  if(data$model == "Assessment") {
    full_df$temp_label = paste0(full_df$source, "-", full_df$obs_type)
    gplt = ggplot(full_df) +
      geom_point(aes(x = Year, y = temp_label, col = temp_label, size = indicator)) +
      guides(colour = "none", size = "none") +
      labs(y = "") +
      facet_wrap(~type, ncol = 1) +
      theme_bw() +
      theme(
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 14)
      )
  } else {

    gplt = ggplot(full_df) +
      geom_point(aes(x = Year, y = forcats::fct_reorder(label, as.integer(type)), col = label, size = indicator)) +
      guides(colour = "none", size = "none") +
      labs(y = "") +
      facet_wrap(~Region, ncol = 1) +
      theme_bw() +
      theme(
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 14)
      )
  }
  return(gplt)
}

#'
#' plot_input_timeblocks
#' @param data list that is passed to the MakeADfun for the TMB model
#' @param survey_labels character vector for each of the n_surveys
#' @return ggplot2 object visualising the time-blocks for selectivities and catchabilities
#' @export
plot_input_timeblocks = function(data, survey_labels = NULL) {

  if(data$model == "Assessment") {
    ## assessment model
    full_df = data.frame(Year = data$years,
                         fixed_sel = data$ll_sel_by_year_indicator +1,
                         fixed_q = data$ll_cpue_q_by_year_indicator + 1,
                         trwl_sel = data$trwl_sel_by_year_indicator +1,
                         srvdomll_sel = data$srv_dom_ll_sel_by_year_indicator +1,
                         srvdomll_q = data$srv_dom_ll_q_by_year_indicator +1,
                         srvjapll_sel = data$srv_jap_fishery_ll_sel_by_year_indicator +1,
                         srvjapll_q = data$srv_jap_fishery_ll_q_by_year_indicator +1,
                         srvnmfs_sel = data$srv_nmfs_trwl_sel_by_year_indicator +1,
                         srvnmfs_q = data$srv_nmfs_trwl_q_by_year_indicator +1,
                         japll_sel = data$srv_jap_fishery_ll_sel_by_year_indicator +1,
                         japll_q = data$srv_jap_fishery_ll_q_by_year_indicator +1)

    full_df_lng = full_df %>% pivot_longer(!Year)
    ## Create labels
    full_df_lng$label = Reduce(c, lapply(strsplit(full_df_lng$name, split = "_"), function(x) {x[1]}))
    full_df_lng$type = Reduce(c, lapply(strsplit(full_df_lng$name, split = "_"), function(x) {x[2]}))

    ##
    full_df_lng = full_df_lng %>% mutate(label =
                                   case_when(label == "fixed"  ~ "Fixed gear fishery",
                                             label == "trwl"  ~ "Trawl gear fishery",
                                             label == "japll"  ~ "Japanese LL survey",
                                             label == "srvdomll"  ~ "Domestic LL survey",
                                             label == "srvjapll"  ~ "Japanese LL fishery",
                                             label == "srvnmfs"  ~ "NMFS survey"),
                                   type = case_when(type == "q" ~ "Catchability",
                                                    type == "sel" ~ "Selectivity")
    )

    gplt = ggplot(full_df_lng) +
      geom_point(aes(x = Year, y = label, col = factor(value), fill = factor(value)), size = 3) +
      facet_wrap(~type, ncol = 1) +
      guides(size = "none", fill = "none") +
      labs(x = "Year", y = "", col = "Time block") +
      theme_bw() +
      theme(axis.text = element_text(size = 13),
            axis.title = element_text(size = 13))
  } else {
    ## Spatial tag-integrated model
    projyears = min(data$years):(max(data$years) + data$n_projections_years)
    full_df = data.frame(Year = projyears, fixed_sel = data$fixed_sel_by_year_indicator +1, trwl_sel = data$trwl_sel_by_year_indicator +1, Movement = data$movement_time_block_indicator +1)

    surveys = paste0("Survey ", 1:data$n_surveys)
    if(!is.null(survey_labels))
      surveys = survey_labels

    srv_sel = reshape2::melt(data$srv_sel_by_year_indicator +1)
    srv_q = reshape2::melt(data$srv_q_by_year_indicator +1)

    srv_q$survey = paste0(surveys[srv_q$Var2], "\nCatchability")
    srv_sel$survey = paste0(surveys[srv_sel$Var2], "\nSelectivity")
    colnames(srv_sel) = colnames(srv_q) = c("year_ndx", "survey_ndx", "time_block", "label")
    srv_q$Year = projyears[srv_q$year_ndx]
    srv_sel$Year = projyears[srv_sel$year_ndx]

    full_df_lng = full_df %>% pivot_longer(!Year)
    colnames(full_df_lng) = c("Year", "label", "time_block")
    full_df_lng = rbind(full_df_lng, srv_q %>% dplyr::select(colnames(full_df_lng)),  srv_sel %>% dplyr::select(colnames(full_df_lng)))

    full_df_lng$time_block = factor(full_df_lng$time_block)
    full_df_lng = full_df_lng %>% mutate(label = case_when(
      label == "fixed_sel" ~ "Fixed\nSelectivity",
      label == "trwl_sel" ~ "Trawl\nSelectivity",
      .default = label

    ))
    gplt = ggplot(full_df_lng) +
      geom_point(aes(x = Year, y = label, col = time_block, fill = time_block), size = 3) +
      guides(size = "none", fill = "none") +
      labs(x = "Year", y = "", col = "Time block") +
      theme_bw() +
      theme(axis.text = element_text(size = 13),
            axis.title = element_text(size = 13))
  }
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
