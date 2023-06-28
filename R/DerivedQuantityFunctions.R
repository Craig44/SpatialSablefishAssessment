#
#' plot_fishing_mortalities
#' @param MLE_report a list that is output from obj$report() usually once an optimsation routine has been done.
#' @param region_key data.frame with colnames area and TMB_ndx for providing real region names to objects
#' @return ggplot2 object that will plot if an observation occurs in a year and region
#' @export

plot_fishing_mortalities = function(MLE_report, region_key = NULL) {
  F_df = get_fishing_mortalities(MLE_report = MLE_report, region_key = region_key)
  gplt = NULL
  if(MLE_report$model_type == 0) {
    ## Assessment
    ## Spatial model
    gplt = ggplot(F_df, aes(x = Year, y= F, col = Fishery, linetype = Fishery)) +
      geom_line(linewidth = 1.1) +
      theme_bw() +
      facet_wrap(~Region)
  } else {
    ## Spatial model
    gplt = ggplot(F_df, aes(x = Year, y= F, col = Fishery, linetype = Fishery)) +
      geom_line(linewidth = 1.1) +
      theme_bw() +
      facet_wrap(~Region)
  }
  return(gplt)
}
#' get_fishing_mortalities
#' @param MLE_report a list that is output from obj$report() usually once an optimsation routine has been done.
#' @param region_key data.frame with colnames area and TMB_ndx for providing real region names to objects
#' @return data frame in long format
#' @export
get_fishing_mortalities = function(MLE_report, region_key = NULL) {
  regions = NULL
  fixed_F = trwl_F = NULL
  if(MLE_report$model_type == 0) {
    regions = paste0("Region ", 1:1)
    fixed_F = matrix(MLE_report$annual_F_ll, nrow = 1)
    trwl_F = matrix(MLE_report$annual_F_trwl, nrow = 1)
  } else {
    regions = paste0("Region ", 1:MLE_report$n_regions)
    fixed_F = MLE_report$annual_F_fixed
    trwl_F = MLE_report$annual_F_trwl
  }
  if(!is.null(region_key))
    regions = region_key$area[region_key$TMB_ndx + 1]

  projyears = min(MLE_report$years):(max(MLE_report$years) + MLE_report$n_projections_years)
  full_f_df = NULL;

  dimnames(fixed_F) = dimnames(trwl_F) = list(regions, projyears)
  molten_fixed_F = reshape2::melt(fixed_F)
  molten_trwl_F = reshape2::melt(trwl_F)

  colnames(molten_fixed_F) = colnames(molten_trwl_F) = c("Region","Year",  "F")
  molten_fixed_F$Fishery = "Fixed gear"
  molten_trwl_F$Fishery = "Trawl gear"
  full_f_df = rbind(molten_fixed_F, molten_trwl_F)

  return(full_f_df);
}
#'
#' get_catches
#' @param MLE_report a list that is output from obj$report() usually once an optimsation routine has been done.
#' @param region_key data.frame with colnames area and TMB_ndx for providing real region names to objects
#' @return ggplot2 object that will plot if an observation occurs in a year and region
#' @export
get_catches = function(MLE_report, region_key = NULL) {
  years = MLE_report$years
  regions = NULL

  if(MLE_report$model_type == 0) {
    regions = paste0("Region ", 1:1)
  } else {
    regions = paste0("Region ", 1:MLE_report$n_regions)
  }
  if(!is.null(region_key))
    regions = region_key$area[region_key$TMB_ndx + 1]

  if(MLE_report$model_type == 0) {
    MLE_report$fixed_fishery_catch = matrix(MLE_report$ll_fishery_catch, nrow = 1)
    MLE_report$trwl_fishery_catch = matrix(MLE_report$trwl_fishery_catch, nrow = 1)
    MLE_report$annual_fixed_catch_pred = matrix(MLE_report$annual_ll_catch_pred, nrow = 1)
    MLE_report$annual_trwl_catch_pred = matrix(MLE_report$annual_trwl_catch_pred, nrow = 1)
  }

  dimnames(MLE_report$fixed_fishery_catch) = dimnames(MLE_report$trwl_fishery_catch) = list(regions, years)
  dimnames(MLE_report$annual_fixed_catch_pred) = dimnames(MLE_report$annual_trwl_catch_pred) = list(regions, years)


  fixed_catch = reshape2::melt(MLE_report$fixed_fishery_catch)
  trwl_catch = reshape2::melt(MLE_report$trwl_fishery_catch)
  colnames(fixed_catch) = colnames(trwl_catch) = c("Region", "Year", "Catch")
  fixed_catch$Fishery = "Fixed gear"
  trwl_catch$Fishery = "Trawl"
  obs_full_df = rbind(fixed_catch, trwl_catch)
  obs_full_df$type = "Observed"

  fixed_catch = reshape2::melt(MLE_report$annual_fixed_catch_pred)
  trwl_catch = reshape2::melt(MLE_report$annual_trwl_catch_pred)
  colnames(fixed_catch) = colnames(trwl_catch) = c("Region", "Year", "Catch")
  fixed_catch$Fishery = "Fixed gear"
  trwl_catch$Fishery = "Trawl"
  pred_full_df = rbind(fixed_catch, trwl_catch)
  pred_full_df$type = "Predicted"

  full_df = rbind(obs_full_df, pred_full_df)

  return(full_df)
}
#'
#' get_movement
#' @param MLE_report obj$report()
#' @param region_key data.frame with colnames area and TMB_ndx for providing real region names to objects
#' @return data frame
#' @export
get_movement = function(MLE_report, region_key = NULL) {
  regions = paste0("Region ", 1:MLE_report$n_regions)
  if(!is.null(region_key))
    regions = region_key$area[region_key$TMB_ndx + 1]


  move_est_df = NULL
  if(MLE_report$apply_fixed_movement == 1) {
    dimnames(MLE_report$fixed_movement_matrix) = list(regions, regions)
    move_est_df = reshape2::melt(MLE_report$fixed_movement_matrix)
  } else {
    dimnames(MLE_report$movement_matrix) = list(regions, regions)
    move_est_df = reshape2::melt(MLE_report$movement_matrix)
  }
  colnames(move_est_df) = c("From","To", "Proportion")
  move_est_df$From = factor(move_est_df$From, levels = rev(regions), ordered = T)
  move_est_df$To = factor(move_est_df$To, levels = rev(regions), ordered = T)
  return(move_est_df)
}

#'
#' plot_movement
#' @param MLE_report obj$report()
#' @param region_key data.frame with colnames area and TMB_ndx for providing real region names to objects
#' @return ggplot2
#' @export
plot_movement = function(MLE_report, region_key = NULL) {
  move_est_df = get_movement(MLE_report, region_key)

  gplt = ggplot(move_est_df, aes(x = To, y = From, fill = Proportion)) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = "red") +
    geom_text(aes(x = To, y = From, label = round(Proportion,2)), color = "black", size = 4) +
    labs(x = "To", y = "From")
  return(gplt)
}
#'
#' get_selectivities
#' @param MLE_report obj$report()
#' @return data frame
#' @export
get_selectivities = function(MLE_report) {
  sel_df = NULL
  if(MLE_report$model_type == 0) {
    sel_df = data.frame(fixed_male = MLE_report$sel_ll_m, fixed_female = MLE_report$sel_ll_f,
                        trawl_male = MLE_report$sel_trwl_m, trawl_female = MLE_report$sel_trwl_f,
                        domsurveyll_male = MLE_report$sel_srv_dom_ll_m, domsurveyll_female = MLE_report$sel_srv_dom_ll_f,
                        japfishery_combined = MLE_report$sel_srv_jap_fishery_ll,
                        japsurveyll_male = MLE_report$sel_srv_jap_ll_m, japsurveyll_female = MLE_report$sel_srv_jap_ll_f,
                        nmfssurveytrwl_male = MLE_report$sel_srv_nmfs_trwl_m, nmfssurveytrwl_female = MLE_report$sel_srv_nmfs_trwl_f,
                        age = MLE_report$ages
    )

  } else {
    sel_df = data.frame(fixed_male = MLE_report$sel_fixed_m, fixed_female = MLE_report$sel_fixed_f,
                        trawl_male = MLE_report$sel_trwl_m, trawl_female = MLE_report$sel_trwl_f,
                        surveyll_male = MLE_report$sel_srv_dom_ll_m, surveyll_female = MLE_report$sel_srv_dom_ll_f,
                        age = MLE_report$ages
    )
  }
  sel_lng_df = sel_df %>% tidyr::pivot_longer(!age)
  sel_lng_df$gear = Reduce(c, lapply(sel_lng_df$name %>% stringr::str_split(pattern = "_"), function(x){x[1]}))
  sel_lng_df$sex = Reduce(c, lapply(sel_lng_df$name %>% stringr::str_split(pattern = "_"), function(x){x[2]}))
  time_block = as.numeric(Reduce(c, lapply(sel_lng_df$sex %>% stringr::str_split(pattern = "\\."), function(x){x[2]})))
  sel_lng_df$sex = Reduce(c, lapply(sel_lng_df$sex %>% stringr::str_split(pattern = "\\."), function(x){x[1]}))
  time_block[is.na(time_block)] = 1
  sel_lng_df$time_block = time_block
  return(sel_lng_df)
}

#'
#' plot_selectivities
#' @param MLE_report obj$report()
#' @return ggplot2
#' @export
plot_selectivities = function(MLE_report) {
  sel_lng_df = get_selectivities(MLE_report)

  gplt = ggplot(sel_lng_df, aes(x = age, y = value, col = sex, linetype = sex)) +
    geom_line(linewidth = 1.1) +
    facet_wrap(~gear) +
    theme_bw() +
    labs(x = "Age", y = "Selectivity", col = "Sex.timeblock", linetype = "Sex.timeblock") +
    theme(
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 14)
    )
  return(gplt)
}


#' get_recruitment
#' @param MLE_report created from obj$report()
#' @param region_key data.frame with colnames area and TMB_ndx for providing real region names to objects
#' @return ggplot2
#' @export
#'
get_recruitment = function(MLE_report, region_key = NULL) {
  if(MLE_report$model_type == 0) {
    regions = paste0("Region ", 1:1)
  } else {
    regions = paste0("Region ", 1:MLE_report$n_regions)
  }
  if(!is.null(region_key))
    regions = region_key$area[region_key$TMB_ndx + 1]
  projyears = min(MLE_report$years):(max(MLE_report$years) + MLE_report$n_projections_years)

  ssbs = NULL
  if(MLE_report$model_type == 0) {
    MLE_report$annual_recruitment = matrix(MLE_report$annual_recruitment, ncol = 1)
    MLE_report$ln_rec_dev = matrix(MLE_report$ln_rec_dev, ncol= 1)
    dimnames(MLE_report$annual_recruitment) = list(projyears, regions)
    dimnames(MLE_report$ln_rec_dev) = list(projyears, regions)
    MLE_report$recruitment_multipliers = MLE_report$ln_rec_dev
    MLE_report$recruitment_multipliers = exp(MLE_report$recruitment_multipliers + 0.5*MLE_report$sigma_R^2)
    names(MLE_report$mean_rec) = regions
    recruit_df = reshape2::melt(MLE_report$annual_recruitment)
    recruit_dev_df = reshape2::melt(MLE_report$ln_rec_dev)
    recruit_ycs_df = reshape2::melt(MLE_report$recruitment_multipliers)
    colnames(recruit_dev_df) =  c("Year", "Region", "Recruitment_deviation")
    colnames(recruit_ycs_df) = c("Year", "Region", "YCS")
  } else {
    dimnames(MLE_report$recruitment_yr) = list(projyears, regions)
    dimnames(MLE_report$recruitment_devs) = list(regions, projyears)
    dimnames(MLE_report$recruitment_multipliers) = list(regions, projyears)

    names(MLE_report$mean_rec) = regions
    recruit_df = reshape2::melt(MLE_report$recruitment_yr)
    recruit_dev_df = reshape2::melt(MLE_report$recruitment_devs)
    recruit_ycs_df = reshape2::melt(MLE_report$recruitment_multipliers)
    colnames(recruit_dev_df) = c("Region", "Year", "Recruitment_deviation")
    colnames(recruit_ycs_df) = c("Region", "Year", "YCS")
  }
  mean_recruit_df = reshape2::melt(as.matrix(MLE_report$mean_rec))
  colnames(recruit_df) = c("Year","Region", "Recruitment")


  recruit_df$Recruitment_devs = recruit_dev_df
  colnames(mean_recruit_df) = c("Region", "row_ndx","Mean_recruitment")
  recruit_df = recruit_df %>% dplyr::inner_join(mean_recruit_df, by = "Region")
  recruit_df = recruit_df %>% dplyr::inner_join(recruit_dev_df, by = c("Region", "Year"))
  recruit_df = recruit_df %>% dplyr::inner_join(recruit_ycs_df, by = c("Region", "Year"))

  return(recruit_df)
}
#' plot_recruitment
#' @param MLE_report created from obj$report()
#' @param region_key data.frame with colnames area and TMB_ndx for providing real region names to objects
#' @return ggplot2
#' @export
#'
plot_recruitment = function(MLE_report, region_key = NULL) {
  recruit_df = get_recruitment(MLE_report, region_key)
  gplt = ggplot(recruit_df) +
    geom_line(aes(x = Year, y = Recruitment), linewidth = 1.1) +
    geom_hline(aes(yintercept = Mean_recruitment), linetype = "dashed", linewidth = 1.1, col = "gray60") +
    facet_wrap(~Region) +
    theme_bw()

  return(gplt)
}

#'
#' get_SSB
#' @param MLE_report a list that is output from obj$report() usually once an optimsation routine has been done.
#' @param region_key data.frame with colnames area and TMB_ndx for providing real region names to objects
#' @param depletion boolean if true we will scale values by Bzero
#' @return ggplot2 object that will plot if an observation occurs in a year and region
#' @export
get_SSB = function(MLE_report, region_key = NULL, depletion = F) {
  years = MLE_report$years
  projyears = min(MLE_report$years):(max(MLE_report$years) + MLE_report$n_projections_years)

  ssbs = NULL
  if(MLE_report$model_type == 0) {
    regions = 1:1
    ssbs = matrix(MLE_report$SSB, ncol = 1)
  } else {
    regions = 1:MLE_report$n_regions
    ssbs = MLE_report$SSB_yr
  }
  if(depletion)
    ssbs = sweep(ssbs, MARGIN = 2, STATS = MLE_report$Bzero, FUN = "/") * 100


  dimnames(ssbs) = list(projyears, regions)
  molten_ssbs = reshape2::melt(ssbs)

  colnames(molten_ssbs) = c("Year", "Region", "SSB")
  if(depletion)
    colnames(molten_ssbs) = c("Year", "Region", "Depletion")

  if(is.null(region_key)) {
    molten_ssbs$Region = paste0("Region ", molten_ssbs$Region)
  } else {
    molten_ssbs$Region = region_key$area[match(molten_ssbs$Region, (region_key$TMB_ndx + 1))]
  }
  return(molten_ssbs);
}

#'
#' plot_SSB
#' @param MLE_report a list that is output from obj$report() usually once an optimsation routine has been done.
#' @param data list that is passed to the MakeADfun for the TMB model
#' @param region_key data.frame with colnames area and TMB_ndx for providing real region names to objects
#' @param depletion boolean if true we will scale values by Bzero
#' @return ggplot2 object that will plot if an observation occurs in a year and region
#' @export
plot_SSB = function(MLE_report, region_key = NULL, depletion = F) {
  molten_ssbs = get_SSB(MLE_report, region_key, depletion = depletion)
  gplt = NULL
  if(depletion) {
    gplt = ggplot(molten_ssbs) +
      geom_line(aes(x = Year, y = Depletion, col = Region), linewidth= 1.1) +
      guides(colour = "none", linewidth = "none") +
      labs(y = "Depletion (SSB/B0 %)") +
      facet_wrap(~Region) +
      theme_bw() +
      ylim(0, NA)
  } else {
    gplt = ggplot(molten_ssbs) +
      geom_line(aes(x = Year, y = SSB, col = Region), linewidth= 1.1) +
      guides(colour = "none", linewidth = "none") +
      labs(y = "Spawning Stock Biomass (SSB)") +
      facet_wrap(~Region) +
      theme_bw() +
      ylim(0, NA)
  }

  return(gplt)
}


#'
#' get_other_derived_quantities
#' @param MLE_report a list that is output from obj$report() usually once an optimsation routine has been done.
#' @param data list that is passed to the MakeADfun for the TMB model
#' @param region_key data.frame with colnames area and TMB_ndx for providing real region names to objects
#' @return A list of
#' @export

get_other_derived_quantities <- function(MLE_report, data, region_key = NULL) {
  if(MLE_report$model_type == 0) {
    cat("Skipping this function not built for the Current Assessment")
    return(NULL)
  }
  Region_lab = paste0("Region ", 1:MLE_report$n_regions)
  if(!is.null(region_key))
    Region_lab =  region_key$area[region_key$TMB_ndx + 1]

  catchability = MLE_report$srv_dom_ll_q
  dimnames(catchability) = list(Region_lab, paste0("Block-", 1:ncol(catchability)))
  molten_catchabilties = reshape2::melt(catchability)
  colnames(molten_catchabilties) = c("Region", "time-block", "q")

  ## Scalar model quantities
  scalar_quants = data.frame(F_init = round(MLE_report$init_F_hist, 3),
                             tag_phi = MLE_report$tag_phi,
                             theta_fixed_catchatlgth = MLE_report$theta_fixed_catchatlgth,
                             theta_fixed_catchatage = MLE_report$theta_fixed_catchatage,
                             theta_trwl_catchatlgth = MLE_report$theta_trwl_catchatlgth,
                             theta_srv_catchatage = MLE_report$theta_srv_dom_ll_catchatage,
                             sigma_R = round(MLE_report$sigma_R,3),
                             sigma_init_age_devs = round(MLE_report$sigma_init_devs,3),
                             catch_sd = round(MLE_report$catch_sd, 3),
                             apply_fixed_movement = MLE_report$apply_fixed_movement,
                             do_recruits_move = data$do_recruits_move,
                             evaluate_tag_likelihood = data$evaluate_tag_likelihood,
                             SR = ifelse(data$SrType == 2, "BH", "NO SR")
                             )

  ## spatial_scalars
  spatial_params = data.frame(Region = Region_lab, Bzero = MLE_report$Bzero, Binit = MLE_report$Binit, Bzero_recent_growth = MLE_report$Bzero_w_recent_growth, R0 = MLE_report$mean_rec)

  ## tag stuff
  tag_reporting_rate = unique(as.numeric(MLE_report$tag_reporting_rate))

  return(list(catchabilities = molten_catchabilties, scalar_quants = scalar_quants, spatial_params = spatial_params, tag_reporting_rate = tag_reporting_rate))

}



#'
#' get_qs
#' @param MLE_report a list that is output from obj$report() usually once an optimsation routine has been done.
#' @return A data frame of catchabilities
#' @export

get_qs <- function(MLE_report) {
  molten_catchabilties = NULL
  if(MLE_report$model_type == 0) {
    ## four catchabilities in this model
    ## they are dynamic in terms of number of time-blocks
    ll_cpue = data.frame(Q = MLE_report$ll_cpue_q, label = "Fixed gear CPUE", time_block = paste0("TimeBlock: ", 1:length(MLE_report$ll_cpue_q)))
    srv_dom_ll = data.frame(Q = MLE_report$srv_dom_ll_q, label = "Domestic LL survey", time_block = paste0("TimeBlock: ", 1:length(MLE_report$srv_dom_ll_q)))
    jap_fishery_cpue = data.frame(Q = MLE_report$srv_jap_fishery_ll_q, label = "Japanese CPUE", time_block = paste0("TimeBlock: ", 1:length(MLE_report$srv_jap_fishery_ll_q)))
    srv_jap_ll = data.frame(Q = MLE_report$srv_jap_ll_q, label = "Japanese LL survey", time_block = paste0("TimeBlock: ", 1:length(MLE_report$srv_jap_ll_q)))
    # combine thiese
    molten_catchabilties = rbind(ll_cpue, srv_dom_ll, jap_fishery_cpue, srv_jap_ll)
  } else {
    Region_lab = paste0("Region ", 1:MLE_report$n_regions)
    if(!is.null(region_key))
      Region_lab =  region_key$area[region_key$TMB_ndx + 1]

    catchability = MLE_report$srv_dom_ll_q
    dimnames(catchability) = list(Region_lab, paste0("Block-", 1:ncol(catchability)))
    molten_catchabilties = reshape2::melt(catchability)
    colnames(molten_catchabilties) = c("Region", "time-block", "q")
  }
  return(molten_catchabilties)
}
