#
#
#
#' plot_fishing_mortalities
#' @param MLE_report a list that is output from obj$report() usually once an optimsation routine has been done.
#' @param region_key data.frame with colnames area and TMB_ndx for providing real region names to objects
#' @return ggplot2 object that will plot if an observation occurs in a year and region
#' @export

plot_fishing_mortalities = function(MLE_report, region_key = NULL) {
  F_df = get_fishing_mortalities(MLE_report = mle_report, region_key = region_key)
  gplt = ggplot(F_df, aes(x = Year, y= F, col = Fishery, linetype = Fishery)) +
    geom_line(linewidth = 1.1) +
    theme_bw() +
    facet_wrap(~Region)
  return(gplt)
}
#' get_fishing_mortalities
#' @param MLE_report a list that is output from obj$report() usually once an optimsation routine has been done.
#' @param region_key data.frame with colnames area and TMB_ndx for providing real region names to objects
#' @return data frame in long format
#' @export
get_fishing_mortalities = function(MLE_report, region_key = NULL) {
  years = MLE_report$years
  regions = 1:MLE_report$n_regions
  fixed_F = MLE_report$annual_F_fixed
  trwl_F = MLE_report$annual_F_trwl

  dimnames(fixed_F) = dimnames(trwl_F) = list(regions, years)
  molten_fixed_F = reshape2::melt(fixed_F)
  molten_trwl_F = reshape2::melt(trwl_F)

  colnames(molten_fixed_F) = colnames(molten_trwl_F) = c("Region","Year",  "F")
  molten_fixed_F$Fishery = "Fixed gear"
  molten_trwl_F$Fishery = "Trawl"
  full_f_df = rbind(molten_fixed_F, molten_trwl_F)
  if(is.null(region_key)) {
    full_f_df$Region = paste0("Region ", full_f_df$Region)
  } else {
    full_f_df$Region = region_key$area[match(full_f_df$Region, (region_key$TMB_ndx + 1))]
  }
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
  regions = 1:MLE_report$n_regions
  dimnames(MLE_report$fixed_fishery_catch) = dimnames(MLE_report$trwl_fishery_catch) = list(regions, years)
  fixed_catch = reshape2::melt(MLE_report$fixed_fishery_catch)
  trwl_catch = reshape2::melt(MLE_report$trwl_fishery_catch)
  colnames(fixed_catch) = colnames(trwl_catch) = c("Region", "Year", "Catch")
  fixed_catch$label = "Fixed gear"
  trwl_catch$label = "Trawl"
  obs_full_df = rbind(fixed_catch, trwl_catch)
  if(is.null(region_key)) {
    obs_full_df$Region = paste0("Region ", obs_full_df$Region)
  } else {
    obs_full_df$Region = region_key$area[match(obs_full_df$Region, (region_key$TMB_ndx + 1))]
  }
  obs_full_df$type = "Observed"

  dimnames(MLE_report$annual_fixed_catch_pred) = dimnames(MLE_report$annual_trwl_catch_pred) = list(regions, years)
  fixed_catch = reshape2::melt(MLE_report$annual_fixed_catch_pred)
  trwl_catch = reshape2::melt(MLE_report$annual_trwl_catch_pred)
  colnames(fixed_catch) = colnames(trwl_catch) = c("Region", "Year", "Catch")
  fixed_catch$label = "Fixed gear"
  trwl_catch$label = "Trawl"
  pred_full_df = rbind(fixed_catch, trwl_catch)
  if(is.null(region_key)) {
    pred_full_df$Region = paste0("Region ", pred_full_df$Region)
  } else {
    pred_full_df$Region = region_key$area[match(pred_full_df$Region, (region_key$TMB_ndx + 1))]
  }
  pred_full_df$type = "Predicted"

  full_df = rbind(obs_full_df, pred_full_df)

  return(full_df)
}
#'
#' plot_movement
#' @param MLE_report obj$report()
#' @param region_key data.frame with colnames area and TMB_ndx for providing real region names to objects
#' @return ggplot2
#' @export
plot_movement = function(MLE_report, region_key = NULL) {
  regions = paste0("Region ", 1:data$n_regions)
  if(!is.null(region_key))
    regions = region_key$area[region_key$TMB_ndx + 1]

  dimnames(MLE_report$movement_matrix) = list(regions, regions)


  move_est_df = reshape2::melt(MLE_report$movement_matrix)
  colnames(move_est_df) = c("From","To", "Proportion")
  move_est_df$From = factor(move_est_df$From, levels = rev(regions), ordered = T)
  move_est_df$To = factor(move_est_df$To, levels = rev(regions), ordered = T)

  gplt = ggplot(move_est_df, aes(x = To, y = From, fill = Proportion)) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = "red") +
    geom_text(aes(x = To, y = From, label = round(Proportion,2)), color = "black", size = 4) +
    labs(x = "To", y = "From")
  return(gplt)
}

#'
#' plot_selectivities
#' @param MLE_report obj$report()
#' @return ggplot2
#' @export
plot_selectivities = function(MLE_report) {
  sel_df = data.frame(fixed_male = MLE_report$sel_fixed_m, fixed_female = MLE_report$sel_fixed_f,
                      trawl_male = MLE_report$sel_trwl_m, trawl_female = MLE_report$sel_trwl_f,
                      surveyll_male = MLE_report$sel_srv_dom_ll_m, surveyll_female = MLE_report$sel_srv_dom_ll_f,
                      age = MLE_report$ages
  )
  sel_lng_df = sel_df %>% pivot_longer(!age)
  sel_lng_df$gear = Reduce(c, lapply(sel_lng_df$name %>% str_split(pattern = "_"), function(x){x[1]}))
  sel_lng_df$sex = Reduce(c, lapply(sel_lng_df$name %>% str_split(pattern = "_"), function(x){x[2]}))
  gplt = ggplot(sel_lng_df, aes(x = age, y = value, col = sex, linetype = sex)) +
    geom_line(linewidth = 1.1) +
    facet_wrap(~gear) +
    theme_bw()
  return(gplt)
}


#' plot_recruitment
#' @param MLE_report created from obj$report()
#' @param region_key data.frame with colnames area and TMB_ndx for providing real region names to objects
#' @return ggplot2
#' @export
#'
plot_recruitment = function(MLE_report, region_key = NULL) {
  regions = paste0("Region ", 1:data$n_regions)
  if(!is.null(region_key))
    regions = region_key$area[region_key$TMB_ndx + 1]

  dimnames(MLE_report$recruitment_yr) = list(data$years, regions)

  recruit_df = reshape2::melt(MLE_report$recruitment_yr)
  colnames(recruit_df) = c("Year","Region", "Recruitment")
  gplt = ggplot(recruit_df, aes(x = Year, y = Recruitment)) +
    geom_line(linewidth = 1.1) +
    facet_wrap(~Region) +
    theme_bw()
  return(gplt)
}
