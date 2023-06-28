#'
#' Test MSE
#'
#'
library(ggplot2)
library(dplyr)
load(system.file("testdata", "MockAssessmentModel.RData",package="SpatialSablefishAssessment"))

names(data)
## change likelihood formulation
if(F){
  data$ll_catchatage_comp_likelihood = 1
  data$ll_catchatlgth_comp_likelihood = 1
  data$srv_dom_ll_lgth_comp_likelihood = 1
  data$srv_dom_ll_age_comp_likelihood = 1
  ## increase sample size
  data$obs_ll_catchatage = data$obs_ll_catchatage * 5
  data$obs_ll_catchatlgth_m = data$obs_ll_catchatlgth_m * 5
  data$obs_ll_catchatlgth_f = data$obs_ll_catchatlgth_f * 5
  data$obs_srv_dom_ll_age = data$obs_srv_dom_ll_age * 5
  data$obs_srv_dom_ll_lgth_m = data$obs_srv_dom_ll_lgth_m * 5
  data$obs_srv_dom_ll_lgth_f = data$obs_srv_dom_ll_lgth_f * 5
  data$obs_srv_jap_ll_age = data$obs_srv_jap_ll_age * 5
  data$obs_srv_nmfs_trwl_age = data$obs_srv_nmfs_trwl_age * 5
}
# simplify selectivity blocks for fixed gear fishery
data$ll_sel_type = c(0)
data$ll_sel_by_year_indicator = rep(0, length(data$ll_sel_by_year_indicator ))
parameters$ln_ll_sel_pars = array(parameters$ln_ll_sel_pars[1,,], dim = c(1,2,2))
data$srv_dom_ll_sel_type = c(0)
data$srv_dom_ll_sel_by_year_indicator = rep(0, length(data$srv_dom_ll_sel_by_year_indicator ))
data$srv_d
parameters$ln_srv_dom_ll_sel_pars = array(parameters$ln_srv_dom_ll_sel_pars[1,,], dim = c(1,2,2))

## Fix delta params.
## common delta params for all sexes and time-blocks per selectivity
fix_these_parameters = unique(names(parameters)[!names(parameters) %in% c("ln_ll_F_avg","ln_ll_F_devs","ln_trwl_F_avg","ln_trwl_F_devs")])#c("ln_init_rec_dev")#c("ln_srv_jap_fishery_ll_sel_pars")
## turn off the delta parameters that we want to estimate
## as based on one value
array_elements_to_exclude  = list()
if(F) {
  # three time-blocks for this selectivity
  array_elements_to_exclude[["ln_ll_sel_pars"]] =
    matrix( c(1,2,2,
              2,2,1,
              2,2,2,
              3,2,1,
              3,2,2), byrow = T, ncol = 3)
} else {
  # One time-blocks for this selectivity
  array_elements_to_exclude[["ln_ll_sel_pars"]] =
    matrix( c(1,2,2), byrow = T, ncol = 3)
}

array_elements_to_exclude[["ln_trwl_sel_pars"]] =
  matrix( c(1,2,2), byrow = T, ncol = 3)
array_elements_to_exclude[["ln_srv_jap_ll_sel_pars"]] =
  matrix( c(1,2,2), byrow = T, ncol = 3)
array_elements_to_exclude[["ln_srv_dom_ll_sel_pars"]] =
  matrix( c(1,2,2), byrow = T, ncol = 3)
## turn off params
na_map = fix_pars(parameters, pars_to_exclude = c(fix_these_parameters, names(array_elements_to_exclude)),array_elements_to_exclude = array_elements_to_exclude)

## set the same and copy pars
same_pars = base_pars = list()
if(F) {
  # three time-blocks for this selectivity
  base_pars[["ln_ll_sel_pars"]] = matrix(rep(c(1,2,1), 5), byrow = T, ncol = 3)
  same_pars[["ln_ll_sel_pars"]] =
    matrix( c(1,2,2,
            2,2,1,
            2,2,2,
            3,2,1,
            3,2,2), byrow = T, ncol = 3)
} else {
  # One time-blocks for this selectivity
  base_pars[["ln_ll_sel_pars"]] = matrix(rep(c(1,2,1), 1), byrow = T, ncol = 3)
  same_pars[["ln_ll_sel_pars"]] =
    matrix( c(1,2,2), byrow = T, ncol = 3)
}

base_pars[["ln_trwl_sel_pars"]] = matrix(c(1,2,1), byrow = T, ncol = 3)
same_pars[["ln_trwl_sel_pars"]] =
  matrix( c(1,2,2), byrow = T, ncol = 3)
base_pars[["ln_srv_jap_ll_sel_pars"]] = matrix(c(1,2,1), byrow = T, ncol = 3)
same_pars[["ln_srv_jap_ll_sel_pars"]] =
  matrix( c(1,2,2), byrow = T, ncol = 3)
base_pars[["ln_srv_dom_ll_sel_pars"]] = matrix(rep(c(1,2,1), 1), byrow = T, ncol = 3)
same_pars[["ln_srv_dom_ll_sel_pars"]] =
  matrix( c(1,2,2), byrow = T, ncol = 3)
# convert array pars to vectors for set_pars_to_be_the_same
new_base = new_copy = list()
for(i in 1:length(base_pars)) {
  for(j in 1:nrow(base_pars[[i]])) {
    this_bse_lst  = evalit(x = paste0("this_bse_lst = list(",names(base_pars)[i]," = get_TMB_vector_from_array(element = base_pars[[i]][j, ], array = parameters[[names(base_pars)[i]]]))"))
    this_cpy_lst  = evalit(x = paste0("this_cpy_lst = list(",names(same_pars)[i]," = get_TMB_vector_from_array(element = same_pars[[i]][j, ], array = parameters[[names(same_pars)[i]]]))"))
    new_base = append(new_base, this_bse_lst)
    new_copy = append(new_copy, this_cpy_lst)
  }
}
na_map = set_pars_to_be_the_same(par_list = parameters, map = na_map, base_parameters = new_base, copy_parameters = new_copy)

fix_these_parameters = unique(names(parameters)[!names(parameters) %in% c("ln_ll_F_avg","ln_ll_F_devs","ln_trwl_F_avg","ln_trwl_F_devs")])#c("ln_init_rec_dev")#c("ln_srv_jap_fishery_ll_sel_pars")
na_map = fix_pars(parameters, pars_to_exclude = c(fix_these_parameters))


## turn off Japanese LF selectivity params
#data$srv_jap_fishery_ll_lgth_indicator = rep(0, length(data$srv_jap_fishery_ll_lgth_indicator))
parameters$ln_init_rec_dev = rep(-0.5*data$sigma_R^2, data$n_init_rec_devs)
# change catch likelihood
data$catch_likelihood = 1
OM <- TMB::MakeADFun(data = data,
                     parameters = parameters,
                     map = na_map,
                     DLL = "SpatialSablefishAssessment_TMBExports", silent  = T)

OM_report = OM$report()
sim_data = OM$simulate(complete = T)
est_lst = list()
est_lst[["OM"]] = sim_data
##

## self test with SSBs
SSBs = depletion = NULL
i = 1
for(i in 1:50) {
  ## simulate OM_data
  sim_data = OM$simulate(complete = T)
  sim_data = convert_simdata_integers(sim_data, data)

  ## estimate
  EM <- TMB::MakeADFun(data = sim_data,
                       parameters = parameters,
                       map = na_map,
                       DLL = "SpatialSablefishAssessment_TMBExports", silent  =T,checkParameterOrder=TRUE)
  #check_gradients(EM)
  first_MLE = nlminb(start = EM$par, objective = EM$fn, gradient  = EM$gr, control = list(eval.max = 1000, iter.max = 1000))
  cat(" par ", names(EM$par)[which.max(EM$gr())], " had largest gradient = ", EM$gr()[which.max(EM$gr())], "\n")
  try_improve = tryCatch(expr =
                           for(i in 1:2) {
                             g = as.numeric(EM$gr(first_MLE$par))
                             h = optimHess(first_MLE$par, fn = EM$fn, gr = EM$gr)
                             first_MLE$par = first_MLE$par - solve(h,g)
                             first_MLE$objective = EM$fn(first_MLE$par)
                           }
                         , error = function(e){e})

  if(inherits(try_improve, "error") | inherits(try_improve, "warning")) {
    cat("didn't converge on simulation ", i, "\n")
    ## investigate furthur
    h = optimHess(first_MLE$par, fn = EM$fn, gr = EM$gr)
    cat("Problem parameter = ", names(first_MLE$par)[which(eigen(h)$values <= 1e-6)], "\n")

  } else {
    # get derived quantities
    mle_report = EM$report(first_MLE$par)
    est_lst[[as.character(i)]] = mle_report
    #plot_catch_fit(mle_report)
  }
}

## retrieve parameters
ssbs = get_multiple_ssbs(est_lst, run_labels = names(est_lst), depletion = F)
ggplot(ssbs, aes(x = Year, y = SSB, col = label, linetype = label)) +
  geom_line(linewidth = 1.1) +
  theme_bw()

catch = get_multiple_catch_fits(est_lst, run_labels = names(est_lst))
ggplot() +
  geom_point(data = catch %>% filter(type == "Observed"), aes(x = Year, y = Catch, shape = label), size = 1.6) +
  geom_line(data = catch %>% filter(type == "Predicted"), aes(x = Year, y = Catch, col = label, linetype = label)) +
  theme_bw() +
  facet_wrap(~Fishery)


