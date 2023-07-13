#'
#' Conduct a self test for the current assessment
#' A useful script for Dan to help investigate
#' his assessment model in more detail
#'
#'
library(ggplot2)
library(SpatialSablefishAssessment)
library(tidyverse)
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
#fix_these_parameters = unique(names(parameters)[!names(parameters) %in% c("ln_ll_F_avg","ln_ll_F_devs","ln_trwl_F_avg","ln_trwl_F_devs")])#c("ln_init_rec_dev")#c("ln_srv_jap_fishery_ll_sel_pars")
fix_these_parameters = c("ln_init_rec_dev")

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

#fix_these_parameters = unique(names(parameters)[!names(parameters) %in% c("ln_mean_rec", "ln_rec_dev", "ln_ll_F_avg","ln_ll_F_devs","ln_trwl_F_avg","ln_trwl_F_devs")])#c("ln_init_rec_dev")#c("ln_srv_jap_fishery_ll_sel_pars")
#na_map = fix_pars(parameters, pars_to_exclude = c(fix_these_parameters))


## turn off Japanese LF selectivity params
#data$srv_jap_fishery_ll_lgth_indicator = rep(0, length(data$srv_jap_fishery_ll_lgth_indicator))
parameters$ln_init_rec_dev = rep(-0.5*data$sigma_R^2, data$n_init_rec_devs)
# change catch likelihood
data$catch_likelihood = 0
OM <- TMB::MakeADFun(data = data,
                     parameters = parameters,
                     map = na_map,
                     DLL = "SpatialSablefishAssessment_TMBExports", silent  = T)
## run simulation
est_lst = OM_lst = list()
OM_lst[["OM"]] = OM$report()
set.seed(123)
n_sims = 20
for(j in 1:n_sims) {
  if(j %% 10 ==  0)
    cat("sim iteration = ", j, "\n")
  ## simulate OM_data
  sim_data = OM$simulate(complete = T)
  sim_data = convert_simdata_integers(sim_data, data)
  ## estimate
  EM <- TMB::MakeADFun(data = sim_data,
                       parameters = parameters,
                       map = na_map,
                       DLL = "SpatialSablefishAssessment_TMBExports", silent  =T,checkParameterOrder=TRUE)

  ## use OM params for comparing EMs
  OM_for_report <- TMB::MakeADFun(data = sim_data,
                                  parameters = parameters,
                                  map = na_map,
                                  DLL = "SpatialSablefishAssessment_TMBExports", silent  = T)
  OM_lst[[as.character(j)]] = OM_for_report$report()
  #check_gradients(EM)
  first_MLE = nlminb(start = EM$par, objective = EM$fn, gradient  = EM$gr, control = list(eval.max = 1000, iter.max = 1000))
  #cat(" par ", names(EM$par)[which.max(EM$gr())], " had largest gradient = ", EM$gr()[which.max(EM$gr())], "\n")
  try_improve = tryCatch(expr =
                           for(i in 1:2) {
                             g = as.numeric(EM$gr(first_MLE$par))
                             h = optimHess(first_MLE$par, fn = EM$fn, gr = EM$gr)
                             first_MLE$par = first_MLE$par - solve(h,g)
                             first_MLE$objective = EM$fn(first_MLE$par)
                           }
                         , error = function(e){e}, warning = function(w){w})

  if(inherits(try_improve, "error") | inherits(try_improve, "warning")) {
    cat("didn't converge on simulation ", j, "\n")
    ## investigate furthur
    h = optimHess(first_MLE$par, fn = EM$fn, gr = EM$gr)
  } else {
    # get derived quantities
    mle_report = EM$report(first_MLE$par)
    est_lst[[as.character(j)]] = mle_report
    #plot_catch_fit(mle_report)
  }
}

## Look at relative error in derived quantitites
ssbs = get_multiple_ssbs(est_lst, run_labels = names(est_lst), depletion = F)
ssbs_OM = get_multiple_ssbs(OM_lst, run_labels = names(OM_lst), depletion = F)
## merge OM into EM to calcualte relative error
ssbs = ssbs %>% left_join(ssbs_OM, by = c("Region", "Year")) %>% mutate(RE = (SSB.x - SSB.y)/SSB.y * 100)
ggplot(ssbs, aes(x = Year, y = RE, group = Year)) +
  geom_boxplot() +
  geom_hline(yintercept = 0, col = "red", linewidth = 1.1, linetype = "dashed") +
  theme_bw()

Fs = get_multiple_Fs(est_lst, run_labels = names(est_lst))
Fs_OM = get_multiple_Fs(OM_lst, run_labels = names(OM_lst))
## merge OM into EM to calcualte relative error
Fs = Fs %>% left_join(Fs_OM, by = c("Region", "Year","Fishery")) %>% mutate(RE = (F.x - F.y)/F.y * 100)
ggplot(Fs, aes(x = Year, y = RE, group = Year)) +
  geom_boxplot() +
  geom_hline(yintercept = 0, col = "red", linewidth = 1.1, linetype = "dashed") +
  theme_bw() +
  facet_wrap(~Fishery)

catch = get_multiple_catch_fits(est_lst, run_labels = names(est_lst))
catch_OM = get_multiple_catch_fits(OM_lst, run_labels = names(OM_lst))
catch = catch %>% left_join(catch_OM, by = c("Region", "Year","Fishery","type")) %>% mutate(RE = (Catch.x - Catch.y)/Catch.y * 100)

ggplot(catch %>% filter(type == "Predicted"), aes(x = Year, y = RE, group = Year)) +
  geom_boxplot() +
  geom_hline(yintercept = 0, col = "red", linewidth = 1.1, linetype = "dashed") +
  theme_bw() +
  facet_wrap(~Fishery)

ggplot() +
  geom_point(data = catch %>% filter(type == "Observed"), aes(x = Year, y = Catch.x), size = 1.6, alpha = 0.3) +
  geom_line(data = catch %>% filter(type == "Predicted"), aes(x = Year, y = Catch.x), linewidth = 0.9, alpha = 0.3) +
  theme_bw() +
  facet_wrap(~Fishery)

## relative index fits
index_df = get_multiple_index_fits(est_lst, run_labels = names(est_lst))
ggplot(data = index_df, aes(x = Year, y = Pearsons_residuals)) +
  geom_point(size= 1.4) +
  facet_wrap(~observation) +
  theme_bw()

ggplot(data = index_df, aes(x = Year, y = Predicted, col = label, linetype = label)) +
  geom_line(linewidth = 1.1) +
  facet_wrap(~observation) +
  theme_bw()


## compare negative loglikelihoods
nlls = get_multiple_nlls(est_lst, run_labels = names(est_lst))
nll_wider = nlls %>% pivot_wider(id_cols = observations, names_from = label, values_from = negloglike) %>% mutate(diff = abs(OM - `1`))
print(nll_wider, n = 24)


#############################################
## Some more model comparisons
#############################################
## save in a named list so we can use the get_multiple accessors
mod_lst = list()
mod_lst[["OM"]] = OM_report
mod_lst[["EM"]] = mle_report
## do a bunch of model comparisons
# get likelihoods
nlls = get_multiple_nlls(mle_ls = mod_lst, run_labels = names(mod_lst))
nll_wider = nlls %>% pivot_wider(id_cols = observations, names_from = label, values_from = negloglike) %>% mutate(diff = abs(OM - EM))
print(nll_wider, n = 24)

## plot initial age
init_age = get_multiple_init_nage(mle_ls = mod_lst, run_labels = names(mod_lst))
ggplot(init_age, aes(x = Age, y = Numbers, col = label, linetype = label)) +
  geom_line(linewidth = 1.1) +
  facet_wrap(~sex) +
  theme_bw()

## plot SSBs
ssbs = get_multiple_ssbs(mod_lst, run_labels = names(mod_lst), depletion = F)
ggplot(ssbs, aes(x = Year, y = SSB, col = label, linetype = label)) +
  geom_line(linewidth = 1.1) +
  theme_bw()
## plot catch
catch = get_multiple_catch_fits(mod_lst, run_labels = names(mod_lst))
ggplot() +
  geom_point(data = catch %>% filter(type == "Observed"), aes(x = Year, y = Catch, shape = label), size = 1.6) +
  geom_line(data = catch %>% filter(type == "Predicted"), aes(x = Year, y = Catch, col = label, linetype = label), linewidth = 0.9) +
  theme_bw() +
  facet_wrap(~Fishery)

## plot recruitment
recruits = get_multiple_recruits(mod_lst, run_labels = names(mod_lst))
ggplot(recruits, aes(x = Year, y = Recruitment, col = label, linetype = label)) +
  geom_line(linewidth = 1.1) +
  theme_bw()
ggplot(recruits, aes(x = Year, y = Recruitment_deviation, col = label, linetype = label)) +
  geom_line(linewidth = 1.1) +
  theme_bw()
## Selectivities
select_df = get_multiple_selectivities(mle_ls = mod_lst, run_labels = names(mod_lst))
ggplot(select_df, aes(x = age, y = value, col = label, linetype = label)) +
  geom_line(linewidth = 1.1) +
  facet_grid(gear~sex) +
  theme_bw()
## plot index_fit
index_df = get_multiple_index_fits(mod_lst, run_labels = names(mod_lst))
ggplot(data = index_df, aes(x = Year, y = Pearsons_residuals, col = label, shape = label)) +
  geom_point(size= 1.4) +
  facet_wrap(~observation) +
  theme_bw()
ggplot(data = index_df, aes(x = Year, y = Predicted, col = label, linetype = label)) +
  geom_line(linewidth = 1.1) +
  facet_wrap(~observation) +
  theme_bw()
ggplot(data = index_df, aes(x = Year, y = Observed, col = label, linetype = label)) +
  geom_line(linewidth = 1.1) +
  facet_wrap(~observation) +
  theme_bw()
## plot AF fits
AF_df = get_multiple_AFs(mle_ls = mod_lst, run_labels = names(mod_lst))
ggplot(AF_df %>% dplyr::filter(label == "EM", observation == "Fixed gear fishery"), aes(x = Age)) +
  geom_point(aes(y = Observed, col = "Observed")) +
  geom_line(aes(y = Predicted, col = "Predicted"), linewidth = 0.9, linetype = "dashed") +
  facet_wrap(~Year) +
  ggtitle("Fixed gear fishery") +
  ylim(0, 20) +
  theme_bw()
ggplot(AF_df %>% dplyr::filter(label == "EM", observation == "Domestic LL survey"), aes(x = Age)) +
  geom_point(aes(y = Observed, col = "Observed")) +
  geom_line(aes(y = Predicted, col = "Predicted"), linewidth = 0.9, linetype = "dashed") +
  facet_wrap(~Year) +
  ggtitle("Domestic LL survey") +
  theme_bw()
## plot LF fits
LF_df = get_multiple_LFs(mle_ls = mod_lst, run_labels = names(mod_lst))

ggplot(LF_df %>% dplyr::filter(label == "EM", observation == "Fixed gear fishery", Sex == "male"), aes(x = Length)) +
  geom_point(aes(y = Observed, col = "Observed")) +
  geom_line(aes(y = Predicted, col = "Predicted"), linewidth = 0.9, linetype = "dashed") +
  facet_wrap(~Year) +
  ggtitle("Fixed gear fishery Male") +
  theme_bw()

ggplot(LF_df %>% dplyr::filter(label == "EM", observation == "Fixed gear fishery", Sex == "female"), aes(x = Length)) +
  geom_point(aes(y = Observed, col = "Observed")) +
  geom_line(aes(y = Predicted, col = "Predicted"), linewidth = 0.9, linetype = "dashed") +
  facet_wrap(~Year) +
  ggtitle("Fixed gear fishery Female") +
  theme_bw()

ggplot(LF_df %>% dplyr::filter(label == "EM", observation == "Trawl gear fishery", Sex == "male"), aes(x = Length)) +
  geom_point(aes(y = Observed, col = "Observed")) +
  geom_line(aes(y = Predicted, col = "Predicted"), linewidth = 0.9, linetype = "dashed") +
  facet_wrap(~Year) +
  ggtitle("Trawl gear fishery Male") +
  theme_bw()

ggplot(LF_df %>% dplyr::filter(label == "EM", observation == "Trawl gear fishery", Sex == "female"), aes(x = Length)) +
  geom_point(aes(y = Observed, col = "Observed")) +
  geom_line(aes(y = Predicted, col = "Predicted"), linewidth = 0.9, linetype = "dashed") +
  facet_wrap(~Year) +
  ggtitle("Trawl gear fishery Female") +
  theme_bw()

ggplot(LF_df %>% dplyr::filter(label == "EM", observation == "Domestic LL survey", Sex == "male"), aes(x = Length)) +
  geom_point(aes(y = Observed, col = "Observed")) +
  geom_line(aes(y = Predicted, col = "Predicted"), linewidth = 0.9, linetype = "dashed") +
  facet_wrap(~Year) +
  ggtitle("Domestic LL survey Male") +
  theme_bw()

ggplot(LF_df %>% dplyr::filter(label == "EM", observation == "Domestic LL survey", Sex == "female"), aes(x = Length)) +
  geom_point(aes(y = Observed, col = "Observed")) +
  geom_line(aes(y = Predicted, col = "Predicted"), linewidth = 0.9, linetype = "dashed") +
  facet_wrap(~Year) +
  ggtitle("Domestic LL survey Female") +
  theme_bw()

