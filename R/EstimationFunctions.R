#
# This script contains functions used to turn estimated paraemters off, and also set parameters to be the same during estimation.
#


#' get_tmb_fixed_effects
#' @description return MLE estaimtes for fixed effects
#' @param obj An optimised list that has been build by MakeAdFun
#' @export
get_tmb_fixed_effects <- function(obj) {
  if (length(obj$env$random) == 0) {
    return(obj$env$last.par.best)
  }
  return(obj$env$last.par.best[-obj$env$random])
}
#' check_gradients
#' @details checks a TMB object for fixed effect parameters that have 0 gradients, suggesting they don't contribute the log likelihood. And you should
#' look into these parameters.
#' @param obj A TMB list that has been built by MakeAdFun
#' @export
#' @return bool and print a message. True is gradients all non zero false suggests there are parameters with zero gradients which is a problem
#'
check_gradients = function(obj) {
  if(sum(obj$gr() == 0) == 0) {
    cat("no parameters had gradients equal to zero.\n")
    return(TRUE)
  } else {
    cat("Found the following parameters with zero gradients:\n\n")
    cat(paste(paste0(names(obj$par)[which(obj$gr() == 0)], ", index = ", which(obj$gr() == 0)), collapse = "\n"))
    return(FALSE)
  }
  return(TRUE)
}

#' rmvnorm_prec
#' @description simulates parameters from the joint precision matrix derived from a TMB objects
#' @param mu vector of MLE both fixed and random effect parameters
#' @param prec precision matrix, derived from sdreport(obj, getJointPrecision = T)
#' @param n.sims integer number of simulations
#' @param random_seed integer seed
#' @export
#' @importFrom Matrix Cholesky solve
rmvnorm_prec <- function(mu, prec, n.sims, random_seed ) {
  set.seed( random_seed )
  z = matrix(stats::rnorm(length(mu) * n.sims), ncol=n.sims)
  L = Cholesky(prec, super=TRUE)
  z = solve(L, z, system = "Lt") ## z = Lt^-1 %*% z
  z = solve(L, z, system = "Pt") ## z = Pt    %*% z
  z = as.matrix(z)
  return(mu + z)
}

#' fix_pars
#' @author C.Marsh
#' @description TMB helper function this function returns a list of factors used in the map argument of the MakeADFun function
#' values with NA will not be estimated.
#' @param par_list a named list that you give to the par argument in the MakeADFun
#' @param pars_to_exclude a vector of strings with names of parameters you want to FIX in the objective object.
#' @param vec_elements_to_exclude a named list (names %in% pars_to_exclude) with number of elements = length(vec_pars_to_adjust). each list element
#' @param array_elements_to_exclude a named list (names %in% pars_to_exclude) with a matrix each row corresponds to an element with the first column being the array row index and second column being the array column index to fix
#' @param existing_map a named list that already contains NAs from previous fix_pars calls
#' contains a vector of elements that we want to exclude from estimation.
#' @return a list of factors used in the MakeADFun function
#' @export
fix_pars <- function(par_list, pars_to_exclude, vec_elements_to_exclude = NULL, array_elements_to_exclude = NULL, existing_map = NULL) {
  if (!any(pars_to_exclude %in% names(par_list))) {
    stop(paste0("The parameters ", paste(pars_to_exclude[!pars_to_exclude %in% names(par_list)],collapse = " ")," in exclusion parameters could not be found in the 'par_list', please sort this out"))
  }
  pars = names(par_list)
  mapped_pars = list();
  if(!is.null(existing_map))
    mapped_pars = existing_map

  if (!is.null(vec_elements_to_exclude)) {
    if (!all(names(vec_elements_to_exclude) %in% pars_to_exclude))
      stop("parameters names in vec_elements_to_exclude, need to also be in pars_to_exclude")
  }
  if (!is.null(array_elements_to_exclude)) {
    if (!all(names(array_elements_to_exclude) %in% pars_to_exclude))
      stop("parameters names in array_elements_to_exclude, need to also be in pars_to_exclude")
  }
  param_factor = 1;
  for(i in 1:length(pars)) {
    if (pars[i] %in% pars_to_exclude) {
      params_in_this_par = par_list[[pars[i]]];
      if (pars[i] %in% names(vec_elements_to_exclude)) {
        include_element_index = c(1:length(params_in_this_par))[-vec_elements_to_exclude[[pars[i]]]]
        params_vals = factor(rep(NA, length(params_in_this_par)), levels = factor(param_factor:(param_factor + length(include_element_index) - 1)))
        params_vals[include_element_index] = factor(param_factor:(param_factor + length(include_element_index) - 1))#, levels = factor(include_element_index))
        param_factor = param_factor + length(include_element_index)
        mapped_pars[[pars[i]]] = params_vals;
      } else if(pars[i] %in% names(array_elements_to_exclude)) {
        elements_to_drop = array_elements_to_exclude[[pars[i]]]
        mapped_vector = rep(NA, length(params_in_this_par))
        first_param_factor = param_factor
        vec_ndx = 1;
        ## TMB converts arrays to vectors down columns (not by rows)
        ## can handle up to 3-dimension arrays
        if(length(dim(params_in_this_par)) == 2) {
          for(col_ndx in 1:ncol(params_in_this_par)) {
            for(row_ndx in 1:nrow(params_in_this_par)) {
              dropping_this_element = F
              for(drop_ndx in 1:nrow(elements_to_drop)) {
                if(all(c(row_ndx, col_ndx) == elements_to_drop[drop_ndx,])) {
                  dropping_this_element = T
                  break;
                }
              }
              if(!dropping_this_element) {
                mapped_vector[vec_ndx] = param_factor
                param_factor = param_factor + 1
              }
              vec_ndx = vec_ndx + 1;
            }
          }
        } else if (length(dim(params_in_this_par)) == 3) {
          counter = 1;
          for(dim3_ndx in 1:dim(params_in_this_par)[3]) {
            for(dim2_ndx in 1:dim(params_in_this_par)[2]) {
              for(dim1_ndx in 1:dim(params_in_this_par)[1]) {
                ## check if we need to drop this value
                dropping_this_element = F
                for(drop_ndx in 1:nrow(elements_to_drop)) {
                  if(all(c(dim1_ndx, dim2_ndx, dim3_ndx) == elements_to_drop[drop_ndx,])) {
                    dropping_this_element = T
                    break;
                  }
                }
                if(!dropping_this_element) {
                  mapped_vector[vec_ndx] = param_factor
                  param_factor = param_factor + 1
                }
                vec_ndx = vec_ndx + 1;
              }
            }
          }
        } else {
          stop("this function can only deal with 2 or 3 dimensional arrays")
        }
        mapped_vector = factor(mapped_vector, levels = first_param_factor:max(mapped_vector, na.rm = T))
        mapped_pars[[pars[i]]] = mapped_vector;
      } else {
        ## exclude entire parameters
        mapped_pars[[pars[i]]] = rep(factor(NA),length(params_in_this_par));
        n_params_to_exclude = nrow(vec_elements_to_exclude[[pars[i]]])
      }
    } else {
      params_in_this_par = par_list[[pars[i]]];
      params_vals = factor(param_factor:(param_factor + length(params_in_this_par) - 1))
      param_factor = param_factor + length(params_in_this_par)
      mapped_pars[[pars[i]]] = params_vals
    }
  }
  return(mapped_pars);
}

#' set_up_parameters utility function to help 'turn off' parameters and share estimated parameters across elements
#' @param data a list of data inputs for the model
#' @param parameters a list of parameter values for the model
#' @param na_map an existing map that has already had parameters turned off, not well tested
#' @param srv_sel_first_param_shared_by_sex bool, whether the first survey selectivity parameter is shared among male and female for all time-blocks
#' @param srv_sel_second_param_shared_by_sex bool, whether the second survey selectivity parameter is shared among male and female for all time-blocks
#' @param fixed_sel_first_shared_by_sex bool, whether the first fixed gear selectivity parameter is shared among male and female for all time-blocks
#' @param fixed_sel_second_shared_by_sex bool, whether the second fixed gear selectivity parameter is shared among male and female for all time-blocks
#' @param trwl_sel_first_shared_by_sex bool, whether the first trawl gear selectivity parameter is shared among male and female for all time-blocks
#' @param trwl_sel_second_shared_by_sex bool, whether the second trawl gear selectivity parameter is shared among male and female for all time-blocks
#' @param recruit_dev_years_not_to_estimate vector of years that are not estimated.
#' @param srv_q_spatial whether regional Q's exist
#' @param tag_reporting_rate how to deal with tag-reporting rate will be ignored if there are no tag-recovery observations.
#' \itemize{
#'   \item `off`: not estimated
#'   \item `ignore`: ignore this parameter
#'   \item `constant`: single value for all years and regions
#'   \item `time`: estimates a tag-reporting rate coefficient for all recovery years which is common across all regions
#'   \item `space`: TODO - not implemented estimates an regional tag-reporting rate which is common across all recovery years
#'   \item `spatio-temporal`: TODO - not implemented - estimates annual and regional tag-reporting rates
#' }
#' @param est_init_F bool whether you want to estimate an initial F
#' @param est_catch_sd bool whether you want to estimate the catch sd parameter
#' @param est_movement bool whether you want to estimate movement parameters
#' @param est_sigma_R bool whether you want to estimate the recruitment sd parameter
#' @param est_sigma_init_dev bool whether you want to estimate the initial-dev sd parameter
#' @return a named list that can be used by the `map` input for the `TMB::MakeADFun` function. NAs mean parameters are not estimated and elements with the same factor level mean they will be estimated with a common coefficient i.e., shared parameters
#' @export
set_up_parameters <- function(data, parameters,
                              na_map = NULL,
                              srv_sel_first_param_shared_by_sex = F,
                              srv_sel_second_param_shared_by_sex = F,
                              fixed_sel_first_shared_by_sex = F,
                              fixed_sel_second_shared_by_sex = F,
                              trwl_sel_first_shared_by_sex = F,
                              trwl_sel_second_shared_by_sex = F,
                              recruit_dev_years_not_to_estimate = NULL,
                              srv_q_spatial = F,
                              tag_reporting_rate = "constant",
                              est_init_F = F,
                              est_catch_sd = F,
                              est_movement = T,
                              est_sigma_R = F,
                              est_sigma_init_dev = F
) {

  if(!tag_reporting_rate %in% c("ignore","constant", "time", "off"))#, "space", "spatio-temporal"))
    stop('tag_reporting_rate: needs to be one of the following "ignore", "constant", "time", "off"')
  #
  map_to_fix = na_map
  parameters_completely_fixed = c()
  vectors_with_elements_fixed = list()
  arrays_with_elements_fixed = list()

  # turn off F-avg and devs
  if(data$F_method == 1)
    parameters_completely_fixed = c(parameters_completely_fixed, c("ln_fixed_F_avg", "ln_fixed_F_devs","ln_trwl_F_avg", "ln_trwl_F_devs"))
  # trun off init F if not estimating it
  if(!est_init_F)
    parameters_completely_fixed = c(parameters_completely_fixed, c("ln_init_F_avg"))
  # trun off sd catch if not estimating it
  if(!est_catch_sd)
    parameters_completely_fixed = c(parameters_completely_fixed, c("ln_catch_sd"))
  if(!est_sigma_R)
    parameters_completely_fixed = c(parameters_completely_fixed, c("ln_sigma_R"))
  if(!est_sigma_init_dev)
    parameters_completely_fixed = c(parameters_completely_fixed, c("ln_sigma_init_devs"))

  # turn off movement estimation parameters
  if(!est_movement)
    parameters_completely_fixed = c(parameters_completely_fixed, c("transformed_movement_pars"))
  # turn off init_rec devs if not applying
  if(data$n_init_rec_devs == 0) {
    parameters_completely_fixed = c(parameters_completely_fixed, c("ln_init_rec_dev"))
    if(!"ln_sigma_init_devs" %in% parameters_completely_fixed)
      parameters_completely_fixed = c(parameters_completely_fixed, c("ln_sigma_init_devs"))
  }
  ## survey catchability regional and annual
  base_q_vals = list()
  copy_q_vals = list()
  if(!srv_q_spatial) {
    ## regionally similar q's
    logis_sel_q = list(logistic_srv_dom_ll_q = expand.grid(1:data$n_regions, 1:ncol(parameters$logistic_srv_dom_ll_q))[-1,])
    arrays_with_elements_fixed[["logistic_srv_dom_ll_q"]] = logis_sel_q$logistic_srv_dom_ll_q
    base_q_vals = rep(list(logistic_srv_dom_ll_q = 1), dim(parameters$logistic_srv_dom_ll_q)[1] * dim(parameters$logistic_srv_dom_ll_q)[2] - 1)
    copy_q_vals = evalit(paste0("list(",paste(paste0("logistic_srv_dom_ll_q = ", 2:(dim(parameters$logistic_srv_dom_ll_q)[1] * dim(parameters$logistic_srv_dom_ll_q)[2])), collapse = ", "),")"))
  }
  ## no tag recovery observations
  base_tag_report_vals = list()
  copy_tag_report_vals = list()
  if(sum(data$tag_recovery_indicator) == 0) {
    parameters_completely_fixed = c(parameters_completely_fixed, c("logistic_tag_reporting_rate", "ln_tag_phi"))
  } else {
    if(tag_reporting_rate == "off") {
      parameters_completely_fixed = c(parameters_completely_fixed, c("logistic_tag_reporting_rate"))
    } else if(tag_reporting_rate == "ignore") {

      # don't do anything

    } else if(tag_reporting_rate == "constant") {
      # all regions and years have the same reporting value
      tag_report_fixd_elements = list(logistic_tag_reporting_rate = expand.grid(1:data$n_regions, 1:ncol(parameters$logistic_tag_reporting_rate))[-1,])
      arrays_with_elements_fixed[["logistic_tag_reporting_rate"]] = tag_report_fixd_elements$logistic_tag_reporting_rate

      base_tag_report_vals = rep(list(logistic_tag_reporting_rate = 1), dim(parameters$logistic_tag_reporting_rate)[1] * dim(parameters$logistic_tag_reporting_rate)[2] - 1)
      copy_tag_report_vals = evalit(paste0("list(",paste(paste0("logistic_tag_reporting_rate = ", 2:(dim(parameters$logistic_tag_reporting_rate)[1] * dim(parameters$logistic_tag_reporting_rate)[2])), collapse = ", "),")"))

    } else if (tag_reporting_rate == "time") {
      # all years have the same reporting value
      tag_report_fixd_elements = list(logistic_tag_reporting_rate = expand.grid(2:data$n_regions, 1:ncol(parameters$logistic_tag_reporting_rate)))
      arrays_with_elements_fixed[["logistic_tag_reporting_rate"]] = tag_report_fixd_elements$logistic_tag_reporting_rate
      ## build copy and base parameter labels
      counter = 1;
      for(i in 1:ncol(parameters$logistic_tag_reporting_rate)) {
        this_bse_lst = rep(list(logistic_tag_reporting_rate = counter), data$n_regions - 1)
        base_tag_report_vals = append(base_tag_report_vals, this_bse_lst)
        this_cpy_lst = sapply((counter + 1):(counter + data$n_regions - 1), FUN = function(x) {
          list(logistic_tag_reporting_rate =  x)
        })
        copy_tag_report_vals = append(copy_tag_report_vals, this_cpy_lst)

        counter = counter + data$n_regions

      }
    }
    ## negative binomial we estimate ln_tag_phi else it should be estimated
    if(data$tag_likelihood != 1)
      parameters_completely_fixed = c(parameters_completely_fixed, c("ln_tag_phi"))
  }

  ## deal with recruit devs
  if(!is.null(recruit_dev_years_not_to_estimate)) {
    yr_ndx = which(data$years %in% recruit_dev_years_not_to_estimate)

    if(data$global_rec_devs == 1) {
      ln_rec_dev_elements_to_fix = matrix(c(rep(1, length(yr_ndx)), yr_ndx), byrow = F, ncol = 2)
    } else {
      ln_rec_dev_elements_to_fix = matrix(c(rep(1:data$n_regions, each = length(yr_ndx)), rep(yr_ndx, data$n_regions)), byrow = F, ncol = 2)
    }
    arrays_with_elements_fixed[["ln_rec_dev"]] = ln_rec_dev_elements_to_fix
  }
  ## are we sharing selectivity parameters
  base_srv_sel_param_vals = list()
  copy_srv_sel_param_vals = list()
  if(srv_sel_first_param_shared_by_sex) {
    ## turn off females delta param will be set the same as male
    srv_sel_param_elements_to_fix = cbind(1:dim(parameters$ln_srv_dom_ll_sel_pars)[1], 1,2)

    ## fix male and female for all time-blocks
    counter = 1
    n_time_blocks = dim(parameters$ln_srv_dom_ll_sel_pars)[1]
    for(t_ndx in 1:n_time_blocks) {
      n_sel_pars_for_this_time_block = ifelse(data$srv_dom_ll_sel_type[t_ndx] != 2, 2, 1);
      this_bse_lst = list(ln_srv_dom_ll_sel_pars = counter)
      this_cpy_lst = list(ln_srv_dom_ll_sel_pars = counter + n_time_blocks * n_sel_pars_for_this_time_block)
      base_srv_sel_param_vals = append(base_srv_sel_param_vals, this_bse_lst)
      copy_srv_sel_param_vals = append(copy_srv_sel_param_vals, this_cpy_lst)
      counter = counter + 1
    }
    arrays_with_elements_fixed[["ln_srv_dom_ll_sel_pars"]] = srv_sel_param_elements_to_fix
  }
  if(srv_sel_second_param_shared_by_sex) {
    ## turn off females delta param will be set the same as male
    srv_sel_param_elements_to_fix = cbind(1:dim(parameters$ln_srv_dom_ll_sel_pars)[1], 2,2)

    ## fix male and female for all time-blocks
    n_time_blocks = dim(parameters$ln_srv_dom_ll_sel_pars)[1]
    counter = n_time_blocks + 1

    for(t_ndx in 1:dim(parameters$ln_srv_dom_ll_sel_pars)[1]) {
      n_sel_pars_for_this_time_block = ifelse(data$srv_dom_ll_sel_type[t_ndx] != 2, 2, 1);
      this_bse_lst = list(ln_srv_dom_ll_sel_pars = counter)
      this_cpy_lst = list(ln_srv_dom_ll_sel_pars = counter + n_time_blocks * n_sel_pars_for_this_time_block)
      base_srv_sel_param_vals = append(base_srv_sel_param_vals, this_bse_lst)
      copy_srv_sel_param_vals = append(copy_srv_sel_param_vals, this_cpy_lst)
      counter = counter + 1
    }
    if(is.null(arrays_with_elements_fixed[["ln_srv_dom_ll_sel_pars"]])) {
      arrays_with_elements_fixed[["ln_srv_dom_ll_sel_pars"]] = srv_sel_param_elements_to_fix
    } else {
      arrays_with_elements_fixed[["ln_srv_dom_ll_sel_pars"]] = rbind(arrays_with_elements_fixed[["ln_srv_dom_ll_sel_pars"]], srv_sel_param_elements_to_fix)
    }
  }
  ## fixed gear fishery
  ## are we sharing selectivity parameters
  base_fixed_sel_param_vals = list()
  copy_fixed_sel_param_vals = list()
  if(fixed_sel_first_shared_by_sex) {
    ## turn off females delta param will be set the same as male
    fixed_sel_param_elements_to_fix = cbind(1:dim(parameters$ln_fixed_sel_pars)[1], 1,2)

    ## fix male and female for all time-blocks
    counter = 1
    n_time_blocks = dim(parameters$ln_fixed_sel_pars)[1]
    for(t_ndx in 1:dim(parameters$ln_fixed_sel_pars)[1]) {
      n_sel_pars_for_this_time_block = ifelse(data$fixed_sel_type[t_ndx] != 2, 2, 1);
      this_bse_lst = list(ln_fixed_sel_pars = counter)
      this_cpy_lst = list(ln_fixed_sel_pars = counter + n_time_blocks * n_sel_pars_for_this_time_block)
      base_fixed_sel_param_vals = append(base_fixed_sel_param_vals, this_bse_lst)
      copy_fixed_sel_param_vals = append(copy_fixed_sel_param_vals, this_cpy_lst)
      counter = counter + 1
    }
    arrays_with_elements_fixed[["ln_fixed_sel_pars"]] = fixed_sel_param_elements_to_fix
  }
  if(fixed_sel_second_shared_by_sex) {
    ## turn off females delta param will be set the same as male
    fixed_sel_param_elements_to_fix = cbind(1:dim(parameters$ln_fixed_sel_pars)[1], 2,2)

    ## fix male and female for all time-blocks
    n_time_blocks = dim(parameters$ln_fixed_sel_pars)[1]
    counter = n_time_blocks + 1
    for(t_ndx in 1:dim(parameters$ln_fixed_sel_pars)[1]) {
      n_sel_pars_for_this_time_block = ifelse(data$fixed_sel_type[t_ndx] != 2, 2, 1);
      this_bse_lst = list(ln_fixed_sel_pars = counter)
      this_cpy_lst = list(ln_fixed_sel_pars = counter + n_time_blocks * n_sel_pars_for_this_time_block)
      base_fixed_sel_param_vals = append(base_fixed_sel_param_vals, this_bse_lst)
      copy_fixed_sel_param_vals = append(copy_fixed_sel_param_vals, this_cpy_lst)
      counter = counter + 1
    }
    if(is.null(arrays_with_elements_fixed[["ln_fixed_sel_pars"]])) {
      arrays_with_elements_fixed[["ln_fixed_sel_pars"]] = fixed_sel_param_elements_to_fix
    } else {
      arrays_with_elements_fixed[["ln_fixed_sel_pars"]] = rbind(arrays_with_elements_fixed[["ln_fixed_sel_pars"]], fixed_sel_param_elements_to_fix)
    }
  }
  ## trwl gear fishery
  ## are we sharing selectivity parameters
  base_trwl_sel_param_vals = list()
  copy_trwl_sel_param_vals = list()
  if(trwl_sel_first_shared_by_sex) {
    ## turn off females delta param will be set the same as male
    trwl_sel_param_elements_to_fix = cbind(1:dim(parameters$ln_trwl_sel_pars)[1], 1,2)

    ## fix male and female for all time-blocks
    counter = 1
    n_time_blocks = dim(parameters$ln_trwl_sel_pars)[1]
    for(t_ndx in 1:dim(parameters$ln_trwl_sel_pars)[1]) {
      n_sel_pars_for_this_time_block = ifelse(data$trwl_sel_type[t_ndx] != 2, 2, 1);
      this_bse_lst = list(ln_trwl_sel_pars = counter)
      this_cpy_lst = list(ln_trwl_sel_pars = counter + n_time_blocks * n_sel_pars_for_this_time_block)
      base_trwl_sel_param_vals = append(base_trwl_sel_param_vals, this_bse_lst)
      copy_trwl_sel_param_vals = append(copy_trwl_sel_param_vals, this_cpy_lst)
      counter = counter + 1
    }
    arrays_with_elements_fixed[["ln_trwl_sel_pars"]] = trwl_sel_param_elements_to_fix
  }
  if(trwl_sel_second_shared_by_sex) {
    ## turn off females delta param will be set the same as male
    trwl_sel_param_elements_to_fix = cbind(1:dim(parameters$ln_trwl_sel_pars)[1], 2,2)

    ## fix male and female for all time-blocks
    n_time_blocks = dim(parameters$ln_trwl_sel_pars)[1]
    counter = n_time_blocks + 1
    for(t_ndx in 1:dim(parameters$ln_trwl_sel_pars)[1]) {
      n_sel_pars_for_this_time_block = ifelse(data$trwl_sel_type[t_ndx] != 2, 2, 1);
      this_bse_lst = list(ln_trwl_sel_pars = counter)
      this_cpy_lst = list(ln_trwl_sel_pars = counter + n_time_blocks * n_sel_pars_for_this_time_block)
      base_trwl_sel_param_vals = append(base_trwl_sel_param_vals, this_bse_lst)
      copy_trwl_sel_param_vals = append(copy_trwl_sel_param_vals, this_cpy_lst)
      counter = counter + 1
    }
    if(is.null(arrays_with_elements_fixed[["ln_trwl_sel_pars"]])) {
      arrays_with_elements_fixed[["ln_trwl_sel_pars"]] = trwl_sel_param_elements_to_fix
    } else {
      arrays_with_elements_fixed[["ln_trwl_sel_pars"]] = rbind(arrays_with_elements_fixed[["ln_trwl_sel_pars"]], trwl_sel_param_elements_to_fix)
    }
  }
  ## initial fix pars
  map_to_fix = fix_pars(parameters, pars_to_exclude = c(parameters_completely_fixed, names(arrays_with_elements_fixed), names(vectors_with_elements_fixed)), vec_elements_to_exclude = vectors_with_elements_fixed, array_elements_to_exclude = arrays_with_elements_fixed, existing_map = na_map)
  ## append all base and copy lists elements
  bse_params = append(base_tag_report_vals, base_srv_sel_param_vals)
  bse_params = append(bse_params, base_fixed_sel_param_vals)
  bse_params = append(bse_params, base_trwl_sel_param_vals)
  bse_params = append(bse_params, base_q_vals)

  cpy_params = append(copy_tag_report_vals, copy_srv_sel_param_vals)
  cpy_params = append(cpy_params, copy_fixed_sel_param_vals)
  cpy_params = append(cpy_params, copy_trwl_sel_param_vals)
  cpy_params = append(cpy_params, copy_q_vals)

  if(length(bse_params) > 0) {
    ## some-times this wont be populated which means we won't be fixing any paramaeters
    if(length(bse_params) != length(cpy_params))
      stop("An error occured: base_parameters has different length to copy_parameters for the 'set_pars_to_be_the_same' function")
    map_to_fix = set_pars_to_be_the_same(par_list = parameters, map = map_to_fix, base_parameters = bse_params, copy_parameters = cpy_params)
  }
  return(map_to_fix)
}
#' set_pars_to_be_the_same
#' @author C.Marsh
#' @description TMB helper function this function returns a list of factors used in the map argument of the MakeADFun function
#' values with the same factor level will be estimated as the same value
#' @details TMB will estimate parameters based on the index specified in by the map argument in MakeADFun
#' so parameters with the same factor in map will be estimated as the same value.
#' NOTE: this only works for within the same parameter variable It doesn't work across parameters variables.
#' @param par_list a named list that you give to the par argument in the MakeADFun
#' @param map a list of factors that has been created by fix_pars(). parameters that you want fixed to other values should be set to NA in this object
#' @param base_parameters a named list (names) each element contains one index that will be used to set the value in copy_parameters
#' @param copy_parameters a named list (names) each element contains one index that will be set equal to the corresponding base_parameters
#' @return a list of factors used in the MakeADFun function
#' @export
set_pars_to_be_the_same <- function(par_list, map, base_parameters, copy_parameters) {
  if(length(base_parameters) != length(copy_parameters))
    stop("the number of elements in base_parameters must be the same as copy_parameters. Please check these")
  if(!inherits(map, "list"))
    stop("map needs to be a list")
  if(!any(names(base_parameters) %in% names(par_list)))
    stop(!paste0("The parameters in base_parameters ", paste(names(base_parameters)[!names(base_parameters) %in% names(par_list)],collapse = " ")," could not be found in the 'par_list', please sort this out"))
  if(!any(names(copy_parameters) %in% names(par_list)))
    stop(!paste0("The parameters in copy_parameters ", paste(names(copy_parameters)[!names(copy_parameters) %in% names(par_list)],collapse = " ")," could not be found in the 'par_list', please sort this out"))
  pars = names(par_list)

  for(i in 1:length(base_parameters)) {
    if(is.na(map[[names(base_parameters)[i]]][base_parameters[[i]]]))
      stop(paste0("In base_parameters for parameter ", names(base_parameters)[i], " at ndx ", base_parameters[[i]], ". We found an NA. This cannot be, please check"))
    if(!is.na(map[[names(copy_parameters)[i]]][copy_parameters[[i]]]))
      stop(paste0("In copy_parameters for parameter ", names(base_parameters)[i], " at ndx ", base_parameters[[i]], ". Was not an NA. This must be an NA value in 'map', please check"))

    temp_copy_parm = map[[names(copy_parameters)[i]]]
    temp_copy_parm = as.numeric(as.character(temp_copy_parm))
    base_value = as.numeric(as.character(map[[names(base_parameters)[i]]][base_parameters[[i]]]))
    temp_copy_parm[copy_parameters[[i]]] = base_value
    lvls = unique(temp_copy_parm[!is.na(temp_copy_parm)])
    map[[names(copy_parameters)[i]]] = factor(temp_copy_parm, levels = lvls)
  }

  return(map);
}


#' pre_optim_sanity_checks run a few sanity checks before estimating a model
#' @author C.Marsh
#' @description This function will do some simple checks such check likelihoods are all finite, no non-zero gradients
#' @param obj an compiled MakeADFun TMB object that contains fn(), gr(), and par functions
#' @return bool true means model passed false means model failed. Will also print a message to screen on what failed.
#' @export

pre_optim_sanity_checks <- function(obj) {
  passed_pre_sanity_checks = TRUE;
  ## get the data
  data = obj$env$data
  if(data$model != "TagIntegrated")
    stop("sanity checks are made for the 'TagIntegrated' model")
  ## get derived quantities
  mle_report = obj$report(obj$env$last.par.best)
  ## check likelihoods are all finite and not NA
  if(!all(is.finite(mle_report$nll))) {
    cat("Found Inf in log-likelihood, you will need to resolve this before optimisation\n")
    passed_pre_sanity_checks = F
  }
  if(!all(!is.na(mle_report$nll))) {
    cat("Found NA in log-likelihood, you will need to resolve this before optimisation\n")
    passed_pre_sanity_checks = F
  }
  ## check gradients
  if(!check_gradients(obj))
    passed_pre_sanity_checks = F
  ## return a pass or fail message
  cat("\n\n");
  if(passed_pre_sanity_checks) {
    cat("Successfully passed pre-sanity checks\n")
    return(TRUE)
  }

  ## else must have failed
  cat("Failed pre-sanity checks, please sort these out before you try and optimise this model\n")
  return(FALSE)
}

#' post_optim_sanity_checks run a few sanity checks before estimating a model
#' @author C.Marsh
#' @description This function will do some simple checks such check likelihoods are all finite, no non-zero gradients
#' @param mle_obj which is TMB object that contains fn(), gr(), and par functions. Should have been optimised
#' @param mle_pars vector of fixed effect parameters that can be passed to fn and gr obj calls
#' @return bool true means model passed false means model failed. Will also print a message to screen on what failed.
#' @export

post_optim_sanity_checks <- function(mle_obj, mle_pars, max_abs_gradient = 0.00001) {
  passed_post_sanity_checks = TRUE;
  ## turn obj to silent so messages don't get lost
  mle_obj$env$silent = F
  mle_obj$env$tracepar = F
  mle_obj$env$tracemgc = F

  ## get the data
  data = mle_obj$env$data
  if(data$model != "TagIntegrated")
    stop("sanity checks are made for the 'TagIntegrated' model")
  ## get derived quantities
  mle_report = mle_obj$report(mle_obj$env$last.par.best)
  ## check likelihoods are all finite and not NA
  if(!all(is.finite(mle_report$nll))) {
    cat("Found Inf in log-likelihood, you will need to resolve this before optimisation\n")
    passed_post_sanity_checks = F
  }
  if(!all(!is.na(mle_report$nll))) {
    cat("Found NA in log-likelihood, you will need to resolve this before optimisation\n")
    passed_post_sanity_checks = F
  }
  ## check max absolute gradients
  max_abs_grad_ndx = which.max(abs(mle_obj$gr(mle_pars)))
  max_abs_grad = abs(mle_obj$gr(mle_pars))[max_abs_grad_ndx]
  if(max_abs_gradient < max_abs_grad) {
    cat("Parameter: ", names(mle_obj$par)[max_abs_grad_ndx], " had absolute gradient = ", max_abs_grad, " which was greater than tolerance ", max_abs_gradient,". This indicates non-convergence.\n")
    passed_post_sanity_checks = F
  }

  ## check hessian not singular
  calculate_hessian = tryCatch(expr = optimHess(mle_pars, fn = mle_obj$fn, gr = mle_obj$gr)
                               , error = function(e){e})

  if(inherits(calculate_hessian,"error")) {
    cat("Hessian calculation: ", calculate_hessian$message,"\n")
    passed_post_sanity_checks = F
  } else {
    ## invert hessian
    if(!is_matrix_invertable(calculate_hessian)) {
      cat("Hessian matrix not invertible: probably due to singularity. Check the eigen values to find possible problematic parameters\n")
      passed_post_sanity_checks = F
    }
    ## look at eigen values
    hess_eigen_vals = eigen(calculate_hessian)
    if(any(hess_eigen_vals$values < 0)) {
      cat("Found negative eigen values for the Hessian matrix. This occurs when the matrix is non-zero, real, symmetric, and not positive semi-definite.\n")
      passed_post_sanity_checks = F
    }
    ## just a message not sure how to flag this one. TODO add a threshold
    cat("parameter ", names(mle_obj$par)[which.min(hess_eigen_vals$values)], " had the smallest eigen value of ", hess_eigen_vals$values[which.min(hess_eigen_vals$values)], "\n\n")
  }

  ## check parameters are not at bounds
  ## rec-devs
  if(any(abs(mle_report$recruitment_multipliers) > 10)) {
    cat("Found recruitment_multipliers that had absolute values greater than 10. This is extreme value (possible convergence at bound) and worthy of furthur investigation.\n")
    passed_post_sanity_checks = F
  }
  ## Survey Catchability
  if(any(mle_report$srv_dom_ll_q > 0.95)){
    cat("Found survey catchability greater than 0.95. This is an extreme value (possible convergence at bound)and worthy of furthur investigation.\n")
    passed_post_sanity_checks = F
  }
  if(any(mle_report$srv_dom_ll_q < 0.00001)){
    cat("Found survey catchability less than 0.00001. This is an extreme value (possible convergence at bound) and worthy of furthur investigation.\n")
    passed_post_sanity_checks = F
  }
  ## tag reporting
  if(any(mle_report$tag_reporting_rate > 0.95)){
    cat("Found tag-reporting rate greater than 0.95. This is an extreme value (possible convergence at bound)and worthy of furthur investigation.\n")
    passed_post_sanity_checks = F
  }
  if(any(mle_report$tag_reporting_rate < 0.00001)){
    cat("Found tag-reporting rate less than 0.00001. This is an extreme value (possible convergence at bound) and worthy of furthur investigation.\n")
    passed_post_sanity_checks = F
  }
  ## Selectivity values
  if(any(mle_report$trwl_sel_pars > 40)){
    cat("Found trawl fishery selectivity pars greater than 40 This is an extreme value (possible convergence at bound) and worthy of furthur investigation.\n")
    passed_post_sanity_checks = F
  }
  if(any(mle_report$fixed_sel_pars > 40)){
    cat("Found fixed fishery selectivity pars greater than 40 This is an extreme value (possible convergence at bound) and worthy of furthur investigation.\n")
    passed_post_sanity_checks = F
  }
  if(any(mle_report$srv_dom_ll_sel_pars > 40)){
    cat("Found survey selectivity pars greater than 40 This is an extreme value (possible convergence at bound) and worthy of furthur investigation.\n")
    passed_post_sanity_checks = F
  }

  ## return a pass or fail message
  cat("\n\n");
  if(passed_post_sanity_checks) {
    cat("Successfully passed post-optim-sanity checks\n")
    return(TRUE)
  }

  ## else must have failed
  cat("Failed post-optim-sanity checks. This suggests your model has not non-converged. It is best practice to explore these problems before evaluating model fits and derived quantities.\n")
  return(FALSE)
}
