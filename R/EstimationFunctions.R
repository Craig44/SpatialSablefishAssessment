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
#' @return character of good news or labels of problem parameters
#'
check_gradients = function(obj) {
  if(sum(obj$gr() == 0) == 0) {
    return("no gradients were equal to zero. Your model passes the first 'sniff test'. You may continue")
  } else {
    return(names(obj$par)[which(obj$gr() == 0)])
  }
  return(NULL)
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

#' contains a vector of elements that we want to exclude from estimation.
#' @return a list of factors used in the MakeADFun function
#' @export
fix_pars <- function(par_list, pars_to_exclude, vec_elements_to_exclude = NULL, array_elements_to_exclude = NULL) {
  if (!any(pars_to_exclude %in% names(par_list))) {
    stop(paste0("The parameters ", paste(pars_to_exclude[!pars_to_exclude %in% names(par_list)],collapse = " ")," in exclusion parameters could not be found in the 'par_list', please sort this out"))
  }
  pars = names(par_list)
  mapped_pars = list();
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
