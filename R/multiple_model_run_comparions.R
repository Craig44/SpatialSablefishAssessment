#' get_multiple_ssbs
#'
#' @param mle_ls list with multiple obj$report() calls
#' @param run_labels vector of strings that are labels for each element in mle_ls
#' @param region_key data.frame with colnames area and TMB_ndx for providing real region names to objects
#' @param depletion boolean if true we will scale values by Bzero
#' @return data frame with SSBs catch fits
#' @export

get_multiple_ssbs <- function(mle_ls, run_labels = NULL, region_key = NULL, depletion = T) {
  if(!is.null(run_labels)) {
    if(length(run_labels) != length(mle_ls))
      stop(paste0("Number of models provided ", length(mle_ls), ", number of run labels ", length(run_labels), " these need to be the same"))
  }
  full_ssb_df = NULL
  for(i in 1:length(mle_ls)) {
    this_ssb = get_SSB(MLE_report = mle_ls[[i]], region_key = region_key, depletion = depletion)
    if(!is.null(run_labels)) {
      this_ssb$label = run_labels[i]
    } else {
      this_ssb$label = i
    }
    full_ssb_df = rbind(full_ssb_df, this_ssb)
  }
  full_ssb_df$label = factor(full_ssb_df$label)
  return(full_ssb_df)
}

#' get_multiple_recruits
#'
#' @param mle_ls list with multiple obj$report() calls
#' @param run_labels vector of strings that are labels for each element in mle_ls
#' @param region_key data.frame with colnames area and TMB_ndx for providing real region names to objects
#' @return data frame with SSBs catch fits
#' @export

get_multiple_recruits <- function(mle_ls, run_labels = NULL, region_key = NULL) {
  if(!is.null(run_labels)) {
    if(length(run_labels) != length(mle_ls))
      stop(paste0("Number of models provided ", length(mle_ls), ", number of run labels ", length(run_labels), " these need to be the same"))
  }
  full_recruit_df = NULL
  for(i in 1:length(mle_ls)) {
    if(is.null(mle_ls[[i]])) {
      cat("report at element ", i, " was null, so skipping\n")
      next;
    }
    this_recruit = get_recruitment(MLE_report = mle_ls[[i]], region_key = region_key)
    if(!is.null(run_labels)) {
      this_recruit$label = run_labels[i]
    } else {
      this_recruit$label = i
    }
    full_recruit_df = rbind(full_recruit_df, this_recruit)
  }
  full_recruit_df$label = factor(full_recruit_df$label)
  return(full_recruit_df)
}

#' get_multiple_catch_fits
#'
#' @param mle_ls list with multiple obj$report() calls
#' @param run_labels vector of strings that are labels for each element in mle_ls
#' @param region_key data.frame with colnames area and TMB_ndx for providing real region names to objects
#' @return data frame with multiple catch fits
#' @export
get_multiple_catch_fits <- function(mle_ls, run_labels = NULL, region_key = NULL) {
  if(!is.null(run_labels)) {
    if(length(run_labels) != length(mle_ls))
      stop(paste0("Number of models provided ", length(mle_ls), ", number of run labels ", length(run_labels), " these need to be the same"))
  }
  full_catch_df = NULL
  for(i in 1:length(mle_ls)) {
    if(is.null(mle_ls[[i]])) {
      cat("report at element ", i, " was null, so skipping\n")
      next;
    }
    this_catch = get_catches(MLE_report = mle_ls[[i]], region_key = region_key)
    if(!is.null(run_labels)) {
      this_catch$label = run_labels[i]
    } else {
      this_catch$label = i
    }
    full_catch_df = rbind(full_catch_df, this_catch)
  }
  full_catch_df$label = factor(full_catch_df$label)
  return(full_catch_df)
}


#' get_multiple_index_fits
#'
#' @param mle_ls list with multiple obj$report() calls
#' @param run_labels vector of strings that are labels for each element in mle_ls
#' @param region_key data.frame with colnames area and TMB_ndx for providing real region names to objects
#' @return data frame with multiple catch fits
#' @export
get_multiple_index_fits <- function(mle_ls, run_labels = NULL, region_key = NULL) {
  if(!is.null(run_labels)) {
    if(length(run_labels) != length(mle_ls))
      stop(paste0("Number of models provided ", length(mle_ls), ", number of run labels ", length(run_labels), " these need to be the same"))
  }
  full_index_df = NULL
  for(i in 1:length(mle_ls)) {
    if(is.null(mle_ls[[i]])) {
      cat("report at element ", i, " was null, so skipping\n")
      next;
    }
    this_catch = get_index(MLE_report = mle_ls[[i]], region_key = region_key)
    if(!is.null(run_labels)) {
      this_catch$label = run_labels[i]
    } else {
      this_catch$label = i
    }
    full_index_df = rbind(full_index_df, this_catch)
  }
  full_index_df$label = factor(full_index_df$label)
  return(full_index_df)
}


#' get_multiple_mean_age_fits
#'
#' @param mle_ls list with multiple obj$report() calls
#' @param run_labels vector of strings that are labels for each element in mle_ls
#' @param region_key data.frame with colnames area and TMB_ndx for providing real region names to objects
#' @return data frame with multiple catch fits
#' @export
get_multiple_mean_age_fits <- function(mle_ls, run_labels = NULL, region_key = NULL) {
  if(!is.null(run_labels)) {
    if(length(run_labels) != length(mle_ls))
      stop(paste0("Number of models provided ", length(mle_ls), ", number of run labels ", length(run_labels), " these need to be the same"))
  }
  full_mean_age_df = NULL
  for(i in 1:length(mle_ls)) {
    if(is.null(mle_ls[[i]])) {
      cat("report at element ", i, " was null, so skipping\n")
      next;
    }
    this_mean_age = get_mean_age(MLE_report = mle_ls[[i]], observation = "all", sex = "both", subset_years = NULL, region_key = region_key)
    if(!is.null(run_labels)) {
      this_mean_age$label = run_labels[i]
    } else {
      this_mean_age$label = i
    }
    full_mean_age_df = rbind(full_mean_age_df, this_mean_age)
  }
  full_mean_age_df$label = factor(full_mean_age_df$label)
  return(full_mean_age_df)
}

#' get_multiple_mean_length_fits
#'
#' @param mle_ls list with multiple obj$report() calls
#' @param run_labels vector of strings that are labels for each element in mle_ls
#' @param region_key data.frame with colnames area and TMB_ndx for providing real region names to objects
#' @return data frame with multiple catch fits
#' @export
get_multiple_mean_length_fits <- function(mle_ls, run_labels = NULL, region_key = NULL) {
  if(!is.null(run_labels)) {
    if(length(run_labels) != length(mle_ls))
      stop(paste0("Number of models provided ", length(mle_ls), ", number of run labels ", length(run_labels), " these need to be the same"))
  }
  full_mean_len_df = NULL
  for(i in 1:length(mle_ls)) {
    if(is.null(mle_ls[[i]])) {
      cat("report at element ", i, " was null, so skipping\n")
      next;
    }
    this_mean_len = get_mean_length(MLE_report = mle_ls[[i]], observation = "all", sex = "both", subset_years = NULL, region_key = region_key)
    if(!is.null(run_labels)) {
      this_mean_len$label = run_labels[i]
    } else {
      this_mean_len$label = i
    }
    full_mean_len_df = rbind(full_mean_len_df, this_mean_len)
  }
  full_mean_len_df$label = factor(full_mean_len_df$label)
  return(full_mean_len_df)
}
#' get_multiple_Fs
#'
#' @param mle_ls list with multiple obj$report() calls
#' @param run_labels vector of strings that are labels for each element in mle_ls
#' @param region_key data.frame with colnames area and TMB_ndx for providing real region names to objects
#' @return data frame with SSBs catch fits
#' @export

get_multiple_Fs <- function(mle_ls, run_labels = NULL, region_key = NULL) {
  if(!is.null(run_labels)) {
    if(length(run_labels) != length(mle_ls))
      stop(paste0("Number of models provided ", length(mle_ls), ", number of run labels ", length(run_labels), " these need to be the same"))
  }
  full_Fs_df = NULL
  for(i in 1:length(mle_ls)) {
    if(is.null(mle_ls[[i]])) {
      cat("report at element ", i, " was null, so skipping\n")
      next;
    }
    this_F = get_fishing_mortalities(MLE_report = mle_ls[[i]], region_key = region_key)
    if(!is.null(run_labels)) {
      this_F$label = run_labels[i]
    } else {
      this_F$label = i
    }
    full_Fs_df = rbind(full_Fs_df, this_F)
  }
  full_Fs_df$label = factor(full_Fs_df$label)
  return(full_Fs_df)
}


#' get_multiple_movements
#'
#' @param mle_ls list with multiple obj$report() calls
#' @param run_labels vector of strings that are labels for each element in mle_ls
#' @param region_key data.frame with colnames area and TMB_ndx for providing real region names to objects
#' @return data frame with SSBs catch fits
#' @export

get_multiple_movements <- function(mle_ls, run_labels = NULL, region_key = NULL) {
  if(!is.null(run_labels)) {
    if(length(run_labels) != length(mle_ls))
      stop(paste0("Number of models provided ", length(mle_ls), ", number of run labels ", length(run_labels), " these need to be the same"))
  }
  full_move_df = NULL
  for(i in 1:length(mle_ls)) {
    if(is.null(mle_ls[[i]])) {
      cat("report at element ", i, " was null, so skipping\n")
      next;
    }
    this_move = get_movement(MLE_report = mle_ls[[i]], region_key = region_key)
    if(!is.null(run_labels)) {
      this_move$label = run_labels[i]
    } else {
      this_move$label = i
    }
    full_move_df = rbind(full_move_df, this_move)
  }
  full_move_df$label = factor(full_move_df$label)
  return(full_move_df)
}

#' get_multiple_selectivities
#'
#' @param mle_ls list with multiple obj$report() calls
#' @param run_labels vector of strings that are labels for each element in mle_ls
#' @param region_key data.frame with colnames area and TMB_ndx for providing real region names to objects
#' @return data frame with SSBs catch fits
#' @export

get_multiple_selectivities <- function(mle_ls, run_labels = NULL) {
  if(!is.null(run_labels)) {
    if(length(run_labels) != length(mle_ls))
      stop(paste0("Number of models provided ", length(mle_ls), ", number of run labels ", length(run_labels), " these need to be the same"))
  }
  full_sel_df = NULL
  for(i in 1:length(mle_ls)) {
    if(is.null(mle_ls[[i]])) {
      cat("report at element ", i, " was null, so skipping\n")
      next;
    }
    this_sel = get_selectivities(MLE_report = mle_ls[[i]])
    if(!is.null(run_labels)) {
      this_sel$label = run_labels[i]
    } else {
      this_sel$label = i
    }
    full_sel_df = rbind(full_sel_df, this_sel)
  }
  full_sel_df$label = factor(full_sel_df$label)
  return(full_sel_df)
}


#' get_multiple_nlls accessors to get multiple negative log-likelihoods
#'
#' @param mle_ls list with multiple obj$report() calls
#' @param run_labels vector of strings that are labels for each element in mle_ls
#' @param region_key data.frame with colnames area and TMB_ndx for providing real region names to objects
#' @return data frame with multiple negative log-likelihoods
#' @export
get_multiple_nlls <- function(mle_ls, run_labels = NULL, region_key = NULL) {
  if(!is.null(run_labels)) {
    if(length(run_labels) != length(mle_ls))
      stop(paste0("Number of models provided ", length(mle_ls), ", number of run labels ", length(run_labels), " these need to be the same"))
  }
  full_nll_df = NULL
  for(i in 1:length(mle_ls)) {
    if(is.null(mle_ls[[i]])) {
      cat("report at element ", i, " was null, so skipping\n")
      next;
    }
    this_nll = get_negloglike(MLE_report = mle_ls[[i]])
    ## add totol
    this_nll = rbind(this_nll, data.frame(negloglike = round(sum(this_nll$negloglike),4), observations = "Total"))
    ## addd label
    if(!is.null(run_labels)) {
      this_nll$label = run_labels[i]
    } else {
      this_nll$label = i
    }
    full_nll_df = rbind(full_nll_df, this_nll)
  }
  full_nll_df$label = factor(full_nll_df$label)
  return(full_nll_df)
}
