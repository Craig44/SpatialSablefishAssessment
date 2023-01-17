#' get_multiple_ssbs
#'
#' @param mle_ls list with multiple obj$report() calls
#' @param run_labels vector of strings that are labels for each element in mle_ls
#' @return data frame with SSBs catch fits
#' @export

get_multiple_ssbs <- function(mle_ls, run_labels = NULL) {
  if(!is.null(run_labels)) {
    if(length(run_labels) != length(mle_ls))
      stop(paste0("Number of models provided ", length(mle_ls), ", number of run labels ", length(run_labels), " these need to be the same"))
  }
  full_ssb_df = NULL
  for(i in 1:length(mle_ls)) {
    this_ssb = get_SSB(MLE_report = mle_ls[[i]], region_key = region_key)
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

#' get_multiple_catch_fits
#'
#' @param mle_ls list with multiple obj$report() calls
#' @param run_labels vector of strings that are labels for each element in mle_ls
#' @return data frame with multiple catch fits
#' @export
get_multiple_catch_fits <- function(mle_ls, run_labels = NULL) {
  if(!is.null(run_labels)) {
    if(length(run_labels) != length(mle_ls))
      stop(paste0("Number of models provided ", length(mle_ls), ", number of run labels ", length(run_labels), " these need to be the same"))
  }
  full_catch_df = NULL
  for(i in 1:length(mle_ls)) {
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

#' get_multiple_nlls accessors to get multiple negative log-likelihoods
#'
#' @param mle_ls list with multiple obj$report() calls
#' @param run_labels vector of strings that are labels for each element in mle_ls
#' @return data frame with multiple negative log-likelihoods
#' @export
get_multiple_nlls <- function(mle_ls, run_labels = NULL) {
  if(!is.null(run_labels)) {
    if(length(run_labels) != length(mle_ls))
      stop(paste0("Number of models provided ", length(mle_ls), ", number of run labels ", length(run_labels), " these need to be the same"))
  }
  full_nll_df = NULL
  for(i in 1:length(mle_ls)) {
    if(is.null(mle_ls[[i]]))
      next;
    this_nll = get_negloglike(MLE_report = mle_ls[[i]])
    if(!is.null(run_labels)) {
      this_nll$label = run_labels[i]
    } else {
      this_nll$label = i
    }
    this_nll = rbind(this_nll, c(sum(this_nll$negloglike), observations = "Total"))
    full_nll_df = rbind(full_nll_df, this_nll)
  }
  full_nll_df$label = factor(full_nll_df$label)
  return(full_nll_df)
}
