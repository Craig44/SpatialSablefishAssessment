% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TagFunctions.R
\name{plot_tag_recovery_obs}
\alias{plot_tag_recovery_obs}
\title{plot_tag_recovery_obs plot predicted and observed tag-recovery obs}
\usage{
plot_tag_recovery_obs(
  MLE_report,
  region_key = NULL,
  release_ndx_to_plot = 1:5,
  sex = "both"
)
}
\arguments{
\item{MLE_report}{a list that is output from obj$report() usually once an optimsation routine has been done.}

\item{region_key}{data.frame with colnames area and TMB_ndx for providing real region names to objects}

\item{release_ndx_to_plot}{vector of integers to create subset plots}
}
\value{
ggplot2
}
\description{
plot_tag_recovery_obs plot predicted and observed tag-recovery obs
}
