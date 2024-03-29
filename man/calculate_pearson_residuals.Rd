% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ObservationFunctions.R
\name{calculate_pearson_residuals}
\alias{calculate_pearson_residuals}
\title{calculate_pearson_residuals}
\usage{
calculate_pearson_residuals(MLE_report, region_key = NULL)
}
\arguments{
\item{MLE_report}{a list that is output from obj$report() usually once an optimsation routine has been done.}

\item{region_key}{data.frame with colnames area and TMB_ndx for providing real region names to objects}
}
\value{
data.frame of pearson residuals for each observation
}
\description{
calculate_pearson_residuals
}
\details{
Take a model fit and calculate Pearson resiudals for all observations.
}
