% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ObservationFunctions.R
\name{get_AF}
\alias{get_AF}
\title{get_AF accessor function to get age-frequency data}
\usage{
get_AF(
  MLE_report,
  observation = "fixed",
  subset_years = NULL,
  sex = "both",
  region_key = NULL
)
}
\arguments{
\item{MLE_report}{a list that is output from obj$report() usually once an optimsation routine has been done.}

\item{observation}{character labeling the observation you want to plot. See below for options
\itemize{
\item \code{fixed}
\item \code{srv_dom_ll}
}}

\item{subset_years}{vector of years to plot it for}

\item{sex}{character that allows users to specify if the want sex specific plots
\itemize{
\item \code{both}
\item \code{male}
\item \code{female}
}}

\item{region_key}{data.frame with colnames area and TMB_ndx for providing real region names to objects}
}
\value{
data frame with age-frequency info
}
\description{
get_AF accessor function to get age-frequency data
}
