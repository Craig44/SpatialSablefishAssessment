% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DerivedQuantityFunctions.R
\name{get_catches}
\alias{get_catches}
\title{get_catches}
\usage{
get_catches(MLE_report, region_key = NULL)
}
\arguments{
\item{MLE_report}{a list that is output from obj$report() usually once an optimsation routine has been done.}

\item{region_key}{data.frame with colnames area and TMB_ndx for providing real region names to objects}
}
\value{
ggplot2 object that will plot if an observation occurs in a year and region
}
\description{
get_catches
}
