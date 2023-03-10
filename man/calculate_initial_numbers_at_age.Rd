% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PartitionFunctions.R
\name{calculate_initial_numbers_at_age}
\alias{calculate_initial_numbers_at_age}
\title{calculate_initial_numbers_at_age}
\usage{
calculate_initial_numbers_at_age(
  n_regions,
  n_ages,
  R0,
  movement_matrix,
  natural_mortality
)
}
\arguments{
\item{n_regions}{integer number of regions}

\item{n_ages}{integer number of ages}

\item{R0}{vector of R0 parameters for each region}

\item{movement_matrix}{movement matrix containing annual movement rates}

\item{natural_mortality}{vector of natural mortality rates for each age}
}
\value{
matrix of numbers at age with regions being rows and ages being the cols
}
\description{
calculate_initial_numbers_at_age
}
