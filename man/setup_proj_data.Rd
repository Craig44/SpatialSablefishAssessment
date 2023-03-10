% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ProjectionFunction.R
\name{setup_proj_data}
\alias{setup_proj_data}
\title{setup_proj_data}
\usage{
setup_proj_data(mle_obj, n_proj_years = 100, future_recruitment = 2)
}
\arguments{
\item{mle_obj}{an optimised obj that has been created by \code{TMB::MakeADFun}}

\item{n_proj_years}{number of years you want to project out to}

\item{future_recruitment}{the type of future recruitment during projection phase}
}
\value{
data list that can be used by \code{TMB::MakeADFun} to make a projection model.
}
\description{
setup_proj_data
}
\details{
this function will extend all the biological objects and selectivity objects so that they have the correct dimensions.
This is done by repeating the value at n_years or dimension for all projection periods
this function will return a data list that you can further modify if you wish.
}
