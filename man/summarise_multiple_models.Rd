% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BookdownFunctions.R
\name{summarise_multiple_models}
\alias{summarise_multiple_models}
\title{summarise_multiple_models
A bookdown function that will summaries a suite of models in individual Sections of a single Rmarkdown file}
\usage{
summarise_multiple_models(
  model_dir,
  bookdown_dir,
  bookdown_labels,
  model_description = ""
)
}
\arguments{
\item{model_dir}{a vector of absolute paths were we can find model output to summarise}

\item{bookdown_dir}{an absolute path where we will build the Bookdown}

\item{bookdown_labels}{a vector of model labels}

\item{model_description}{A string must be single quotes '' that describes all the models preferably in bullet points that are pointed at by model_dir}
}
\value{
Creates a suite of Rmarkdown files in bookdown_dir and compiles an R markdown file
}
\description{
summarise_multiple_models
A bookdown function that will summaries a suite of models in individual Sections of a single Rmarkdown file
}
\details{
this function expects the following RDS objects to be in each \code{model_dir} for each model, data.RDS region_key.RDS, parameters.RDS, mle_report.RDS, sd_report.RDS, mle_optim.RDS, map_fixed_pars.RDS

Start comparing model fits

Start comparing model quantities
}
