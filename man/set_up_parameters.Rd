% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/EstimationFunctions.R
\name{set_up_parameters}
\alias{set_up_parameters}
\title{set_up_parameters utility function to help 'turn off' parameters and share estimated parameters across elements}
\usage{
set_up_parameters(
  data,
  parameters,
  na_map = NULL,
  srv_sel_first_param_shared_by_sex = F,
  srv_sel_second_param_shared_by_sex = F,
  srv_sel_third_param_shared_by_sex = F,
  fixed_sel_first_shared_by_sex = F,
  fixed_sel_second_shared_by_sex = F,
  fixed_sel_third_shared_by_sex = F,
  trwl_sel_first_shared_by_sex = F,
  trwl_sel_second_shared_by_sex = F,
  trwl_sel_third_shared_by_sex = F,
  recruit_dev_years_not_to_estimate = NULL,
  srv_q_spatial = F,
  tag_reporting_rate = "constant",
  est_init_F = F,
  est_catch_sd = F,
  est_movement = T,
  est_sigma_R = F,
  est_sigma_init_dev = F,
  est_fixed_AF_theta = F,
  est_fixed_LF_theta = F,
  est_trwl_LF_theta = F,
  est_srv_AF_theta = F,
  est_prop_male_recruit = "off",
  est_SR_pars = F,
  common_survey_selex = NULL
)
}
\arguments{
\item{data}{a list of data inputs for the model}

\item{parameters}{a list of parameter values for the model}

\item{na_map}{an existing map that has already had parameters turned off, not well tested}

\item{srv_sel_first_param_shared_by_sex}{vector of bool specifying whether the first survey selectivity parameter is shared among male and female for all time-blocks. There needs to be a bool for each of the n_surveys}

\item{srv_sel_second_param_shared_by_sex}{vector of bool specifying whether the second survey selectivity parameter is shared among male and female for all time-blocks. There needs to be a bool for each of the n_surveys}

\item{srv_sel_third_param_shared_by_sex}{bool, whether the third survey selectivity parameter is shared among male and female for all time-blocks}

\item{fixed_sel_first_shared_by_sex}{bool, whether the first fixed gear selectivity parameter is shared among male and female for all time-blocks}

\item{fixed_sel_second_shared_by_sex}{bool, whether the second fixed gear selectivity parameter is shared among male and female for all time-blocks}

\item{fixed_sel_third_shared_by_sex}{bool, whether the third fixed gear selectivity parameter is shared among male and female for all time-blocks}

\item{trwl_sel_first_shared_by_sex}{bool, whether the first trawl gear selectivity parameter is shared among male and female for all time-blocks}

\item{trwl_sel_second_shared_by_sex}{bool, whether the second trawl gear selectivity parameter is shared among male and female for all time-blocks}

\item{trwl_sel_third_shared_by_sex}{bool, whether the third trawl gear selectivity parameter is shared among male and female for all time-blocks}

\item{recruit_dev_years_not_to_estimate}{vector of years that are not estimated. if \code{data$rec_devs_sum_to_zero} don't include the last year as this assumes n-1 parameters}

\item{srv_q_spatial}{whether regional Q's exist}

\item{tag_reporting_rate}{how to deal with tag-reporting rate will be ignored if there are no tag-recovery observations.
\itemize{
\item \code{off}: not estimated
\item \code{ignore}: ignore this parameter
\item \code{constant}: single value for all years and regions
\item  numeric values indicate the start of a time-block: estimates a tag-reporting rate coefficient for all recovery years which is common across all regions
\item \code{space}: TODO - not implemented estimates an regional tag-reporting rate which is common across all recovery years
\item \code{spatio-temporal}: TODO - not implemented - estimates annual and regional tag-reporting rates
}}

\item{est_init_F}{bool whether you want to estimate an initial F}

\item{est_catch_sd}{bool whether you want to estimate the catch sd parameter}

\item{est_movement}{bool whether you want to estimate movement parameters}

\item{est_sigma_R}{bool whether you want to estimate the recruitment sd parameter}

\item{est_sigma_init_dev}{bool whether you want to estimate the initial-dev sd parameter}

\item{est_fixed_AF_theta}{bool whether you want to estimate the theta parameter for fixed gear AF}

\item{est_fixed_LF_theta}{bool whether you want to estimate the theta parameter for fixed gear LF}

\item{est_trwl_LF_theta}{bool whether you want to estimate the theta parameter for trawl gear LF}

\item{est_srv_AF_theta}{vector of bool specifying whether you want to estimate the theta parameter for survey AF likelihood. length is n_suveys}

\item{est_prop_male_recruit}{vector of years that indicate time-blocks or one of the following strings
\itemize{
\item \code{off}: not estimated
\item \code{constant}: single value for all years
}}

\item{est_SR_pars}{if data$SrType == 2 you can estimate the steepness parameter or have it fixed}

\item{common_survey_selex}{a list each element of the list contains a vector that indicates the survey index that have a common selectivity}
}
\value{
a named list that can be used by the \code{map} input for the \code{TMB::MakeADFun} function. NAs mean parameters are not estimated and elements with the same factor level mean they will be estimated with a common coefficient i.e., shared parameters
}
\description{
set_up_parameters utility function to help 'turn off' parameters and share estimated parameters across elements
}
