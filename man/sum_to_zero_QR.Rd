% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/UtilityFunctions.R
\name{sum_to_zero_QR}
\alias{sum_to_zero_QR}
\title{sum_to_zero_QR
take a vector of unconstrained values length (N - 1) and derive a vector of length N that sum = 0 using the QR method
see here https://discourse.mc-stan.org/t/test-soft-vs-hard-sum-to-zero-constrain-choosing-the-right-prior-for-soft-constrain/3884}
\usage{
sum_to_zero_QR(x_raw)
}
\arguments{
\item{x_raw}{vector of unconstrained values length N - 1}
}
\value{
a vector length N that sums = 0
}
\description{
sum_to_zero_QR
take a vector of unconstrained values length (N - 1) and derive a vector of length N that sum = 0 using the QR method
see here https://discourse.mc-stan.org/t/test-soft-vs-hard-sum-to-zero-constrain-choosing-the-right-prior-for-soft-constrain/3884
}
