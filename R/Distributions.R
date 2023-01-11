#' ddirichmult
#' Return the PDF for the dirichlet-multinomial distribution based on \insertCite{thorson2017model}{SpatialSablefishAssessment}
#' @details \insertCite{thorson2017model}{SpatialSablefishAssessment} deviates from the classic formulation by pulling out (x!) of the denominator from the right hand product
#' and moving it to the left denominator has been validated against extraDistr::ddirmnom
#' @param obs, vector of observed compositions assumes sum(obs) = 1
#' @param est, vector of fiited compositions assumes sum(est) = 1
#' @param beta, the variance inflation coefficient (beta = n * theta)
#' @param n, the total number of samples in the available data (which is restricted to any non-negative real number),
#' @param log_it return log pdf
#' @return pdf
#' @references
#' \insertAllCited{}
#' @export
ddirichmult = function(obs, beta, n, est, log_it = F) {
  if(abs(sum(obs) - 1) > 0.001)
    stop("obs needs to sum to 1")
  if(abs(sum(est) - 1) > 0.001)
    stop("est needs to sum to 1")

  loglike <- lgamma(n + 1) + lgamma(beta) - lgamma(n + beta) +
    sum(lgamma(n * obs + beta * est) - lgamma(beta * est)) - sum(lgamma(n * obs + 1))

  val <- ifelse(log_it == T, loglike, exp(loglike))

  return(val)
}

#' dmultinom_upd
#' an updated version of the multinomial so that we can test TMB's multinomial call
#' @param obs vector of observed values sum = N-effective
#' @param prob vector of fitted compositions assumes sum(prob) = 1
#' @param log return log of of the pdf
#' @return PDF or log of the pdf
#' @export
dmultinom_upd = function(x, prob,log = F) {
  if(abs(sum(prob) - 1) > 0.001)
    stop("prob needs to sum to 1")
  n = sum(x)

  xp1 = x +(1);
  logres = lgamma(n + (1)) - sum(lgamma(xp1)) + sum((x*log(prob)));
  if(log)
    return(logres);

  return(exp(logres));
}

#' log_cv Calculate the CV of the lognormal distribution based on \deqn{cv = \sqrt{e^{\sigma^2} - 1}}
#'
#' @param sigma The standard deviation of the lognormal distribution
#' @return The the cv
#' @export
log_cv = function(sigma) {
  cv = sqrt(exp(sigma^2) - 1)
  return(cv)
}


#' log_sigma
#' @description Calculate the sigma of the lognormal distribution based on \deqn{\sigma = \sqrt{log(cv^2 + 1)}}
#' @param cv The CV (note this is in proportion not percentage) of the lognormal distribution
#' @return The the sigma
#' @export

log_sigma = function(cv) {
  sigma = sqrt(log(cv^2 + 1))
  return(sigma)
}


#' lognormal_CI
#' @description Calculate upper and lower bounds for the lognormal distribution based with expectation and cv
#' @param cv The standard deviation of the lognormal distribution
#' @param expectation The expectation of the distribution, this is not the mu parameter because for the lognormal distribution the expectation is not the mu parameter.
#' @param CI level of confidence (units are proportions not percentage i.e. 0.95 for 95 CI)
#' @export
#' @return a list with upper and lower elements
lognormal_CI <- function(expectation, sigma ,CI = 0.95) {
  mu = (log(expectation) - 0.5*(sigma^2))
  zscore = abs(stats::qnorm((1 - CI)/2))
  U_CI = exp(mu + zscore * sigma)
  L_CI = exp(mu - zscore * sigma)
  return(list("upper" = U_CI, "lower" = L_CI))
}




