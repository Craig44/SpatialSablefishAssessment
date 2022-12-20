#' Utility function for executing a string in R
#'
#' @author Craig Marsh
#' @param x string
#' @return will exectute the string
#' @export
#'
evalit <- function(x) {
  result <- eval(parse(text = x))
  return(result)
}
#' zerofun
#'
#' @param x value to check if is zero
#' @param delta how small the value is before we alter it
#' @return adjusted x if x < delta
#' @export
zerofun = function(x, delta) {
  if (x >= delta)
    return(x);

  return (delta / (2.0 - (x / delta)));
}
#' zerofun_v vectorised version of zerofun
#'
#' @param x value to check if is zero
#' @param delta how small the value is before we alter it
#' @return adjusted x if x < delta
#' @export

zerofun_v = Vectorize(zerofun)

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

#' gm_mean calculate geometric mean when you have na's in vector
#' @param x vector of strictly positive values
#' @param na.rm boolean on whether to ignore NA's
#' @return geometric mean with na's ignored
#' @export
#'
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

#' bound_unit constrains Y from -inf -> inf to be between -1 -> 1
#' @param Y scalar range [-inf, inf]
#' @export
#' @return X to be between [-1,1]
bound_unit = function(Y) {
  return(Y / sqrt(1.0 + Y * Y))
}

#' inv_bound_unit constrains Y from -inf -> inf to be between -1 -> 1
#' @param X scalar range [-1,1]
#' @export
#' @return Y to be between [-inf, inf]
inv_bound_unit = function(X) {
  return(sqrt((X*X) / (1 - X*X)) * ifelse(X < 0, -1, 1))
}
#' logit bounds X which is between 0-1 to -inf -> inf based on the logit transformation
#' equivalent to qlogis(X)
#' @param X scalar range [0,1]
#' @export
#' @return Y to be between [-inf, inf]
logit = function(X) {
  log(X / (1 - X))
}
#' invlogit Inverse logit transformation, equivalent to plogis(Y)
#' @param Y scalar between [-inf, inf]
#' @export
#' @return X between [0,1]
invlogit<- function(Y) {
  1/(1 + exp(-Y))
}
#' simplex
#' takes a unit vector (or a vector that will be scaled by the mean) of length K and converts to unconstrained K - 1 vector
#' @param xk vector of length K - 1 takes unconstrained values
#' @param sum_to_one whether to rescale xk so it sums to one
#' @return vector of length K - 1 unconstrained values
#' @export
#'
simplex <- function (xk, sum_to_one = TRUE)  {
  zk = vector()
  if (!sum_to_one) {
    xk = xk/sum(xk)
  }
  else {
    if (abs(sum(xk) - 1) > 0.001)
      stop("xk needs to sum = 1, otherwise speify sum_to_one = TRUE")
  }
  K = length(xk)
  zk[1] = xk[1]/(1)
  for (k in 2:(K - 1)) {
    zk[k] = xk[k]/(1 - sum(xk[1:(k - 1)]))
  }
  yk = stats::qlogis(zk) - log(1/(K - 1:(K - 1)))
  return(yk)
}
#' restoresimplex
#' takes an unconstrained (each element is unbounded betweeen -Inf and Inf) vector of length K - 1, and calculates a unit vector of length K
#' @param yk vector of length K - 1 takes unconstrained values
#' @return xk unit vector of length K
#' @export
#'
restoresimplex <- function (yk) {
  K = length(yk) + 1
  zk = stats::plogis(yk + log(1/(K - 1:(K - 1))))
  xk = vector()
  xk[1] = zk[1]
  for (k in 2:(K - 1)) {
    xk[k] = (1 - sum(xk[1:(k - 1)])) * zk[k]
  }
  xk[K] = 1 - sum(xk)
  return(xk)
}

#' vonbert applies the Von Bertalanffy age-length relationship
#' @param age take an age look in globe scope for the rest of parameters
#' @param L_inf asympototic length
#' @param K growth rate parameter
#' @param t0 age for length = 0
#' @export
#' @return mean length at age
vonbert <- function(age,K,L_inf,t0) {
  return(L_inf * (1-exp(-K*(age -t0))))
}
#' logit_general bounds X which is between [lb,ub] to -inf -> inf based on the logit transformation
#' @param X scalar range [lb,ub]
#' @param ub upper bound for X
#' @param lb lower bound for X
#' @export
#' @return Y to be between [-inf, inf]
logit_general = function(X, lb, ub) {
  X1 = (X - lb) / (ub - lb)
  log(X1/(1 - X1))
}
#' invlogit_general bounds X which is between -inf -> inf to [lb,ub] based on the logit transformation
#' @param Y scalar range [-inf, inf]
#' @param ub upper bound for X
#' @param lb lower bound for X
#' @export
#' @return X to be between [lb,ub]
invlogit_general = function(Y, lb, ub) {
  Y1 = 1 / (1 + exp(-Y))
  lb + (ub - lb)*Y1
}
#' is_matrix_invertable
#' @description helper function to see if a matrix is invertable
#' @param m an n x n matrix
#' @export
#' @return bool
#'  \itemize{
#'   \item false: not invertable
#'   \item true: is invertable
#' }
#' @examples
#'\dontrun{
#' x <- matrix(rep(1,25),nc=5)          # singular
#' y <- matrix(1+1e-10*rnorm(25),nc=5)  # very nearly singular matrix
#' z <- 0.001*diag(1,5)                 # non-singular, but very smalll determinant
#' is_matrix_invertable(x)
#' # [1] FALSE
#' is_matrix_invertable(y)
#' # [1] TRUE
#' is_matrix_invertable(z)
#' # [1] TRUE
#' }
is_matrix_invertable <- function(m) {
  any("matrix" %in% class(try(solve(m),silent=TRUE)))
}

#' every_nth
#' @description return a vector of boolean that has every nth value = true
#' @param n every 'n' values will equal TRUE
#' @export
#' @return bool

every_nth = function(n) {
  return(function(x) {x[c(TRUE, rep(FALSE, n - 1))]})
}
## A function to check if the covariance matrix is positive definite
#' is_positive_definite
#' @description helper function to see if a matrix is positive definite
#' @param m an n x n matrix
#' @param tol real value. Eigen values less than this tolerance value will fail the check.
#' @export
#' @return bool
#'  \itemize{
#'   \item false: not a positive definite matrix
#'   \item true: is a positive definite matrix
#' }
is_positive_definite <- function (m, tol = 1e-6)  {
  ## check it is symetric
  if(!isSymmetric(m, tol = tol)) {
    message("Matrix 'm' was not symmetric. Derived using 'isSymmetric(m, tol = tol)'")
    return(FALSE)
  }
  ## get eigen values
  eS <- eigen(m, symmetric = TRUE)
  ev <- eS$values
  ## check with tolerance
  n <- nrow(m)
  for (i in 1:n) {
    if (abs(ev[i]) < tol) {
      ev[i] <- 0
    }
  }
  if (any(ev <= 0)) {
    return(FALSE)
  }
  return(TRUE)
}

#' is_constant
#' is a vector constant
#' @param x vector
#' @param tol the magnitude that flags a difference
is_constant <- function (x, tol = .Machine$double.eps) {
  abs(max(x) - min(x)) < tol
}

# should NA be returned by a convergence diagnostic?
should_return_NA <- function(x) {
  anyNA(x) || any(!is.finite(x)) || is_constant(x)
}

#' unpaste
#' @export
#' @return unpasted string
unpaste <- function (string, sep)  {
  return(unlist(strsplit(string, split = sep)))
}
