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
#' @param Y scalar range -inf, inf
#' @export
#' @return X to be between -1,1
bound_unit = function(Y) {
  return(Y / sqrt(1.0 + Y * Y))
}

#' inv_bound_unit constrains Y from -inf -> inf to be between -1 -> 1
#' @param X scalar range -1,1
#' @export
#' @return Y to be between -inf, inf
inv_bound_unit = function(X) {
  return(sqrt((X*X) / (1 - X*X)) * ifelse(X < 0, -1, 1))
}
#' logit bounds X which is between 0-1 to -inf -> inf based on the logit transformation
#' equivalent to qlogis(X)
#' @param X scalar range [0,1]
#' @export
#' @return Y to be between -inf, inf
logit = function(X) {
  log(X / (1 - X))
}
#' invlogit Inverse logit transformation, equivalent to plogis(Y)
#' @param Y scalar between -inf, inf
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

#' Q_sum_to_zero_QR
#' calculate QR vector for sum to zero hard constraint
#' @param N integer number of elements in vector that we want to sum = 0
#' @return a vector length N * 2 that will be used by sum_to_zero_QR
#' @export
Q_sum_to_zero_QR <- function(N) {
  Q_r = vector(length = N * 2);

  for(i in 1:(N - 1)) {
    Q_r[i] = -sqrt((N-i)/(N-i+1.0));
    Q_r[i+N] = 1.0 / sqrt((N-i) * (N-i+1));
  }
  return (Q_r);
}

#' sum_to_zero_QR
#' take a vector of unconstrained values length (N - 1) and derive a vector of length N that sum = 0 using the QR method
#' see here https://discourse.mc-stan.org/t/test-soft-vs-hard-sum-to-zero-constrain-choosing-the-right-prior-for-soft-constrain/3884
#' @param x_raw vector of unconstrained values length N - 1
#' @return a vector length N that sums = 0
#' @export
sum_to_zero_QR <- function(x_raw) {
  N = length(x_raw) + 1;
  Q_r = Q_sum_to_zero_QR(N);
  x = vector(length = N) ;
  x_aux = 0;

  for(i in 1:(N-1)){
    x[i] = x_aux + x_raw[i] * Q_r[i];
    x_aux = x_aux + x_raw[i] * Q_r[i+N];
  }
  x[N] = x_aux;
  return(x);
}
#'
#' pow
#' utility function to mimic TMB code and likelihood evaluations
#' @param val base value
#' @param exponent exponent to take to the power
#' @return val to the exponent
#' @export
pow <- function(val,exponent) {val^exponent}

#'
#' posfun
#' utility function to mimic TMB code and likelihood evaluations
#' @param x a value to check if its greater than epx
#' @param eps a small value to check if x is greater than
#' @param pen a value which will get incremented by the penalty of the posfun
#' @return a modified value of x which is larger than eps in a differential manor
#' @export
posfun <- function(x, eps, pen = NULL) {
  xp = -(x/eps-1);
  if(x >= eps) {
    return(x)
  }
  return(eps*(1/(1+xp+pow(xp,2)+pow(xp,3)+pow(xp,4)+pow(xp,5))))
}


#' extend_vec_last_val
#' @param vector a vector to add n elements to, all of which have the same value as the last value in vector
#' @param n number of elements to concatenate
#' @return a vector with length(vector) + n all the last n values have the same value as length(vector)
#' @export
extend_vec_last_val <- function(vector, n) {
  new_vec = c(vector, rep(vector[length(vector)], n))
  return(new_vec)
}

#' extend_2darray
#' @param array_2d an array with 2 dimensions or matrix
#' @param n number of times to concatenate the last element in the 2 dimension
#' @param colwise boolean whether to append columns (TRUE) or rows (FALSE)
#' @return a two dimensional array that has had the either the first or second dimension concatenated
#' @export
extend_2darray <- function(array_2d, n, colwise= TRUE) {
  if(colwise) {
    return(cbind(array_2d, replicate(array_2d[,ncol(array_2d)], n = n)))
  } else {
    return(rbind(array_2d, replicate(array_2d[nrow(array_2d),], n = n)))
  }
  return(NULL)
}
#' extend_3darray_last_dim
#' @param array_3d an array with 3 dimensions
#' @param n number of times to concatenate the last element in the 3 dimension
#' @return a three dimensional array that has had the 3dimension concatenated
#' @export
extend_3darray_last_dim <- function(array_3d, n) {
  new_3d_array = abind(array_3d, replicate(array_3d[,,dim(array_3d)[3]], n = n), along = 3)
  return(new_3d_array)
}

#' convert_simdata_integers
#' @param sim_data an array with 3 dimensions
#' @param OM_data data that was passed to `MakeADFun`
#' @details this function should be used if you get the following error from `MakeADFun` `Error in getParameterOrder(data, parameters, new.env(), DLL = DLL) : NOT A VECTOR!`
#' @return sim_data list with integers converted
#' @export
convert_simdata_integers <-function(sim_data, OM_data) {
  if(sim_data$model == "Assessment") {
    sim_data$n_regions = OM_data$n_regions
    sim_data$n_projections_years = OM_data$n_projections_years
    sim_data$do_projection = OM_data$do_projection
    sim_data$global_rec_devs = OM_data$global_rec_devs
    sim_data$ll_catchatage_covar_structure = OM_data$ll_catchatage_covar_structure
    sim_data$ll_catchatage_comp_likelihood = OM_data$ll_catchatage_covar_structure
    sim_data$ll_catchatlgth_covar_structure = OM_data$ll_catchatlgth_covar_structure
    sim_data$ll_catchatlgth_comp_likelihood = OM_data$ll_catchatlgth_comp_likelihood
    sim_data$trwl_catchatlgth_covar_structure = OM_data$trwl_catchatlgth_covar_structure
    sim_data$trwl_catchatlgth_comp_likelihood = OM_data$trwl_catchatlgth_comp_likelihood
    sim_data$dom_ll_bio_likelihood = OM_data$dom_ll_bio_likelihood
    sim_data$jap_ll_bio_likelihood = OM_data$jap_ll_bio_likelihood
    sim_data$ll_cpue_likelihood = OM_data$ll_cpue_likelihood
    sim_data$srv_dom_ll_age_covar_structure = OM_data$srv_dom_ll_age_covar_structure
    sim_data$srv_dom_ll_age_comp_likelihood = OM_data$srv_dom_ll_age_comp_likelihood
    sim_data$srv_dom_ll_lgth_covar_structure = OM_data$srv_dom_ll_lgth_covar_structure
    sim_data$srv_dom_ll_lgth_comp_likelihood = OM_data$srv_dom_ll_lgth_comp_likelihood
    sim_data$srv_jap_ll_age_covar_structure = OM_data$srv_jap_ll_age_covar_structure
    sim_data$srv_jap_ll_age_comp_likelihood = OM_data$srv_jap_ll_age_comp_likelihood
    sim_data$srv_jap_ll_lgth_covar_structure = OM_data$srv_jap_ll_lgth_covar_structure
    sim_data$srv_jap_ll_lgth_comp_likelihood = OM_data$srv_jap_ll_lgth_comp_likelihood
    sim_data$srv_nmfs_trwl_age_covar_structure = OM_data$srv_nmfs_trwl_age_covar_structure
    sim_data$srv_nmfs_trwl_age_comp_likelihood = OM_data$srv_nmfs_trwl_age_comp_likelihood
    sim_data$srv_nmfs_trwl_lgth_covar_structure = OM_data$srv_nmfs_trwl_lgth_covar_structure
    sim_data$srv_nmfs_trwl_lgth_comp_likelihood = OM_data$srv_nmfs_trwl_lgth_comp_likelihood
    sim_data$nmfs_trwl_bio_likelihood = OM_data$nmfs_trwl_bio_likelihood
    sim_data$jap_fishery_ll_bio_likelihood = OM_data$jap_fishery_ll_bio_likelihood
    sim_data$srv_jap_fishery_ll_lgth_covar_structure = OM_data$srv_jap_fishery_ll_lgth_covar_structure
    sim_data$srv_jap_fishery_ll_lgth_comp_likelihood = OM_data$srv_jap_fishery_ll_lgth_comp_likelihood
    sim_data$catch_likelihood = OM_data$catch_likelihood
  }
  return(sim_data)
}
