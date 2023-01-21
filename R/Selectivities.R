#' logis Logistic selectivity
#' @param x age or length to evaluate selectivity at
#' @param x50 value of x where selectivity is equal to 0.5
#' @param a95 difference in x when the selectivity is equal to 0.5 and 0.95
#' @return selectivity along x
#' @export
logis<- function(x, x50, a95, amax = 1.0) {
  amax/(1+19^((x50-x)/a95))
}

#' logis_alt an alternative logistic selectivity
#' @param x age or length to evaluate selectivity at
#' @param x50 value of x where selectivity is equal to 0.5
#' @param delta delta parameter for the
#' @return selectivity along x
#' @export
logis_alt = function (x, x50, delta, amax = 1.0) {
  return(amax/(1 + exp(-delta * (x - x50))))
}

#' double_normal a double normal selectivity
#' @param x age or length to evaluate selectivity at
#' @param x50 value of x where selectivity is equal to 0.5
#' @param delta delta parameter for the
#' @return selectivity along x
#' @export
double_normal = function (x, x50, delta, amax = 1.0) {
  square <- function(val) {val*val}
  return (pow(x / x50, x50 / (0.5*(sqrt(square(x50)+4*square(delta)) - x50)))*exp((x50 - x)/(0.5*(sqrt(square(x50)+4*square(delta))-x50))));
}

#' constant_sel constant selectivity
#' @param C values
#' @param bins could be age or length
#' @return constant selectivty
#' @export
constant_sel<- function(C, bins) {
  if(C<0| C>1) {
    warning("C must be bound between 0 and 1")
  }
  return(rep(C, length(bins)))
}

#' knife_sel Knife Edge selectivity
#' @param E age at inflection
#' @param bins could be age or length
#' @return selectivty for all ages
#' @export
knife_sel<- function(E, bins) {
  result <- ifelse(bins < E, 0, 1)
  return(result)
}

#' inv_logis_sel Inverse Logisitic selectivity
#' @export
#' @param bins could be age or length
#' @param a50 bin which equates to selectivy being at 0.5
#' @param ato95 bins between the selectivty 0.5 and 0.95
#' @return selectivity for each age.
inv_logis_sel<- function(bins, a50, ato95) {
  1 - 1/(1+19^((a50 - bins)/ato95))
}

#' exp_sel Exponential selectivity
#' @export
#' @param bins could be age or length
#' @param lambda rate decrease
#' @return selectivity for each age.
exp_sel <- function(bins,lambda) {
  exp(-bins*lambda)
}


#' d_norm_sel Double normal selectivity
#' @export
#' @param bins could be age or length
#' @param mu mean value selectivity at one
#' @param sig_l sigma for the left hand curve
#' @param sig_r sigma for the right hand curve
#' @return selectivity for each age
d_norm_sel<- function(bins, mu, sig_l,sig_r) {
  store<- vector()
  for( i in 1:length(bins)) {
    if( bins[i] <= mu) {
      store[i]<- 2^-((bins[i]-mu)/sig_l)^2
    } else {
      store[i]<- 2^-((bins[i]-mu)/sig_r)^2
    }
  }
  return(store)
}


#' d_exp_sel Double Exponential selectivity
#' @export
#' @param bins could be age or length
#' @param x_1 reference bin for the left hand point
#' @param x_2 reference bin for the right hand point
#' @param x_0 reference bin for middle point
#' @param y_1 selectivity at left hand bin
#' @param y_2 selectivity at right hand bin
#' @param y_0 selectivty at middle bin
d_exp_sel <- function(bins, x_1, x_2, x_0, y_0, y_1, y_2)
{
  store<- vector()
  for( i in 1:length(bins)) {
    if( bins[i] <= x_0) {
      store[i]<- min(1, y_0*(y_1/y_0)^((bins[i]-x_0)/(x_1-x_0)))
    } else {
      store[i]<- min(1, y_0*(y_2/y_0)^((bins[i]-x_0)/(x_2-x_0)))
    }
  }
  return(store)
}

#' double_norm_ss3 Double normal selectivity from SS3
#' @param bin_vec could be age or length
#' @param p1 peak parameter beginning size (or age) for the plateau (in cm or year).
#' @param p2 width of plateau, as logistic between peak and maximum length (or age)
#' @param p3 Ascending width: parameter value is ln(width).
#' @param p4 Descending width: parameter value is ln(width).
#' @param p5 y0 logistic transformed selectivity at first bin
#' @param p6 y1 logistic transformed selectivity at last bin
#' @export
double_norm_ss3 <- function(bin_vec, p1, p2, p3, p4, p5, p6) {
  max_x_val <- max(bin_vec)
  # from above, sent by Adam Langley on 2022-02-22
  p1trans <- p1
  p2trans <- p1trans + 1 + (0.99 + max_x_val - p1trans - 1)/(1 + exp(-1.0 * p2))
  p3trans <- exp(p3)
  p4trans <- exp(p4)
  p5trans <- 1/(1 + exp(-1.0 * p5))
  p6trans <- 1/(1 + exp(-1.0 * p6))
  midbin <- bin_vec
  asc <- exp(-((midbin - p1trans)^2/p3trans))
  asc.scaled <- (p5trans + (1 - p5trans) * (asc - 0)/(1 - 0))
  desc <- exp(-((midbin - p2trans)^2/p4trans))
  stj <- exp(-((40 - p2trans)^2/p4trans))
  des.scaled <- (1 + (p6trans - 1) * (desc - 1) /(stj - 1))
  join1 <- 1/(1 + exp(-(20 * (midbin - p1trans)/(1 + abs(midbin - p1trans)))))
  join2 <- 1/(1 + exp(-(20 * (midbin - p2trans)/(1 + abs(midbin - p2trans)))))
  selex <- asc.scaled * (1 - join1) + join1 * (1 * (1 - join2) + des.scaled * join2)
  selex[1] <- exp(p5)
  return(selex)
}
