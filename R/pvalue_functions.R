# Functions for computing truncated chi p-value

#' Trigonometry for solving 1-dimensional problem
#'
#' @param num Numerator term
#' @param den Denominator term
#' @param w Weight of chosen group
#' @param tol Numerical tolerance
#' @return Upper and lower limits for truncated chi
#' @references \url{http://arxiv.org/abs/1405.3920}
#' @seealso \code{\link{tchi_cdf}}
tchi_trig_solution = function(num, den, w, tol=1.e-10) {

  norm_num = sqrt(sum(num^2))
  norm_den = sqrt(sum(den^2))

  # Next two cases: effectively no truncation
  if (norm_den == 0) {
    return(c(0, Inf))
  }

  if ((norm_num / norm_den) < tol) {
    return(c(0, Inf))
  }

  Ctheta = sum(num * den) / (norm_num * norm_den)
  # R might believe 1 > 1 and return NaN in sqrt() below
  # so we restrict to (-1, 1) first
  Ctheta = min(max(Ctheta, -1), 1)
  Stheta = sqrt(1-Ctheta^2)
  theta = acos(Ctheta)
  Sphi = norm_den * Stheta / w

  # Again, no truncation:
  if (Sphi > 1) {
    return(c(0, Inf))
  }

  phi1 = asin(Sphi)
  phi2 = pi - phi1

  V1 = norm_num * cos(phi1) / (w - norm_den * cos(theta-phi1))
  V2 = norm_num * cos(phi2) / (w - norm_den * cos(theta-phi2))

  if (norm_den < w) {
    # No upper truncation
    return(c(max(V1,V2), Inf))
  }

  return(c(min(V1,V2), max(V1,V2)))
}


#' Truncated CDF transform
#'
#' @param L Lambda, the maximum of \deqn{\| X_g^T y\|_2}
#' @param lower_bound Truncation lower bound
#' @param upper_bound Truncation upper bound
#' @param cvar Conditional \deqn{\sigma}
#' @param k Rank of maximizing group
#' @return Truncated chi p-value
#' @references \url{http://arxiv.org/abs/1405.3920}
#' @seealso \code{\link{tchi_trig_solution}}
tchi_cdf <- function(L, lower_bound, upper_bound, cvar, k) {

  L2 = L^2
  ub2 = upper_bound^2
  lb2 = lower_bound^2
  first.term = pchisq(ub2/cvar, k, lower.tail=TRUE)

  if (first.term == 1) {
    num = pchisq(L2/cvar, k, lower.tail=FALSE, log.p=TRUE)
    den = pchisq(lb2/cvar, k, lower.tail=FALSE, log.p=TRUE)
    value = exp(num - den)

  } else {
    # Is this numerically unstable?
    num = first.term - pchisq(L2/cvar, k, lower.tail=TRUE)
    den = first.term - pchisq(lb2/cvar, k, lower.tail=TRUE)
    value = num/den
  }

  return(value)
}
