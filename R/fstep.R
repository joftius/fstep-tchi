# Forward stepwise implementation with groups, weights, and Kac-Rice p-values

#' Add one group to the model
#'
#' @param x Design matrix
#' @param y Response vector
#' @param index Group membership indicator of length p
#' @return Index \code{imax} of added group, residualized data, truncated chi p-value
fstep_add1 <- function(x, y, index, ...) UseMethod("fstep_add1")

fstep_add1.default <- function(x, y, index, ...) {
  U = t(x)
  labels = unique(index)
  terms = lapply(labels, function(i) {
    inds = which(index == i)
    return(sum((U[inds, ] %*% y)^2))
  })
  imax = which.max(terms)
  return(list(imax = imax))
}
