# Forward stepwise implementation with groups, weights, and Kac-Rice p-values

#' Add one group to the model
#'
#' @param x Design matrix
#' @param y Response vector
#' @param index Group membership indicator of length p
#' @param inactive Indices of inactive groups
#' @param wt Pre-specified weight for each group
#' @param Sigma Error covariance matrix
#' @param mc.cores Number of cores for parallel computing
#' @return Index \code{imax} of added group, residualized data, truncated chi p-value
fstep_add1 <- function(x, y, index, inactive, wt, Sigma, ...) UseMethod("fstep_add1")

fstep_add1.default <- function(x, y, index, inactive, wt, Sigma, mc.cores = 1, ...) {

  p <- ncol(x)
  Xy <- t(x) %*% y

  if (missing(index)) index <- 1:p
  labels <- unique(index)
  if (missing(wt)) wt <- c(rep(1, length(labels)))
  if (missing(inactive)) inactive <- labels
  if (missing(Sigma)) {
    Sigma <- 1
    warning("Error covariance unspecified, using diag(1)")
  }

  terms <- lapply(inactive, function(i) {
    inds <- which(index == i)
    return(sum((Xy[inds])^2)/wt[i]^2)
  })

  terms.maxind <- which.max(terms)
  imax <- inactive[terms.maxind]
  maxinds <- which(index == imax)
  L <- sqrt(terms[[terms.maxind]])
  w_max <- wt[imax]
  Xy_max <- Xy[maxinds]
  Xmax <- x[, maxinds]

  active.inds <- which(!index %in% inactive)
  inactive.inds <- which(index %in% setdiff(inactive, imax))
  if (length(inactive) < length(labels)) {
    Xmax <- lm(Xmax ~ x[, active.inds])$residuals
  }

  Xmax_ncol <- length(maxinds)
  if (is.matrix(Xmax)) {
    k <- rankMatrix(Xmax)[1]
  } else {
    k <- 1
  }

  eta <- rep(0, p)
  eta[maxinds] <- Xy_max / (sqrt(sum(Xy_max^2)) * w_max)

  if (Xmax_ncol > 1) {
    # Group of size > 1, non-trivial gradient term
    tangent.space <- Null(Xy_max)
    V <- matrix(0, ncol = Xmax_ncol, nrow = p)
    V[maxinds, ] <- tangent.space
    XV <- x %*% V[, -ncol(V)]
    XV <- cbind(XV, rep(0, nrow(XV)))
    XXy <- Xmax %*% Xy_max

    if (is.matrix(Sigma)) {
      SXV <- Sigma %*% XV
    } else {
      SXV <- Sigma*XV
    }

    P <- SXV %*% ginv(t(XV) %*% SXV) %*% t(XV)
    OP <- diag(rep(1, ncol(P))) - P

    if (is.matrix(Sigma)) {
      H <- OP %*% Sigma
    } else {
      H <- OP*Sigma
    }

    temp_den <- w_max * sqrt(sum(Xy_max^2))
    Xeta <- H %*% XXy
    conditional_var <- t(XXy) %*% Xeta
    Xeta <- Xeta / temp_den
    conditional_var <- diag(conditional_var) / temp_den^2

  } else {
    # Group of size 1, no gradient term
    Xeta <- Xmax / w_max * sign(Xy_max)
    conditional_var <- sum(Xeta^2)
  }

  if (conditional_var <= 0) {
    p.value <- 1
    lower_bound <- NaN
    upper_bound <- NaN
  } else {
    Xeta <- Xeta / conditional_var
    b <- t(x) %*% Xeta
    a <- Xy - b * L

    inactive_new <- setdiff(inactive, imax)
    if (mc.cores > 1) {
      trig_solutions <- mclapply(inactive_new, function(i) {
        inds <- which(index == i)
        w <- wt[i]
        return(tchi_trig_solution(num=a[inds], den=b[inds], w=w))
      }, mc.cores = mc.cores)

    } else {
      trig_solutions <- lapply(inactive_new, function(i) {
        inds <- which(index == i)
        w <- wt[i]
        return(tchi_trig_solution(num=a[inds], den=b[inds], w=w))
      })
    }

    trig_solutions <- do.call(cbind, trig_solutions)
    V_lower <- trig_solutions[1, ]
    V_upper <- trig_solutions[2, ]

    if (length(V_lower) >= 1) {
      V_lower <- V_lower[!is.nan(V_lower)]
      lower_bound <- max(V_lower)
    } else {
      lower_bound <- 0
    }
    if (length(V_upper) >= 1) {
      V_upper <- V_upper[!is.nan(V_upper)]
      upper_bound <- min(V_upper)
    } else {
      upper_bound <- Inf
    }

    p.value <- tchi_cdf(L, lower_bound, upper_bound, conditional_var, k)
  }

  ## if ((p.value < 0) | (p.value > 1)) {
  ##   print(V_lower)
  ##   print(lower_bound)
  ##   print(L)
  ##   print(upper_bound)
  ##   print(V_upper)
  ##   print("-----------------------------------------------")
  ## }

  return(list(imax=imax, p.value=p.value, k=k, conditional_var=conditional_var, L=L, lower_bound=lower_bound, upper_bound=upper_bound))
}


#' Forward stepwise, computing tchi p-values along the path
#'
#' @param x Design matrix
#' @param y Response vector
#' @param index Group membership indicator of length p
#' @param wt Pre-specified weight for each group
#' @param Sigma Error covariance matrix
#' @param steps Maximum number of steps for forward stepwise
#' @param normalize Should the design matrix be normalized?
#' @param mc.cores Number of cores for parallel computing
#' @return Index \code{imax} of added group, residualized data, truncated chi p-value
fstep <- function(x, y, index, wt, Sigma, steps, normalize = TRUE, mc.cores = 1, ...) UseMethod("fstep")

fstep.default <- function(x, y, index, wt, Sigma, steps, normalize = TRUE, mc.cores = 1, verbose=FALSE, ...) {

  p <- ncol(x)
  n <- nrow(x)

  if (missing(index)) index <- 1:p
  labels <- unique(index)
  if (missing(wt)) wt <- c(rep(1, length(labels)))
  if (missing(Sigma)) {
    Sigma <- 1
    warning("Error covariance unspecified, using diag(1)")
  }
  if (missing(steps)) steps <- min(n, p) - 1
  inactive <- labels

  pvals <- c()
  active <- c()
  x.update <- x
  if (normalize) {
    for (i in labels) {
      inds <- which(index == i)
      frob.norm <- sqrt(sum(x[, inds]^2))
      if (frob.norm > 0) {
        x.update[, inds] <- x.update[, inds] / frob.norm
      }
    }
  }

  y.update <- y
  y.last <- y

  output <- data.frame(imax=integer(), p.value=numeric(), k=integer(), conditional_var=numeric(), L=numeric(), lower_bound=numeric(), upper_bound=numeric())
  for (step in 1:steps) {
    added <- fstep_add1(x.update, y.update, index, inactive, wt, Sigma, ...)

    imax <- added$imax

    inactive <- setdiff(inactive, imax)

    active.inds <- which(index %in% active)
    inactive.inds <- which(!index %in% active)
    imax.inds <- which(index == imax)
    if (length(active) > 0) {
      x.imax <- lm(x.update[, imax.inds] ~ x.update[, active.inds])$residuals
    } else {
      x.imax <- x.update[, imax.inds]
    }
    active <- union(active, imax)
    x.update[, imax.inds] <- x.imax
    #P <- diag(rep(1, n)) - x.imax %*% ginv(x.imax)

    #y.update <- P %*% y.update
    y.update <- lm(y.update ~ x.imax)$residuals
    #x.update[, inactive.inds] <- P %*% x.update[, inactive.inds]

    added$RSSdrop <- sum((y.last - y.update)^2)
    y.last <- y.update

    output <- rbind(output, data.frame(added))
    if (verbose) print(added)
  }

  value <- list(variable=output$imax, p.value=output$p.value, log=output)
  class(value) <- "fstep"
  invisible(value)
}

