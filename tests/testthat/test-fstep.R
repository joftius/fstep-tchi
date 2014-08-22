context("Forward stepwise")


test_that("fstep adds the correct group", {
  x <- matrix(rep(1, 25), nrow=5)
  x[lower.tri(x)] <- 0
  y <- as.numeric(1:5 >= 3) #+ .1*rnorm(5)
  index = 1:5

  expect_that(fstep_add1(x, y, index)$imax, equals(5))
  expect_that(fstep_add1(x, y, index, inactive=1:3)$imax, equals(3)) # with inactive

  index = c(1, 2, 3, 3, 4)
  expect_that(fstep_add1(x, y, index)$imax, equals(4)) # with groups

  weights = c(1, 1, 1, 3)
  expect_that(fstep_add1(x, y, index, weights)$imax, equals(3)) # and with weights

  x <- toeplitz(.8^c(1:50))[1:25, ]
  y <- rep(1, 25)
  expect_that(fstep_add1(x, y)$imax, equals(13))
  }
)

test_that("p-values lie in the unit interval", {
  x <- matrix(rep(1,100), nrow=10)
  x[lower.tri(x)] <- 0
  y <- as.numeric(1:10 >= 5)
  fit = fstep(x, y, steps = 8)
  expect_that(min(fit$p.value) >= 0, is_true())
  expect_that(max(fit$p.value) <= 1, is_true())

  #y <- rnorm(10)
  y <- c(rep(-1,5), rep(1,5))
  fit = fstep(x, y, steps = 8)
  expect_that(min(fit$p.value) >= 0, is_true())
  expect_that(max(fit$p.value) <= 1, is_true())

  x <- toeplitz(.8^c(1:50))[1:25,]
  y <- as.numeric(1:25 >= 10)
  fit = fstep(x, y, steps = 20)
  expect_that(min(fit$p.value) >= 0, is_true())
  expect_that(max(fit$p.value) <= 1, is_true())

  y <- 10 * rnorm(25)
  fit = fstep(x, y, steps = 20)
  expect_that(min(fit$p.value) >= 0, is_true())
  expect_that(max(fit$p.value) <= 1, is_true())

  x <- matrix(rnorm(100*5), nrow=100)
  y <- 2*x[, 1] - x[, 2] + rnorm(100)
  fit = fstep(x, y, steps = 4)
  expect_that(min(fit$p.value) >= 0, is_true())
  expect_that(max(fit$p.value) <= 1, is_true())

  y <- sample(c(-1,0,1), size=100, replace=TRUE)
  fit = fstep(x, y, steps = 4)
  expect_that(min(fit$p.value) >= 0, is_true())
  expect_that(max(fit$p.value) <= 1, is_true())
  }
)
