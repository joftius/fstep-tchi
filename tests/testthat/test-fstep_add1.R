context("Forward stepwise")


test_that("fstep_add1 adds the correct group", {
  x <- matrix(rep(1, 25), nrow=5)
  x[lower.tri(x)] <- 0
  y <- as.numeric(1:5 >= 3)
  index = 1:5

  expect_that(fstep_add1(x, y, index)$imax, equals(5))
  expect_that(fstep_add1(x, y, index, inactive=1:3)$imax, equals(3)) # with inactive
  index = c(1, 2, 3, 3, 4)
  expect_that(fstep_add1(x, y, index)$imax, equals(4)) # with groups
  weights = c(1, 1, 1, 3)
  expect_that(fstep_add1(x, y, index, weights)$imax, equals(3)) # and with weights
  }
)
