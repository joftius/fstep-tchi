context("Forward stepwise")

test_that("fstep_add1 adds the correct group", {
  x <- matrix(rep(1, 25), nrow=5)
  x[lower.tri(x)] <- 0
  y <- as.numeric(1:5 == 5)
  index = 1:5

  expect_that(fstep_add1(x, y, index)$imax, equals(5))
  y[2] = -1
  index = c(1, 2, 3, 3, 4)
  expect_that(fstep_add1(x, y, index)$imax, equals(3))
  y[5] = 3
  expect_that(fstep_add1(x, y, index)$imax, equals(4))
  }
)
