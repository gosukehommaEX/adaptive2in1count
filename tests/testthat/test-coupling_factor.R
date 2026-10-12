# Reference values were computed independently in Python by numerical
# integration over the follow-up distribution, without the package.

test_that("complete follow-up leaves the frailty term alone", {
  expect_equal(coupling_factor(2, 1.3, 0.7, 0.8, 0.4, window1 = 0.5),
               0.4 / (0.7 * 0.8), tolerance = 1e-10)
})

test_that("the factor matches an independent evaluation under dropout", {
  expect_equal(coupling_factor(2, 1.3, 0.7, 0.8, 0.4, window1 = 0.5,
                               dropout = 0.4,
                               dropout_shape = "exponential"),
               0.7595851069, tolerance = 1e-7)
  expect_equal(coupling_factor(2, 1.3, 0.7, 0.8, 0.4, window1 = 0.5,
                               dropout = 0.3, dropout_shape = "uniform"),
               0.7421596659, tolerance = 1e-7)
  expect_equal(coupling_factor(10, 10, 1, 1, 1, dropout = 0.4,
                               dropout_shape = "exponential"),
               1.0316831864, tolerance = 1e-7)
})

test_that("no shared frailty leaves the estimators uncorrelated", {
  for (shape in c("uniform", "exponential")) {
    expect_equal(coupling_factor(2, 1.3, 0.7, 0.8, 0, window1 = 1,
                                 dropout = 0.6, dropout_shape = shape), 0)
  }
})

test_that("one endpoint gives the canonical covariance", {
  value <- coupling_factor(2, 1.3, 0.7, 0.8, 0.4, window1 = 0.5,
                           dropout = 0.4, dropout_shape = "exponential",
                           endpoints = "shared")
  expect_equal(value, 2.3848378809, tolerance = 1e-7)
  expect_equal(value, 1 / information_factor(1.3, 0.8, 1, 1, 0.4,
                                             "exponential"),
               tolerance = 1e-8)
  expect_equal(coupling_factor(2, 1.3, 0.7, 0.8, 0.4, window1 = 0.5,
                               dropout = 0.3, dropout_shape = "uniform",
                               endpoints = "shared"),
               2.2515846776, tolerance = 1e-7)
})

test_that("the factor obeys the Cauchy-Schwarz bound", {
  grid <- expand.grid(rate = c(0.5, 2, 10), dispersion = c(0.3, 1, 5),
                      dropout = c(0, 0.3, 0.6),
                      shape = c("uniform", "exponential"),
                      stringsAsFactors = FALSE)
  for (i in seq_len(nrow(grid))) {
    z <- grid[i, ]
    value <- coupling_factor(z$rate, z$rate, z$dispersion, z$dispersion,
                             z$dispersion, window1 = 1, dropout = z$dropout,
                             dropout_shape = z$shape)
    bound <- 1 / information_factor(z$rate, z$dispersion, 1, 1, z$dropout,
                                    z$shape)
    expect_lte(value, bound * (1 + 1e-8))
    expect_gte(value, 0)
  }
})

test_that("dropout raises the factor above the frailty term", {
  # Both weighted exposures increase with the follow-up period, so their
  # product has an expectation above the product of their expectations.
  for (shape in c("uniform", "exponential")) {
    value <- coupling_factor(2, 1.3, 0.7, 0.8, 0.4, window1 = 0.5,
                             dropout = 0.4, dropout_shape = shape)
    expect_gt(value, 0.4 / (0.7 * 0.8))
  }
})

test_that("invalid arguments are rejected", {
  expect_error(coupling_factor(2, 1.3, 0.7, 0.8, 0.9),
               "smaller of the two dispersions")
  expect_error(coupling_factor(2, 1.3, 0.7, 0.8, 0.4, window1 = 2),
               "not exceeding")
})
