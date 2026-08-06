test_that("the convolution of two negative binomials is negative binomial", {
  settings <- list(c(10, 20, 15, 30), c(4, 2, 6, 3), c(30, 5, 60, 10))
  for (s in settings) {
    ymax <- nb_support_max(s[1] + s[3], s[2] + s[4], tol = 1e-13)
    f <- nb_count_pmf(s[1], s[2], ymax)
    g <- nb_count_pmf(s[3], s[4], ymax)
    expect_equal(convolve_pmf(f, g, ymax), nb_count_pmf(s[1] + s[3],
                                                        s[2] + s[4], ymax),
                 tolerance = 1e-8)
  }
})

test_that("the mean of the convolution is the sum of the means", {
  ymax <- 400
  f <- nb_count_pmf(4, 6, ymax)
  g <- nb_count_pmf(9, 11, ymax)
  expect_equal(sum(convolve_pmf(f, g, ymax) * seq(0, ymax)), 17,
               tolerance = 1e-7)
})

test_that("convolving with a point mass at zero returns the input", {
  f <- nb_count_pmf(3, 5, 80)
  point <- c(1, rep(0, 80))
  expect_equal(convolve_pmf(f, point, 80), f, tolerance = 1e-12)
})

test_that("the convolution is symmetric in its two arguments", {
  f <- nb_count_pmf(3, 5, 120)
  g <- nb_count_pmf(7, 9, 120)
  expect_equal(convolve_pmf(f, g, 120), convolve_pmf(g, f, 120),
               tolerance = 1e-12)
})

test_that("the convolution is associative", {
  f <- nb_count_pmf(2, 3, 150)
  g <- nb_count_pmf(3, 4, 150)
  h <- nb_count_pmf(4, 5, 150)
  left <- convolve_pmf(convolve_pmf(f, g, 150), h, 150)
  right <- convolve_pmf(f, convolve_pmf(g, h, 150), 150)
  expect_equal(left, right, tolerance = 1e-10)
})

test_that("truncation never adds probability mass", {
  f <- nb_count_pmf(3, 20, 200)
  g <- nb_count_pmf(3, 20, 200)
  expect_lte(sum(convolve_pmf(f, g, 30)), 1)
  expect_length(convolve_pmf(f, g, 30), 31)
})

test_that("invalid arguments are rejected", {
  expect_error(convolve_pmf(numeric(0), 1, 10), "non-empty")
  expect_error(convolve_pmf(1, numeric(0), 10), "non-empty")
  expect_error(convolve_pmf(1, 1, -2), "non-negative")
})
