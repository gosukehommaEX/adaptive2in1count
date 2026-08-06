test_that("independence factorises the probability", {
  for (a in c(-1.5, 0, 2)) {
    for (b in c(-1, 0.5, 1.8)) {
      expect_equal(bivariate_normal_lower(a, b, 0), pnorm(a) * pnorm(b),
                   tolerance = 1e-9)
    }
  }
})

test_that("perfect correlation gives the smaller marginal", {
  expect_equal(bivariate_normal_lower(0.4, 1.2, 1), pnorm(0.4),
               tolerance = 1e-8)
  expect_equal(bivariate_normal_lower(-0.5, -2, 1), pnorm(-2),
               tolerance = 1e-8)
})

test_that("perfect negative correlation gives the Frechet lower bound", {
  expect_equal(bivariate_normal_lower(1, 1, -1), 2 * pnorm(1) - 1,
               tolerance = 1e-8)
  expect_equal(bivariate_normal_lower(-1, -1, -1), 0, tolerance = 1e-8)
})

test_that("the probability increases with the correlation at the median", {
  values <- vapply(c(-0.8, -0.4, 0, 0.4, 0.8), function(r) {
    bivariate_normal_lower(0, 0, r)
  }, numeric(1))
  expect_true(all(diff(values) > 0))
  expect_equal(values[3], 0.25, tolerance = 1e-9)
})

test_that("the known value at the median matches the arcsine formula", {
  for (rho in c(-0.6, -0.2, 0.3, 0.75)) {
    expect_equal(bivariate_normal_lower(0, 0, rho),
                 0.25 + asin(rho) / (2 * pi), tolerance = 1e-8)
  }
})

test_that("infinite limits reduce to a marginal", {
  expect_equal(bivariate_normal_lower(Inf, 0.3, 0.5), pnorm(0.3),
               tolerance = 1e-9)
  expect_equal(bivariate_normal_lower(0.3, Inf, 0.5), pnorm(0.3),
               tolerance = 1e-9)
  expect_equal(bivariate_normal_lower(-Inf, 0.3, 0.5), 0)
})

test_that("the probability never leaves the unit interval", {
  grid <- expand.grid(a = c(-3, -1, 0, 1, 3), b = c(-3, -1, 0, 1, 3),
                      rho = c(-0.9, -0.3, 0.3, 0.9))
  for (i in seq_len(nrow(grid))) {
    value <- bivariate_normal_lower(grid$a[i], grid$b[i], grid$rho[i])
    expect_gte(value, 0)
    expect_lte(value, 1)
  }
})

test_that("invalid arguments are rejected", {
  expect_error(bivariate_normal_lower(0, 0, 1.5), "between -1 and 1")
})
