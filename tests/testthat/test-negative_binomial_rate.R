test_that("equal exposure reproduces the moment estimator", {
  set.seed(1)
  count <- matrix(rpois(200, 3), nrow = 20)
  exposure <- matrix(1, nrow = 20, ncol = 10)
  out <- negative_binomial_rate(count, exposure, dispersion = 0.8)
  expect_equal(out$rate, rowSums(count) / 10, tolerance = 1e-10)
})

test_that("the estimate does not depend on the dispersion at equal exposure", {
  set.seed(2)
  count <- matrix(rpois(150, 2), nrow = 15)
  exposure <- matrix(0.5, nrow = 15, ncol = 10)
  first <- negative_binomial_rate(count, exposure, dispersion = 0.3)$rate
  second <- negative_binomial_rate(count, exposure, dispersion = 5)$rate
  expect_equal(first, second, tolerance = 1e-9)
})

test_that("the score equation is solved at the returned rate", {
  set.seed(3)
  exposure <- matrix(runif(120, 0.2, 1), nrow = 12)
  count <- matrix(rpois(120, 2 * exposure), nrow = 12)
  dispersion <- 0.7
  out <- negative_binomial_rate(count, exposure, dispersion)
  score <- rowSums(count / out$rate -
                     (count + dispersion) * exposure /
                     (dispersion + out$rate * exposure))
  expect_true(all(abs(score) < 1e-6))
})

test_that("the information matches the sum over patients", {
  exposure <- matrix(c(1, 0.5, 1, 0.25), nrow = 1)
  count <- matrix(c(2, 1, 3, 0), nrow = 1)
  dispersion <- 0.8
  out <- negative_binomial_rate(count, exposure, dispersion)
  mean_count <- out$rate * exposure
  expect_equal(out$information,
               sum(dispersion * mean_count / (dispersion + mean_count)),
               tolerance = 1e-10)
})

test_that("the information matches the design quantity under complete follow-up", {
  size <- 400
  rate <- 1.3
  dispersion <- 0.8
  exposure <- matrix(1, nrow = 1, ncol = size)
  count <- matrix(round(rate * exposure), nrow = 1)
  out <- negative_binomial_rate(count, exposure, dispersion)
  expect_equal(out$information / size,
               information_factor(out$rate, dispersion, 1, 1),
               tolerance = 1e-8)
})

test_that("a group without events is reported as a zero rate", {
  out <- negative_binomial_rate(matrix(0, nrow = 2, ncol = 5),
                                matrix(1, nrow = 2, ncol = 5),
                                dispersion = 0.8)
  expect_equal(out$rate, c(0, 0))
})

test_that("invalid arguments are rejected", {
  expect_error(negative_binomial_rate(matrix(1, 2, 2), matrix(1, 2, 3), 0.8),
               "same shape")
  expect_error(negative_binomial_rate(matrix(1, 2, 2), matrix(1, 2, 2), 0),
               "positive")
})
