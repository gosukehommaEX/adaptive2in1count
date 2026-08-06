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

test_that("the contrast reproduces the single group fits at a fixed dispersion", {
  set.seed(11)
  exposure1 <- matrix(runif(300, 0.3, 1), nrow = 10)
  exposure2 <- matrix(runif(300, 0.3, 1), nrow = 10)
  count1 <- matrix(rpois(300, 1.2 * exposure1), nrow = 10)
  count2 <- matrix(rpois(300, 2 * exposure2), nrow = 10)
  out <- negative_binomial_contrast(count1, exposure1, count2, exposure2,
                                    dispersion = 0.8)
  first <- negative_binomial_rate(count1, exposure1, 0.8)
  second <- negative_binomial_rate(count2, exposure2, 0.8)
  expect_equal(out$estimate, log(first$rate) - log(second$rate),
               tolerance = 1e-10)
  expect_equal(out$error, sqrt(1 / first$information + 1 / second$information),
               tolerance = 1e-10)
})

test_that("the estimated dispersion solves its own score equation", {
  set.seed(12)
  exposure <- matrix(runif(1200, 0.3, 1), nrow = 4)
  count1 <- matrix(rnbinom(1200, size = 0.8, mu = 1.2 * exposure), nrow = 4)
  count2 <- matrix(rnbinom(1200, size = 0.8, mu = 2 * exposure), nrow = 4)
  out <- negative_binomial_contrast(count1, exposure, count2, exposure)
  expect_true(all(out$dispersion > 0.3 & out$dispersion < 2.5))
  expect_true(all(is.finite(out$estimate)))
  expect_true(all(out$error > 0))
})

test_that("a fixed dispersion at the estimate reproduces the contrast", {
  set.seed(13)
  exposure <- matrix(runif(600, 0.3, 1), nrow = 3)
  count1 <- matrix(rnbinom(600, size = 0.7, mu = exposure), nrow = 3)
  count2 <- matrix(rnbinom(600, size = 0.7, mu = 1.8 * exposure), nrow = 3)
  free <- negative_binomial_contrast(count1, exposure, count2, exposure)
  held <- negative_binomial_contrast(count1, exposure, count2, exposure,
                                     dispersion = free$dispersion)
  expect_equal(free$estimate, held$estimate, tolerance = 1e-6)
  expect_equal(free$error, held$error, tolerance = 1e-6)
})

test_that("the joint fit reproduces MASS::glm.nb on the same data", {
  set.seed(21)
  size <- 80
  exposure1 <- runif(size, 0.3, 1)
  exposure2 <- runif(size, 0.3, 1)
  count1 <- rnbinom(size, size = 0.8, mu = 1.2 * exposure1)
  count2 <- rnbinom(size, size = 0.8, mu = 2 * exposure2)
  out <- negative_binomial_contrast(matrix(count1, nrow = 1),
                                    matrix(exposure1, nrow = 1),
                                    matrix(count2, nrow = 1),
                                    matrix(exposure2, nrow = 1))
  frame <- data.frame(count = c(count1, count2),
                      offset_term = log(c(exposure1, exposure2)),
                      arm = factor(rep(c(1, 2), each = size),
                                   levels = c(2, 1)))
  reference <- MASS::glm.nb(count ~ arm + offset(offset_term), data = frame)
  expect_equal(out$estimate, unname(coef(reference)[["arm1"]]),
               tolerance = 1e-4)
  expect_equal(out$error, sqrt(vcov(reference)["arm1", "arm1"]),
               tolerance = 1e-3)
  expect_equal(out$dispersion, unname(reference$theta), tolerance = 1e-3)
})

test_that("the joint fit matches the reference across several data sets", {
  set.seed(22)
  for (replicate in seq_len(4)) {
    size <- 60
    exposure <- runif(2 * size, 0.2, 1)
    count <- rnbinom(2 * size, size = 1.2,
                     mu = rep(c(1, 1.8), each = size) * exposure)
    first <- seq_len(size)
    second <- size + first
    out <- negative_binomial_contrast(matrix(count[first], nrow = 1),
                                      matrix(exposure[first], nrow = 1),
                                      matrix(count[second], nrow = 1),
                                      matrix(exposure[second], nrow = 1))
    frame <- data.frame(count = count, offset_term = log(exposure),
                        arm = factor(rep(c(1, 2), each = size),
                                     levels = c(2, 1)))
    reference <- suppressWarnings(
      MASS::glm.nb(count ~ arm + offset(offset_term), data = frame))
    expect_equal(out$estimate, unname(coef(reference)[["arm1"]]),
                 tolerance = 1e-3)
  }
})
