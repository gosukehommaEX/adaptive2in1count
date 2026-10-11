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
