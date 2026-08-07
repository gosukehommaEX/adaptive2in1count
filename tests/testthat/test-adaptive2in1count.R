design <- function(...) {
  arguments <- list(rate1 = 2.6, rate2 = 1.3, dispersion1 = 0.8,
                    dispersion2 = 0.8, shared = 0.4, n1 = 60, n2 = 90,
                    n3 = 300, window1 = 0.5, window2 = 1)
  do.call(adaptive2in1count, modifyList(arguments, list(...)))
}

test_that("the rejection probability is the sum of its two parts", {
  fit <- design(cutpoint = c(0.4, 0.7, 1, 1.6))
  expect_equal(fit$results$rejection,
               fit$results$rejection_small + fit$results$rejection_large,
               tolerance = 1e-10)
})

test_that("a cutpoint that never expands attains the nominal level", {
  fit <- design(cutpoint = 1e-8, alpha = 0.025)
  expect_equal(fit$results$go_probability, 0, tolerance = 1e-8)
  expect_equal(fit$results$rejection, 0.025, tolerance = 1e-6)
  expect_equal(fit$results$rejection_large, 0, tolerance = 1e-8)
})

test_that("a cutpoint that always expands attains the nominal level", {
  fit <- design(cutpoint = 1e8, alpha = 0.025)
  expect_equal(fit$results$go_probability, 1, tolerance = 1e-8)
  expect_equal(fit$results$rejection, 0.025, tolerance = 1e-6)
  expect_equal(fit$results$rejection_small, 0, tolerance = 1e-8)
})

test_that("the type I error rate never exceeds the nominal level", {
  grid <- expand.grid(dispersion = c(0.3, 0.8, 2), dropout = c(0, 0.2, 0.5),
                      endpoint = c("phase2", "phase3"),
                      stringsAsFactors = FALSE)
  cutpoints <- c(0.2, 0.4, 0.6, 0.8, 1, 1.4, 2, 4)
  for (i in seq_len(nrow(grid))) {
    shape <- if (grid$dropout[i] == 0) "complete" else "exponential"
    fit <- design(cutpoint = cutpoints, dispersion1 = grid$dispersion[i],
                  dispersion2 = grid$dispersion[i],
                  shared = grid$dispersion[i], dropout = grid$dropout[i],
                  dropout_shape = shape, small_endpoint = grid$endpoint[i])
    expect_true(all(fit$results$rejection <= 0.025 + 1e-9))
  }
})

test_that("the level is respected at other significance levels", {
  for (alpha in c(0.005, 0.01, 0.05, 0.1)) {
    fit <- design(cutpoint = c(0.3, 0.6, 1, 2), alpha = alpha,
                  dropout = 0.3, dropout_shape = "exponential")
    expect_true(all(fit$results$rejection <= alpha + 1e-9))
  }
})

test_that("the probability of expansion increases with the cutpoint", {
  fit <- design(cutpoint = c(0.2, 0.4, 0.6, 0.9, 1.3, 2))
  expect_true(all(diff(fit$results$go_probability) > 0))
  expect_true(all(diff(fit$results$expected_size) > 0))
})

test_that("the expected sample size lies between the two trial sizes", {
  fit <- design(cutpoint = c(0.3, 0.8, 2))
  expect_true(all(fit$results$expected_size >= 180))
  expect_true(all(fit$results$expected_size <= 600))
})

test_that("a stronger treatment effect gives a larger rejection probability", {
  null <- design(cutpoint = 0.8)
  medium <- design(cutpoint = 0.8, ratio1 = 0.8, ratio2 = 0.85)
  strong <- design(cutpoint = 0.8, ratio1 = 0.6, ratio2 = 0.7)
  expect_lt(null$results$rejection, medium$results$rejection)
  expect_lt(medium$results$rejection, strong$results$rejection)
  expect_lt(null$results$go_probability, strong$results$go_probability)
})

test_that("dropout lowers the power", {
  values <- vapply(c(0.02, 0.2, 0.4, 0.6), function(d) {
    design(cutpoint = 0.8, ratio1 = 0.65, ratio2 = 0.75, dropout = d,
           dropout_shape = "exponential")$results$rejection
  }, numeric(1))
  expect_true(all(diff(values) < 0))
})

test_that("a negligible dropout proportion reproduces complete follow-up", {
  complete <- design(cutpoint = 0.8, ratio1 = 0.65, ratio2 = 0.75)
  nearly <- design(cutpoint = 0.8, ratio1 = 0.65, ratio2 = 0.75,
                   dropout = 1e-9, dropout_shape = "exponential")
  expect_equal(complete$results$rejection, nearly$results$rejection,
               tolerance = 1e-5)
})

test_that("the exponential period is more severe than the uniform period", {
  for (dropout in c(0.2, 0.4)) {
    uniform <- design(cutpoint = 0.8, ratio1 = 0.65, ratio2 = 0.75,
                      dropout = dropout, dropout_shape = "uniform")
    exponential <- design(cutpoint = 0.8, ratio1 = 0.65, ratio2 = 0.75,
                          dropout = dropout, dropout_shape = "exponential")
    expect_lt(exponential$results$rejection, uniform$results$rejection)
  }
})

test_that("either random period at no dropout reproduces complete follow-up", {
  reference <- design(cutpoint = c(0.6, 1), ratio1 = 0.65, ratio2 = 0.75)
  for (shape in c("uniform", "exponential")) {
    matched <- design(cutpoint = c(0.6, 1), ratio1 = 0.65, ratio2 = 0.75,
                      dropout = 0, dropout_shape = shape)
    expect_equal(matched$results$rejection, reference$results$rejection,
                 tolerance = 1e-10)
  }
})

test_that("a larger large trial gives more power once the trial expands", {
  values <- vapply(c(200, 400, 800), function(size) {
    design(cutpoint = 4, n3 = size, ratio1 = 0.65,
           ratio2 = 0.75)$results$rejection
  }, numeric(1))
  expect_true(all(diff(values) > 0))
})

test_that("unequal allocation is accepted and controls the type I error rate", {
  for (allocation in list(c(2, 1), c(1, 2), c(3, 1))) {
    fit <- design(cutpoint = c(0.5, 1, 2), n1 = 30 * allocation,
                  n2 = 45 * allocation, n3 = 150 * allocation)
    expect_equal(fit$n1, 30 * allocation)
    expect_true(all(fit$results$rejection <= 0.025 + 1e-9))
  }
})

test_that("the returned object carries its inputs and class", {
  fit <- design(cutpoint = 0.8, alpha = 0.05, ratio1 = 0.7, dropout = 0.3,
                dropout_shape = "exponential", small_endpoint = "phase3")
  expect_s3_class(fit, "adaptive2in1count")
  expect_equal(fit$alpha, 0.05)
  expect_equal(fit$ratio1, 0.7)
  expect_equal(fit$dropout, c(0.3, 0.3))
  expect_equal(fit$small_endpoint, "phase3")
  expect_named(fit$components$variance, c("interim", "small", "large"))
})

test_that("invalid arguments are rejected", {
  expect_error(design(cutpoint = 0.8, rate1 = -1), "positive")
  expect_error(design(cutpoint = 0.8, shared = 2),
               "smaller of the two dispersions")
  expect_error(design(cutpoint = 0.8, n1 = 200), "non-decreasing")
  expect_error(design(cutpoint = -1), "positive values")
  expect_error(design(cutpoint = 0.8, alpha = 1.2), "between 0 and 1")
  expect_error(design(cutpoint = 0.8, dropout_shape = "weibull"),
               "should be one of")
})
