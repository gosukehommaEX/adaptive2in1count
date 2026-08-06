base <- list(rate1 = 2.6, rate2 = 1.3, dispersion1 = 0.8, dispersion2 = 0.8,
             shared = 0.4, n1 = 60, n2 = 90, n3 = 300, window1 = 0.5,
             window2 = 1)

test_that("the cutpoint reproduces the target probability of expansion", {
  for (target in c(0.05, 0.2, 0.5, 0.8, 0.95)) {
    for (ratio1 in c(0.6, 1)) {
      out <- do.call(adaptive2in1count_cutpoint,
                     c(base, list(target = target, ratio1 = ratio1)))
      fit <- do.call(adaptive2in1count,
                     c(base, list(cutpoint = out$cutpoint, ratio1 = ratio1)))
      expect_equal(fit$results$go_probability, target, tolerance = 1e-8)
    }
  }
})

test_that("a larger target gives a larger cutpoint", {
  cutpoints <- vapply(c(0.1, 0.3, 0.5, 0.7, 0.9), function(t) {
    do.call(adaptive2in1count_cutpoint, c(base, list(target = t)))$cutpoint
  }, numeric(1))
  expect_true(all(diff(cutpoints) > 0))
})

test_that("a probability of one half sits at the true rate ratio", {
  out <- do.call(adaptive2in1count_cutpoint,
                 c(base, list(target = 0.5, ratio1 = 0.65)))
  expect_equal(out$cutpoint, 0.65, tolerance = 1e-8)
})

test_that("dropout widens the cutpoint needed for a high target", {
  complete <- do.call(adaptive2in1count_cutpoint,
                      c(base, list(target = 0.8, ratio1 = 0.65)))
  reduced <- do.call(adaptive2in1count_cutpoint,
                     c(base, list(target = 0.8, ratio1 = 0.65, dropout = 0.4,
                                  dropout_shape = "exponential")))
  expect_gt(reduced$cutpoint, complete$cutpoint)
})

test_that("invalid arguments are rejected", {
  expect_error(do.call(adaptive2in1count_cutpoint, c(base, list(target = 0))),
               "between 0 and 1")
  expect_error(do.call(adaptive2in1count_cutpoint, c(base, list(target = 1))),
               "between 0 and 1")
})
