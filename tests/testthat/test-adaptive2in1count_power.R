curve_at <- function(...) {
  arguments <- list(rate1 = 2.6, rate2 = 1.3, dispersion1 = 0.8,
                    dispersion2 = 0.8, shared = 0.4, n1 = 60, n2 = 90,
                    n3 = 300, window1 = 0.5, window2 = 1)
  do.call(adaptive2in1count_power, modifyList(arguments, list(...)))
}

test_that("the power is the sum of its two parts", {
  curve <- curve_at(cutpoint = c(0.6, 1), ratio1 = c(0.6, 0.8),
                    ratio2 = c(0.7, 0.85))
  expect_equal(curve$results$power,
               curve$results$power_small + curve$results$power_large,
               tolerance = 1e-10)
  expect_equal(nrow(curve$results), 8L)
})

test_that("the power increases as the Phase 3 effect strengthens", {
  curve <- curve_at(cutpoint = 10, ratio1 = 0.7,
                    ratio2 = c(0.6, 0.7, 0.85, 1))
  ordered <- curve$results[order(curve$results$ratio2), ]
  expect_true(all(diff(ordered$power) < 0))
})

test_that("the probability of expansion falls as the Phase 2 effect weakens", {
  curve <- curve_at(cutpoint = 0.8, ratio1 = c(0.5, 0.7, 0.9, 1), ratio2 = 0.75)
  ordered <- curve$results[order(curve$results$ratio1), ]
  expect_true(all(diff(ordered$go_probability) < 0))
})

test_that("the probability of expansion ignores the Phase 3 effect", {
  curve <- curve_at(cutpoint = 0.8, ratio1 = 0.7, ratio2 = c(0.6, 0.8, 1))
  expect_equal(diff(range(curve$results$go_probability)), 0, tolerance = 1e-12)
})

test_that("a cutpoint that always expands leaves no power in the small trial", {
  curve <- curve_at(cutpoint = 1e8, ratio1 = 0.7, ratio2 = 0.75)
  expect_equal(curve$results$go_probability, 1, tolerance = 1e-8)
  expect_equal(curve$results$power_small, 0, tolerance = 1e-8)
})

test_that("dropout lowers the power in every scenario", {
  complete <- curve_at(cutpoint = 0.8, ratio1 = c(0.6, 0.8),
                       ratio2 = c(0.7, 0.85))
  reduced <- curve_at(cutpoint = 0.8, ratio1 = c(0.6, 0.8),
                      ratio2 = c(0.7, 0.85), dropout = 0.4,
                      dropout_shape = "exponential")
  expect_true(all(reduced$results$power < complete$results$power))
})

test_that("the power exceeds the type I error rate", {
  null <- curve_at(cutpoint = 0.8, ratio1 = 1, ratio2 = 1)
  alternative <- curve_at(cutpoint = 0.8, ratio1 = 0.6, ratio2 = 0.7)
  expect_gt(alternative$results$power, null$results$power)
  expect_lte(null$results$power, 0.025 + 1e-9)
  expect_s3_class(alternative, "adaptive2in1count_power")
})

test_that("invalid arguments are rejected", {
  expect_error(curve_at(cutpoint = 0.8, ratio1 = -1, ratio2 = 0.8),
               "positive values")
  expect_error(curve_at(cutpoint = 0.8, ratio1 = 0.8, ratio2 = 0),
               "positive values")
})
