test_that("the maximum type I error rate never exceeds the nominal level", {
  size <- adaptive2in1count_typeI(rate1 = c(1.5, 2.6, 4), rate2 = c(0.8, 1.3, 2),
                                  dispersion1 = 0.8, dispersion2 = 0.8,
                                  shared = 0.4, n1 = 60, n2 = 90, n3 = 300,
                                  cutpoint = c(0.3, 0.6, 1, 2, 4),
                                  window1 = 0.5, window2 = 1, dropout = 0.3,
                                  dropout_shape = "exponential")
  expect_true(all(size$results$type1 <= 0.025 + 1e-9))
  expect_lte(size$overall$type1, 0.025 + 1e-9)
  expect_equal(nrow(size$results), 45L)
  expect_equal(nrow(size$maximum), 5L)
})

test_that("the level is respected at other significance levels", {
  for (alpha in c(0.01, 0.05)) {
    size <- adaptive2in1count_typeI(rate1 = c(1.5, 4), rate2 = c(0.8, 2),
                                    dispersion1 = 0.8, dispersion2 = 0.8,
                                    shared = 0.4, n1 = 60, n2 = 90, n3 = 300,
                                    cutpoint = c(0.5, 1, 3), alpha = alpha,
                                    window1 = 0.5, window2 = 1)
    expect_lte(size$overall$type1, alpha + 1e-9)
  }
})

test_that("the maximum agrees with direct calls", {
  grid <- expand.grid(rate1 = c(1.5, 4), rate2 = c(0.8, 2))
  direct <- vapply(seq_len(nrow(grid)), function(i) {
    adaptive2in1count(rate1 = grid$rate1[i], rate2 = grid$rate2[i],
                      dispersion1 = 0.8, dispersion2 = 0.8, shared = 0.4,
                      n1 = 60, n2 = 90, n3 = 300, cutpoint = 0.7,
                      window1 = 0.5, window2 = 1)$results$rejection
  }, numeric(1))
  size <- adaptive2in1count_typeI(rate1 = c(1.5, 4), rate2 = c(0.8, 2),
                                  dispersion1 = 0.8, dispersion2 = 0.8,
                                  shared = 0.4, n1 = 60, n2 = 90, n3 = 300,
                                  cutpoint = 0.7, window1 = 0.5, window2 = 1)
  expect_equal(sort(size$results$type1), sort(direct), tolerance = 1e-12)
  expect_equal(size$overall$type1, max(direct), tolerance = 1e-12)
})

test_that("the cutpoints are reported in increasing order", {
  size <- adaptive2in1count_typeI(rate1 = 2.6, rate2 = 1.3, dispersion1 = 0.8,
                                  dispersion2 = 0.8, shared = 0.4, n1 = 60,
                                  n2 = 90, n3 = 300, cutpoint = c(2, 0.5, 1),
                                  window1 = 0.5, window2 = 1)
  expect_equal(size$maximum$cutpoint, c(0.5, 1, 2))
  expect_s3_class(size, "adaptive2in1count_typeI")
})

test_that("a Phase 3 small trial also controls the type I error rate", {
  size <- adaptive2in1count_typeI(rate1 = c(1.5, 4), rate2 = c(0.8, 2),
                                  dispersion1 = 0.8, dispersion2 = 0.8,
                                  shared = 0.4, n1 = 60, n2 = 90, n3 = 300,
                                  cutpoint = c(0.4, 1, 3), window1 = 0.5,
                                  window2 = 1, dropout = 0.4,
                                  dropout_shape = "exponential",
                                  small_endpoint = "phase3")
  expect_lte(size$overall$type1, 0.025 + 1e-9)
})

test_that("invalid arguments are rejected", {
  expect_error(adaptive2in1count_typeI(-1, 1.3, 0.8, 0.8, 0.4, 60, 90, 300,
                                       0.8), "positive values")
  expect_error(adaptive2in1count_typeI(2.6, 0, 0.8, 0.8, 0.4, 60, 90, 300,
                                       0.8), "positive values")
})
