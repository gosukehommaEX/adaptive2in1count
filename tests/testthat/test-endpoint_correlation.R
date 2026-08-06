test_that("complete follow-up without a shared frailty gives independence", {
  expect_equal(endpoint_correlation(2.6, 1.3, 0.8, 0.8, 0, 0.5, 1), 0,
               tolerance = 1e-10)
})

test_that("a random follow-up period correlates the endpoints on its own", {
  for (shape in c("uniform", "exponential")) {
    value <- endpoint_correlation(2.6, 1.3, 0.8, 0.8, shared = 0,
                                  window1 = 0.5, window2 = 1, dropout = 0.4,
                                  dropout_shape = shape)
    expect_gt(value, 0)
  }
})

test_that("the correlation increases with the shared frailty shape", {
  values <- vapply(c(0, 0.2, 0.5, 0.8), function(s) {
    endpoint_correlation(2.6, 1.3, 0.8, 0.8, s, 0.5, 1)
  }, numeric(1))
  expect_true(all(diff(values) > 0))
})

test_that("the correlation increases with the dropout proportion", {
  values <- vapply(c(0.02, 0.2, 0.4, 0.6), function(d) {
    endpoint_correlation(2.6, 1.3, 0.8, 0.8, 0.4, 0.5, 1, d, "exponential")
  }, numeric(1))
  expect_true(all(diff(values) > 0))
})

test_that("the correlation stays within the unit interval", {
  grid <- expand.grid(rate1 = c(0.3, 2.6), rate2 = c(0.3, 2.6),
                      dispersion = c(0.3, 1, 3), dropout = c(0, 0.3, 0.6))
  for (i in seq_len(nrow(grid))) {
    shape <- if (grid$dropout[i] == 0) "complete" else "exponential"
    value <- endpoint_correlation(grid$rate1[i], grid$rate2[i],
                                  grid$dispersion[i], grid$dispersion[i],
                                  grid$dispersion[i], 0.5, 1, grid$dropout[i],
                                  shape)
    expect_gte(value, 0)
    expect_lte(value, 1)
  }
})

test_that("identical endpoints under complete follow-up are fully correlated", {
  value <- endpoint_correlation(rate1 = 1.3, rate2 = 1.3, dispersion1 = 0.8,
                                dispersion2 = 0.8, shared = 0.8,
                                window1 = 1, window2 = 1)
  expect_gt(value, 0.5)
  expect_lte(value, 1)
})

test_that("invalid arguments are rejected", {
  expect_error(endpoint_correlation(2.6, 1.3, 0.8, 0.8, 2, 0.5, 1),
               "smaller of the two dispersions")
})
