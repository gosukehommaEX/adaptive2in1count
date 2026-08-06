base <- list(rate1 = 2.6, rate2 = 1.3, dispersion1 = 0.8, dispersion2 = 0.8,
             shared = 0.4, n1 = 60, n2 = 90, n3 = 300, window1 = 0.5,
             window2 = 1)

test_that("the requirement is met across a grid", {
  grid <- expand.grid(dispersion = c(0.3, 0.8, 2), dropout = c(0, 0.3, 0.6),
                      n3 = c(90, 300, 900), endpoint = c("phase2", "phase3"),
                      stringsAsFactors = FALSE)
  for (i in seq_len(nrow(grid))) {
    shape <- if (grid$dropout[i] == 0) "complete" else "exponential"
    rho <- adaptive2in1count_correlation(
      rate1 = 2.6, rate2 = 1.3, dispersion1 = grid$dispersion[i],
      dispersion2 = grid$dispersion[i], shared = grid$dispersion[i], n1 = 60,
      n2 = 90, n3 = grid$n3[i], window1 = 0.5, window2 = 1,
      dropout = grid$dropout[i], dropout_shape = shape,
      small_endpoint = grid$endpoint[i])
    expect_true(rho$holds)
  }
})

test_that("a Phase 2 small trial gives the nesting correlation", {
  rho <- do.call(adaptive2in1count_correlation, base)
  expect_equal(rho$interim_small, sqrt(60 / 90), tolerance = 1e-9)
})

test_that("an interim at the full small trial gives a correlation of one", {
  rho <- do.call(adaptive2in1count_correlation,
                 modifyList(base, list(n1 = 90)))
  expect_equal(rho$interim_small, 1, tolerance = 1e-9)
  expect_true(rho$holds)
})

test_that("independent endpoints under complete follow-up decouple", {
  rho <- do.call(adaptive2in1count_correlation,
                 modifyList(base, list(shared = 0)))
  expect_equal(rho$interim_large, 0, tolerance = 1e-10)
  expect_equal(rho$endpoint, 0, tolerance = 1e-10)
  expect_true(rho$holds)
})

test_that("dropout couples the endpoints without a shared frailty", {
  rho <- do.call(adaptive2in1count_correlation,
                 modifyList(base, list(shared = 0, dropout = 0.4,
                                       dropout_shape = "exponential")))
  expect_gt(rho$interim_large, 0)
  expect_gt(rho$endpoint, 0)
  expect_true(rho$holds)
})

test_that("the correlation with the large trial falls as it grows", {
  values <- vapply(c(150, 300, 600, 1200), function(size) {
    do.call(adaptive2in1count_correlation,
            modifyList(base, list(n3 = size)))$interim_large
  }, numeric(1))
  expect_true(all(diff(values) < 0))
})

test_that("unequal allocation is accepted", {
  rho <- do.call(adaptive2in1count_correlation,
                 modifyList(base, list(n1 = c(80, 40), n2 = c(120, 60),
                                       n3 = c(400, 200))))
  expect_true(rho$holds)
  expect_s3_class(rho, "adaptive2in1count_correlation")
})

test_that("invalid arguments are rejected", {
  expect_error(do.call(adaptive2in1count_correlation,
                       modifyList(base, list(n1 = 200))), "non-decreasing")
})
