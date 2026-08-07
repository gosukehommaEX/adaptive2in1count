# Closed forms for the exposures of the two endpoints, used as the reference
# throughout this file. The uniform period is the two point mixture of Tang
# (2015): a proportion d leaves at a time uniform on the window and the rest are
# followed to the end.
uniform_moments <- function(window1, window2, dropout) {
  complete <- c(first = window1, second = window2, cross = window1 * window2,
                square1 = window1^2, square2 = window2^2)
  early <- c(
    first = window1 - window1^2 / (2 * window2),
    second = window2 / 2,
    cross = (window1^3 / 3 + window1 * (window2^2 - window1^2) / 2) / window2,
    square1 = (window1^3 / 3 + window1^2 * (window2 - window1)) / window2,
    square2 = window2^2 / 3)
  (1 - dropout) * complete + dropout * early
}

test_that("complete follow-up gives deterministic exposures", {
  out <- follow_up_moments(window1 = 0.5, window2 = 1, shape = "complete")
  expect_equal(out$first, 0.5)
  expect_equal(out$second, 1)
  expect_equal(out$cross, 0.5)
  expect_equal(out$square1, 0.25)
  expect_equal(out$square2, 1)
  expect_equal(out$inflation, 1)
})

test_that("the uniform period matches the closed form of the mixture", {
  for (dropout in c(0.2, 0.5, 1)) {
    for (windows in list(c(0.5, 2), c(0.25, 1), c(1, 1))) {
      out <- follow_up_moments(windows[1], windows[2], dropout, "uniform")
      reference <- uniform_moments(windows[1], windows[2], dropout)
      expect_equal(out$first, reference[["first"]], tolerance = 1e-7)
      expect_equal(out$second, reference[["second"]], tolerance = 1e-7)
      expect_equal(out$cross, reference[["cross"]], tolerance = 1e-7)
      expect_equal(out$square1, reference[["square1"]], tolerance = 1e-7)
      expect_equal(out$square2, reference[["square2"]], tolerance = 1e-7)
    }
  }
})

test_that("the exponential period matches its closed form", {
  window1 <- 0.5
  window2 <- 2
  dropout <- 0.4
  lambda <- -log(1 - dropout) / window2
  out <- follow_up_moments(window1, window2, dropout, "exponential")
  expect_equal(out$first, (1 - exp(-lambda * window1)) / lambda,
               tolerance = 1e-7)
  expect_equal(out$second, (1 - exp(-lambda * window2)) / lambda,
               tolerance = 1e-7)
  expect_equal(out$square2,
               2 / lambda^2 *
                 (1 - exp(-lambda * window2) * (1 + lambda * window2)),
               tolerance = 1e-7)
})

test_that("a random follow-up period inflates the product moment", {
  for (shape in c("uniform", "exponential")) {
    out <- follow_up_moments(0.5, 1, dropout = 0.4, shape = shape)
    expect_gt(out$inflation, 1)
    expect_gt(out$cross, out$first * out$second)
  }
})

test_that("the inflation increases with the dropout proportion", {
  for (shape in c("uniform", "exponential")) {
    values <- vapply(c(0.05, 0.2, 0.4, 0.6), function(d) {
      follow_up_moments(0.5, 1, dropout = d, shape = shape)$inflation
    }, numeric(1))
    expect_true(all(diff(values) > 0))
  }
})

test_that("a dropout proportion of zero reduces to complete follow-up", {
  reference <- follow_up_moments(0.5, 1, shape = "complete")
  for (shape in c("uniform", "exponential")) {
    out <- follow_up_moments(0.5, 1, dropout = 0, shape = shape)
    expect_equal(out$inflation, 1, tolerance = 1e-10)
    expect_equal(out$cross, reference$cross, tolerance = 1e-10)
  }
})

test_that("equal windows make the cross moment the second moment", {
  for (shape in c("uniform", "exponential")) {
    out <- follow_up_moments(1, 1, dropout = 0.3, shape = shape)
    expect_equal(out$cross, out$square2, tolerance = 1e-8)
    expect_equal(out$first, out$second, tolerance = 1e-8)
  }
})

test_that("the second moments dominate the squared means", {
  for (shape in c("uniform", "exponential")) {
    out <- follow_up_moments(0.5, 1, dropout = 0.3, shape = shape)
    expect_gt(out$square1, out$first^2)
    expect_gt(out$square2, out$second^2)
  }
})

test_that("invalid arguments are rejected", {
  expect_error(follow_up_moments(window1 = 2, window2 = 1), "not exceeding")
  expect_error(follow_up_moments(window1 = 0, window2 = 1), "positive")
})
