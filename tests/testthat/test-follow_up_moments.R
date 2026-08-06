test_that("complete follow-up gives deterministic exposures", {
  out <- follow_up_moments(window1 = 0.5, window2 = 1, shape = "complete")
  expect_equal(out$first, 0.5)
  expect_equal(out$second, 1)
  expect_equal(out$cross, 0.5)
  expect_equal(out$inflation, 1)
})

test_that("the uniform follow-up period matches the closed form", {
  window1 <- 0.5
  window2 <- 2
  out <- follow_up_moments(window1, window2, shape = "uniform")
  expect_equal(out$first, window1 - window1^2 / (2 * window2),
               tolerance = 1e-7)
  expect_equal(out$second, window2 / 2, tolerance = 1e-7)
  expected_cross <- (window1^3 / 3 +
                       window1 * (window2^2 - window1^2) / 2) / window2
  expect_equal(out$cross, expected_cross, tolerance = 1e-7)
})

test_that("a random follow-up period inflates the product moment", {
  for (shape in c("uniform", "exponential")) {
    out <- follow_up_moments(0.5, 1, dropout = 0.4, shape = shape)
    expect_gt(out$inflation, 1)
    expect_gt(out$cross, out$first * out$second)
  }
})

test_that("the inflation increases with the dropout proportion", {
  values <- vapply(c(0.05, 0.2, 0.4, 0.6), function(d) {
    follow_up_moments(0.5, 1, dropout = d, shape = "exponential")$inflation
  }, numeric(1))
  expect_true(all(diff(values) > 0))
})

test_that("equal windows make the cross moment the second moment", {
  out <- follow_up_moments(1, 1, dropout = 0.3, shape = "exponential")
  expect_equal(out$cross, out$square2, tolerance = 1e-8)
  expect_equal(out$first, out$second, tolerance = 1e-8)
})

test_that("the second moments dominate the squared means", {
  out <- follow_up_moments(0.5, 1, dropout = 0.3, shape = "exponential")
  expect_gt(out$square1, out$first^2)
  expect_gt(out$square2, out$second^2)
})

test_that("invalid arguments are rejected", {
  expect_error(follow_up_moments(window1 = 2, window2 = 1), "not exceeding")
  expect_error(follow_up_moments(window1 = 0, window2 = 1), "positive")
})
