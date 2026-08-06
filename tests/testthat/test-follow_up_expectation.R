test_that("complete follow-up evaluates the function at the window", {
  expect_equal(follow_up_expectation(function(t) t, window = 2,
                                     shape = "complete"), 2)
  expect_equal(follow_up_expectation(function(t) t^2, window = 3,
                                     shape = "complete"), 9)
})

test_that("the uniform follow-up period has the expected moments", {
  window <- 2
  expect_equal(follow_up_expectation(function(t) t, window, shape = "uniform"),
               window / 2, tolerance = 1e-8)
  expect_equal(follow_up_expectation(function(t) t^2, window,
                                     shape = "uniform"),
               window^2 / 3, tolerance = 1e-8)
})

test_that("the truncated exponential follow-up period attains the dropout", {
  for (window in c(0.5, 1, 2)) {
    for (dropout in c(0.1, 0.3, 0.5)) {
      leaving <- follow_up_expectation(function(t) as.numeric(t < window),
                                       window, dropout, "exponential")
      expect_equal(leaving, dropout, tolerance = 1e-6)
    }
  }
})

test_that("the mean follow-up period decreases with the dropout proportion", {
  values <- vapply(c(0.05, 0.2, 0.4, 0.6), function(d) {
    follow_up_expectation(function(t) t, window = 1, dropout = d,
                          shape = "exponential")
  }, numeric(1))
  expect_true(all(diff(values) < 0))
  expect_true(all(values < 1))
})

test_that("a small dropout proportion approaches complete follow-up", {
  expect_equal(follow_up_expectation(function(t) t, window = 1,
                                     dropout = 1e-8, shape = "exponential"),
               1, tolerance = 1e-6)
})

test_that("invalid arguments are rejected", {
  expect_error(follow_up_expectation(function(t) t, window = 0), "positive")
  expect_error(follow_up_expectation(function(t) t, window = 1, dropout = 0,
                                     shape = "exponential"), "between 0")
  expect_error(follow_up_expectation(function(t) t, window = 1, dropout = 1,
                                     shape = "exponential"), "between 0")
  expect_error(follow_up_expectation(function(t) t, window = 1,
                                     shape = "weibull"), "should be one of")
})
