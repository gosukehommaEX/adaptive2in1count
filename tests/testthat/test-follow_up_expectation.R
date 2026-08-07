test_that("complete follow-up evaluates the function at the window", {
  expect_equal(follow_up_expectation(function(t) t, window = 2,
                                     shape = "complete"), 2)
  expect_equal(follow_up_expectation(function(t) t^2, window = 3,
                                     shape = "complete"), 9)
})

test_that("a dropout proportion of zero reduces to complete follow-up", {
  for (shape in c("uniform", "exponential")) {
    expect_equal(follow_up_expectation(function(t) t, window = 2, dropout = 0,
                                       shape = shape), 2)
    expect_equal(follow_up_expectation(function(t) t^2, window = 2,
                                       dropout = 0, shape = shape), 4)
  }
})

test_that("the uniform period is a mixture with the complete period", {
  window <- 2
  for (dropout in c(0.2, 0.5, 0.8, 1)) {
    expect_equal(follow_up_expectation(function(t) t, window, dropout,
                                       "uniform"),
                 (1 - dropout) * window + dropout * window / 2,
                 tolerance = 1e-8)
    expect_equal(follow_up_expectation(function(t) t^2, window, dropout,
                                       "uniform"),
                 (1 - dropout) * window^2 + dropout * window^2 / 3,
                 tolerance = 1e-8)
  }
})

test_that("both random periods attain the stated dropout proportion", {
  for (window in c(0.5, 1, 2)) {
    for (dropout in c(0.1, 0.3, 0.5)) {
      for (shape in c("uniform", "exponential")) {
        leaving <- follow_up_expectation(function(t) as.numeric(t < window),
                                         window, dropout, shape)
        expect_equal(leaving, dropout, tolerance = 1e-6)
      }
    }
  }
})

test_that("the mean period decreases with the dropout proportion", {
  for (shape in c("uniform", "exponential")) {
    values <- vapply(c(0.05, 0.2, 0.4, 0.6), function(d) {
      follow_up_expectation(function(t) t, window = 1, dropout = d,
                            shape = shape)
    }, numeric(1))
    expect_true(all(diff(values) < 0))
    expect_true(all(values < 1))
  }
})

test_that("the exponential period is shorter than the uniform period", {
  for (dropout in c(0.1, 0.3, 0.5, 0.7)) {
    uniform <- follow_up_expectation(function(t) t, 1, dropout, "uniform")
    exponential <- follow_up_expectation(function(t) t, 1, dropout,
                                         "exponential")
    expect_lt(exponential, uniform)
  }
})

test_that("a small dropout proportion approaches complete follow-up", {
  for (shape in c("uniform", "exponential")) {
    expect_equal(follow_up_expectation(function(t) t, window = 1,
                                       dropout = 1e-8, shape = shape),
                 1, tolerance = 1e-6)
  }
})

test_that("invalid arguments are rejected", {
  expect_error(follow_up_expectation(function(t) t, window = 0), "positive")
  expect_error(follow_up_expectation(function(t) t, window = 1, dropout = 1,
                                     shape = "exponential"),
               "strictly below one")
  expect_error(follow_up_expectation(function(t) t, window = 1, dropout = -0.1,
                                     shape = "uniform"), "between 0 and 1")
  expect_error(follow_up_expectation(function(t) t, window = 1, dropout = 1.5,
                                     shape = "uniform"), "between 0 and 1")
  expect_error(follow_up_expectation(function(t) t, window = 1,
                                     shape = "weibull"), "should be one of")
})
