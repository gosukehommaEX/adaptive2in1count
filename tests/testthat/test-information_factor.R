test_that("complete follow-up matches the closed form", {
  for (rate in c(0.5, 1.3, 4)) {
    for (dispersion in c(0.4, 0.8, 2)) {
      for (window in c(0.5, 1, 2)) {
        mean_count <- rate * window
        expect_equal(information_factor(rate, dispersion, window, window),
                     dispersion * mean_count / (dispersion + mean_count),
                     tolerance = 1e-10)
      }
    }
  }
})

test_that("the reciprocal reproduces the familiar variance", {
  rate <- 1.3
  dispersion <- 0.8
  window <- 1
  size <- 200
  variance <- 1 / (size * information_factor(rate, dispersion, window, window))
  expect_equal(variance, 1 / (size * rate * window) + 1 / (size * dispersion),
               tolerance = 1e-10)
})

test_that("dropout reduces the information", {
  complete <- information_factor(1.3, 0.8, 1, 1)
  for (shape in c("uniform", "exponential")) {
    for (dropout in c(0.1, 0.3, 0.5)) {
      expect_lt(information_factor(1.3, 0.8, 1, 1, dropout, shape), complete)
    }
    values <- vapply(c(0.1, 0.3, 0.5, 0.7), function(d) {
      information_factor(1.3, 0.8, 1, 1, d, shape)
    }, numeric(1))
    expect_true(all(diff(values) < 0))
  }
})

test_that("the uniform follow-up period matches its closed form", {
  rate <- 1.3
  dispersion <- 0.8
  window <- 1
  mean_count <- rate * window
  complete <- dispersion * mean_count / (dispersion + mean_count)
  early <- dispersion -
    dispersion^2 / mean_count * log(1 + mean_count / dispersion)
  for (dropout in c(0.2, 0.5, 1)) {
    expect_equal(information_factor(rate, dispersion, window, window, dropout,
                                    "uniform"),
                 (1 - dropout) * complete + dropout * early, tolerance = 1e-7)
  }
})

test_that("a dropout proportion of zero reduces to complete follow-up", {
  reference <- information_factor(1.3, 0.8, 0.5, 1)
  for (shape in c("uniform", "exponential")) {
    expect_equal(information_factor(1.3, 0.8, 0.5, 1, 0, shape), reference,
                 tolerance = 1e-10)
  }
})

test_that("the exponential period carries less information than the uniform", {
  for (dropout in c(0.1, 0.3, 0.5)) {
    expect_lt(information_factor(1.3, 0.8, 1, 1, dropout, "exponential"),
              information_factor(1.3, 0.8, 1, 1, dropout, "uniform"))
  }
})

test_that("a shorter window carries less information", {
  expect_lt(information_factor(1.3, 0.8, 0.5, 1),
            information_factor(1.3, 0.8, 1, 1))
})

test_that("the information is bounded by the dispersion", {
  for (rate in c(0.2, 2, 50)) {
    expect_lt(information_factor(rate, 0.8, 1, 1), 0.8)
  }
})

test_that("invalid arguments are rejected", {
  expect_error(information_factor(0, 0.8, 1, 1), "positive")
  expect_error(information_factor(1.3, 0, 1, 1), "positive")
  expect_error(information_factor(1.3, 0.8, 2, 1), "must not exceed")
})
