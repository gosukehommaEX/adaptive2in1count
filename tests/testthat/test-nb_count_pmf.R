test_that("the probability mass function integrates to one", {
  for (dispersion in c(0.3, 0.5, 1, 2.5, 10)) {
    for (mean in c(0.5, 3, 20)) {
      ymax <- nb_support_max(dispersion, mean)
      expect_equal(sum(nb_count_pmf(dispersion, mean, ymax)), 1,
                   tolerance = 1e-9)
    }
  }
})

test_that("the first two moments match the negative binomial model", {
  for (dispersion in c(0.4, 1.5, 6)) {
    for (mean in c(1, 4, 15)) {
      ymax <- nb_support_max(dispersion, mean, tol = 1e-13)
      pmf <- nb_count_pmf(dispersion, mean, ymax)
      support <- seq(0, ymax)
      expect_equal(sum(pmf * support), mean, tolerance = 1e-7)
      variance <- sum(pmf * support^2) - sum(pmf * support)^2
      expect_equal(variance, mean * (1 + mean / dispersion), tolerance = 1e-5)
    }
  }
})

test_that("the probability mass function agrees with dnbinom", {
  for (dispersion in c(0.25, 2.5, 40)) {
    for (mean in c(0.7, 7)) {
      pmf <- nb_count_pmf(dispersion, mean, 80)
      reference <- stats::dnbinom(seq(0, 80), size = dispersion, mu = mean)
      expect_equal(pmf, reference, tolerance = 1e-11)
    }
  }
})

test_that("a large dispersion approaches the Poisson distribution", {
  pmf <- nb_count_pmf(dispersion = 1e6, mean = 4, ymax = 40)
  expect_equal(pmf, stats::dpois(seq(0, 40), 4), tolerance = 1e-4)
})

test_that("the probability mass decays beyond the mode", {
  pmf <- nb_count_pmf(dispersion = 8, mean = 6, ymax = 60)
  tail <- pmf[seq(21, 61)]
  expect_true(all(diff(tail) <= 0))
})

test_that("a degenerate dispersion or mean gives a point mass at zero", {
  expect_equal(nb_count_pmf(0, 3, 5), c(1, rep(0, 5)))
  expect_equal(nb_count_pmf(2, 0, 5), c(1, rep(0, 5)))
  expect_equal(nb_count_pmf(-1, 3, 3), c(1, rep(0, 3)))
})

test_that("the returned vector has the requested length", {
  for (ymax in c(0, 1, 10, 100)) {
    expect_length(nb_count_pmf(1.5, 3, ymax), ymax + 1)
  }
})

test_that("invalid arguments are rejected", {
  expect_error(nb_count_pmf(c(1, 2), 3, 10), "single numeric")
  expect_error(nb_count_pmf(1, c(3, 4), 10), "single numeric")
  expect_error(nb_count_pmf(1, 3, -1), "non-negative")
  expect_error(nb_count_pmf(NA_real_, 3, 10), "single numeric")
})
