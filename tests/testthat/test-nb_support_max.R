test_that("the discarded upper tail is below the tolerance", {
  settings <- expand.grid(dispersion = c(0.3, 1, 5, 40), mean = c(0.5, 4, 60))
  for (i in seq_len(nrow(settings))) {
    dispersion <- settings$dispersion[i]
    mean <- settings$mean[i]
    ymax <- nb_support_max(dispersion, mean, tol = 1e-11)
    expect_lt(1 - sum(nb_count_pmf(dispersion, mean, ymax)), 1e-11)
  }
})

test_that("the support is the smallest one that meets the tolerance", {
  dispersion <- 2
  mean <- 10
  tol <- 1e-9
  ymax <- nb_support_max(dispersion, mean, tol = tol)
  expect_gte(1 - sum(nb_count_pmf(dispersion, mean, ymax - 1)), tol)
})

test_that("a tighter tolerance never shortens the support", {
  loose <- nb_support_max(25, 150, tol = 1e-8)
  tight <- nb_support_max(25, 150, tol = 1e-12)
  expect_gte(tight, loose)
})

test_that("the support grows with the mean and with the overdispersion", {
  expect_gt(nb_support_max(5, 60), nb_support_max(5, 6))
  expect_gt(nb_support_max(0.3, 20), nb_support_max(5, 20))
})

test_that("the cap is respected", {
  expect_lte(nb_support_max(0.1, 5, tol = 1e-11, cap = 50L), 50L)
})

test_that("a degenerate distribution has a support of zero", {
  expect_equal(nb_support_max(0, 5), 0L)
  expect_equal(nb_support_max(5, 0), 0L)
})

test_that("an integer is returned", {
  expect_type(nb_support_max(2, 10), "integer")
})

test_that("invalid arguments are rejected", {
  expect_error(nb_support_max("a", 5), "single numeric")
  expect_error(nb_support_max(1, c(5, 6)), "single numeric")
  expect_error(nb_support_max(1, 5, tol = 0), "positive")
})
