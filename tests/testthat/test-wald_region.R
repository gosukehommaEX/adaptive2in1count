test_that("the region contains only zeros and ones of the right shape", {
  for (smax in c(1, 30, 60)) {
    region <- wald_region(80, 80, 0.5, smax, 0.025)
    expect_true(all(region %in% c(0, 1)))
    expect_equal(dim(region), c(smax + 1L, smax + 1L))
  }
})

test_that("the region points in the direction of a lower investigational rate", {
  region <- wald_region(80, 80, 0.5, 60, 0.025)
  expect_equal(region[1, 61], 1)
  expect_equal(region[61, 1], 0)
  expect_equal(sum(diag(region)), 0)
})

test_that("the region grows with the significance level", {
  levels <- c(0.005, 0.01, 0.025, 0.05, 0.1)
  for (i in seq_len(length(levels) - 1)) {
    expect_true(all(wald_region(80, 80, 0.5, 60, levels[i]) <=
                      wald_region(80, 80, 0.5, 60, levels[i + 1])))
  }
})

test_that("the size of the Wald test stays near the nominal level", {
  settings <- expand.grid(size = c(20, 40, 80), dispersion = c(0.4, 1, 2),
                          rate = c(0.5, 2))
  for (i in seq_len(nrow(settings))) {
    size <- settings$size[i]
    dispersion <- settings$dispersion[i]
    rate <- settings$rate[i]
    smax <- nb_support_max(size * dispersion, size * rate)
    region <- wald_region(size, size, dispersion, smax, 0.025)
    pmf <- nb_count_pmf(size * dispersion, size * rate, smax)
    actual <- sum(outer(pmf, pmf) * region)
    expect_gt(actual, 0.01)
    expect_lt(actual, 0.04)
  }
})

test_that("the Wald test is larger than the exact conditional test", {
  settings <- expand.grid(size = c(20, 40, 80), dispersion = c(0.4, 1, 2),
                          rate = c(0.5, 2))
  for (i in seq_len(nrow(settings))) {
    size <- settings$size[i]
    dispersion <- settings$dispersion[i]
    rate <- settings$rate[i]
    smax <- nb_support_max(size * dispersion, size * rate)
    pmf <- nb_count_pmf(size * dispersion, size * rate, smax)
    joint <- outer(pmf, pmf)
    wald <- sum(joint * wald_region(size, size, dispersion, smax, 0.025))
    exact <- sum(joint * exact_conditional_region(size * dispersion,
                                                  size * dispersion, smax,
                                                  0.025))
    expect_gt(wald, exact)
  }
})

test_that("unequal allocation is accepted", {
  region <- wald_region(60, 30, 0.8, 80, 0.025)
  expect_equal(dim(region), c(81L, 81L))
  expect_true(all(region %in% c(0, 1)))
})

test_that("invalid arguments are rejected", {
  expect_error(wald_region(-1, 80, 0.5, 60, 0.025), "positive")
  expect_error(wald_region(80, 0, 0.5, 60, 0.025), "positive")
  expect_error(wald_region(80, 80, 0, 60, 0.025), "positive")
  expect_error(wald_region(80, 80, 0.5, -2, 0.025), "non-negative")
  expect_error(wald_region(80, 80, 0.5, 60, 0), "between 0 and 1")
})
