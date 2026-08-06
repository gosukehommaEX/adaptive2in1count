test_that("the region contains only zeros and ones of the right shape", {
  for (smax in c(0, 1, 40, 80)) {
    region <- exact_conditional_region(40, 40, smax, 0.025)
    expect_true(all(region %in% c(0, 1)))
    expect_equal(dim(region), c(smax + 1L, smax + 1L))
  }
})

test_that("the test size never exceeds the nominal level at any rate", {
  settings <- expand.grid(size = c(20, 40, 80), dispersion = c(0.4, 1, 2),
                          rate = c(0.5, 2))
  for (i in seq_len(nrow(settings))) {
    shape <- settings$size[i] * settings$dispersion[i]
    mean <- settings$size[i] * settings$rate[i]
    smax <- nb_support_max(shape, mean)
    region <- exact_conditional_region(shape, shape, smax, 0.025)
    pmf <- nb_count_pmf(shape, mean, smax)
    expect_lte(sum(outer(pmf, pmf) * region), 0.025)
  }
})

test_that("the test size never exceeds the level under unequal allocation", {
  dispersion <- 0.8
  rate <- 1.5
  for (sizes in list(c(30, 15), c(20, 40), c(45, 15))) {
    shape <- sizes * dispersion
    smax <- max(nb_support_max(shape[1], sizes[1] * rate),
                nb_support_max(shape[2], sizes[2] * rate))
    region <- exact_conditional_region(shape[1], shape[2], smax, 0.025)
    pmf1 <- nb_count_pmf(shape[1], sizes[1] * rate, smax)
    pmf2 <- nb_count_pmf(shape[2], sizes[2] * rate, smax)
    expect_lte(sum(outer(pmf1, pmf2) * region), 0.025)
  }
})

test_that("the size never exceeds the level at other significance levels", {
  shape <- 30 * 0.8
  smax <- nb_support_max(shape, 30 * 1.5)
  pmf <- nb_count_pmf(shape, 30 * 1.5, smax)
  for (alpha in c(0.005, 0.01, 0.05, 0.1)) {
    region <- exact_conditional_region(shape, shape, smax, alpha)
    expect_lte(sum(outer(pmf, pmf) * region), alpha)
  }
})

test_that("the region grows with the significance level", {
  levels <- c(0.005, 0.01, 0.025, 0.05, 0.1)
  sizes <- vapply(levels, function(a) sum(exact_conditional_region(30, 30, 60, a)),
                  numeric(1))
  expect_true(all(diff(sizes) >= 0))
  for (i in seq_len(length(levels) - 1)) {
    expect_true(all(exact_conditional_region(30, 30, 60, levels[i]) <=
                      exact_conditional_region(30, 30, 60, levels[i + 1])))
  }
})

test_that("the region points in the direction of a lower investigational rate", {
  region <- exact_conditional_region(40, 40, 60, 0.025)
  expect_equal(region[1, 61], 1)
  expect_equal(region[61, 1], 0)
  expect_equal(sum(diag(region)), 0)
})

test_that("no rejection is possible when the total is too small", {
  region <- exact_conditional_region(40, 40, 60, 0.025)
  expect_equal(sum(region[seq(1, 4), seq(1, 4)]), 0)
})

test_that("invalid arguments are rejected", {
  expect_error(exact_conditional_region(0, 40, 60, 0.025), "positive")
  expect_error(exact_conditional_region(40, -1, 60, 0.025), "positive")
  expect_error(exact_conditional_region(40, 40, -1, 0.025), "non-negative")
  expect_error(exact_conditional_region(40, 40, 60, 1.5), "between 0 and 1")
  expect_error(exact_conditional_region(40, 40, 60, 0), "between 0 and 1")
})
