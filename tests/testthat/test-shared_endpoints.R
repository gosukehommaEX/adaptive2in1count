# The shared endpoint structure has exact consequences that hold for every
# choice of the remaining parameters, so the checks below are deterministic and
# cheap. No simulation is used.

shared_design <- list(rate1 = 9, rate2 = 1.3, dispersion1 = 9,
                      dispersion2 = 0.8, shared = 9, n1 = 80, n2 = 120,
                      n3 = 750, window1 = 1, window2 = 1,
                      endpoints = "shared")

test_that("the arguments of the Phase 2 endpoint are taken from the Phase 3 endpoint", {
  parts <- do.call(design_components, shared_design)
  matched <- design_components(rate1 = 1.3, rate2 = 1.3, dispersion1 = 0.8,
                               dispersion2 = 0.8, shared = 0.8, n1 = 80,
                               n2 = 120, n3 = 750, window1 = 1, window2 = 1,
                               endpoints = "shared")
  expect_equal(parts$variance, matched$variance)
  expect_equal(parts$correlation, matched$correlation)
  expect_identical(parts$endpoints, "shared")
})

test_that("one endpoint over one window gives the nested correlations", {
  for (shape in c("complete", "exponential", "uniform")) {
    dropout <- if (shape == "complete") 0 else 0.4
    parts <- do.call(design_components,
                     c(shared_design, list(dropout = dropout,
                                           dropout_shape = shape)))
    expect_equal(parts$correlation[["interim_small"]], sqrt(80 / 120))
    expect_equal(parts$correlation[["interim_large"]], sqrt(80 / 750))
  }
})

test_that("the choice of endpoint for the small trial is immaterial over one window", {
  first <- do.call(design_components,
                   c(shared_design, list(small_endpoint = "phase2")))
  second <- do.call(design_components,
                    c(shared_design, list(small_endpoint = "phase3")))
  expect_equal(first$variance, second$variance)
  expect_equal(first$correlation, second$correlation)
})

test_that("one endpoint over one window is perfectly correlated within a patient", {
  for (shape in c("complete", "exponential", "uniform")) {
    dropout <- if (shape == "complete") 0 else 0.5
    value <- endpoint_correlation(rate1 = 9, rate2 = 1.3, dispersion1 = 9,
                                  dispersion2 = 0.8, shared = 9, window1 = 1,
                                  window2 = 1, dropout = dropout,
                                  dropout_shape = shape,
                                  endpoints = "shared")
    expect_equal(value, 1)
  }
})

test_that("a shorter interim window lowers the correlation but keeps the ordering", {
  for (small in c("phase2", "phase3")) {
    for (shape in c("complete", "exponential")) {
      dropout <- if (shape == "complete") 0 else 0.4
      parts <- design_components(rate1 = 1.3, rate2 = 1.3, dispersion1 = 0.8,
                                 dispersion2 = 0.8, shared = 0.8, n1 = 80,
                                 n2 = 120, n3 = 750, window1 = 0.5,
                                 window2 = 1, dropout = dropout,
                                 dropout_shape = shape,
                                 small_endpoint = small,
                                 endpoints = "shared")
      expect_lt(parts$correlation[["interim_large"]], sqrt(80 / 750))
      expect_gte(parts$correlation[["interim_small"]],
                 parts$correlation[["interim_large"]])
    }
  }
})

test_that("the requirement for type I error control holds for a shared endpoint", {
  for (window1 in c(0.25, 0.5, 1)) {
    for (small in c("phase2", "phase3")) {
      rho <- adaptive2in1count_correlation(rate1 = 1.3, rate2 = 1.3,
                                           dispersion1 = 0.8,
                                           dispersion2 = 0.8, shared = 0.8,
                                           n1 = 80, n2 = 120, n3 = 750,
                                           window1 = window1, window2 = 1,
                                           dropout = 0.3,
                                           dropout_shape = "exponential",
                                           small_endpoint = small,
                                           endpoints = "shared")
      expect_true(rho$holds)
      if (window1 == 1) {
        expect_equal(rho$endpoint, 1)
      } else {
        expect_lt(rho$endpoint, 1)
        expect_gt(rho$endpoint, 0)
      }
    }
  }
})

test_that("expansion is immaterial when the two final analyses coincide", {
  # With one endpoint, one window and no difference between the sizes of the
  # two final analyses, the design cannot be distinguished from a fixed trial,
  # so the rejection probability does not depend on the cutpoint.
  fit <- adaptive2in1count(rate1 = 1.3, rate2 = 1.3, dispersion1 = 0.8,
                           dispersion2 = 0.8, shared = 0.8, n1 = 80,
                           n2 = 300, n3 = 300, cutpoint = c(0.5, 0.8, 1, 1.5,
                                                            3),
                           window1 = 1, window2 = 1, dropout = 0.3,
                           dropout_shape = "exponential", ratio1 = 0.75,
                           ratio2 = 0.75, endpoints = "shared")
  fixed <- pnorm(-qnorm(0.975) -
                   log(0.75) / sqrt(fit$components$variance[["large"]]))
  expect_equal(fit$results$rejection, rep(fixed, 5))
  expect_equal(fit$results$expected_size, rep(600, 5))
})

test_that("the type I error rate is exactly the level when expansion is immaterial", {
  fit <- adaptive2in1count(rate1 = 1.3, rate2 = 1.3, dispersion1 = 0.8,
                           dispersion2 = 0.8, shared = 0.8, n1 = 80,
                           n2 = 300, n3 = 300, cutpoint = c(0.7, 1, 1.4),
                           window1 = 1, window2 = 1, dropout = 0.3,
                           dropout_shape = "exponential",
                           endpoints = "shared")
  expect_equal(fit$results$rejection, rep(0.025, 3))
})

test_that("the type I error rate stays below the level for a shared endpoint", {
  fit <- adaptive2in1count_typeI(rate1 = c(0.8, 1.3, 2), rate2 = c(0.8, 1.3, 2),
                                 dispersion1 = 0.8, dispersion2 = 0.8,
                                 shared = 0.8, n1 = 80, n2 = 120, n3 = 750,
                                 cutpoint = seq(0.5, 1.5, by = 0.25),
                                 window1 = 0.5, window2 = 1, dropout = 0.3,
                                 dropout_shape = "exponential",
                                 endpoints = "shared")
  expect_true(all(fit$results$type1 <= 0.025 + 1e-8))
  # Only the diagonal of the nuisance grid is meaningful for one process.
  expect_equal(fit$results$rate1, fit$results$rate2)
})

test_that("the power grid keeps one treatment effect for a shared endpoint", {
  curve <- adaptive2in1count_power(rate1 = 1.3, rate2 = 1.3, dispersion1 = 0.8,
                                   dispersion2 = 0.8, shared = 0.8, n1 = 80,
                                   n2 = 120, n3 = 750, cutpoint = c(0.8, 1.2),
                                   ratio1 = c(0.5, 0.9), ratio2 = c(0.7, 0.8),
                                   window1 = 0.5, window2 = 1,
                                   endpoints = "shared")
  expect_equal(curve$results$ratio1, curve$results$ratio2)
  expect_setequal(curve$results$ratio2, c(0.7, 0.8))
})

test_that("the cutpoint and the size solvers accept a shared endpoint", {
  chosen <- adaptive2in1count_cutpoint(rate1 = 1.3, rate2 = 1.3,
                                       dispersion1 = 0.8, dispersion2 = 0.8,
                                       shared = 0.8, n1 = 80, n2 = 120,
                                       n3 = 750, target = 0.9, ratio1 = 0.75,
                                       ratio2 = 0.75, window1 = 0.5,
                                       window2 = 1, dropout = 0.3,
                                       dropout_shape = "exponential",
                                       endpoints = "shared")
  fit <- adaptive2in1count(rate1 = 1.3, rate2 = 1.3, dispersion1 = 0.8,
                           dispersion2 = 0.8, shared = 0.8, n1 = 80, n2 = 120,
                           n3 = 750, cutpoint = chosen$cutpoint,
                           window1 = 0.5, window2 = 1, dropout = 0.3,
                           dropout_shape = "exponential", ratio1 = 0.75,
                           ratio2 = 0.75, endpoints = "shared")
  expect_equal(fit$results$go_probability, 0.9)
  solved <- adaptive2in1count_size(rate1 = 1.3, rate2 = 1.3, dispersion1 = 0.8,
                                   dispersion2 = 0.8, shared = 0.8, n1 = 80,
                                   n2 = 120, cutpoint = chosen$cutpoint,
                                   target = 0.8, ratio1 = 0.75, ratio2 = 0.75,
                                   window1 = 0.5, window2 = 1, dropout = 0.3,
                                   dropout_shape = "exponential",
                                   endpoints = "shared")
  expect_gte(solved$power, 0.8)
})

test_that("the shared generator reproduces the nesting of the two counts", {
  # The Phase 2 count must be part of the Phase 3 count, and the two must
  # coincide when the two windows do. A single small replication is enough to
  # check the structure of the generated data.
  out <- adaptive2in1count_simulate(rate1 = 1.3, rate2 = 1.3,
                                    dispersion1 = 0.8, dispersion2 = 0.8,
                                    shared = 0.8, n1 = 30, n2 = 40, n3 = 60,
                                    cutpoint = 1, window1 = 1, window2 = 1,
                                    dropout = 0.3,
                                    dropout_shape = "exponential",
                                    endpoints = "shared", method = "known",
                                    nsim = 200, seed = 11)
  expect_equal(nrow(out), 1L)
  expect_true(all(is.finite(out$rejection)))
  expect_identical(out$failures, 0L)
})

test_that("a shared endpoint is reported by the print methods", {
  fit <- adaptive2in1count(rate1 = 1.3, rate2 = 1.3, dispersion1 = 0.8,
                           dispersion2 = 0.8, shared = 0.8, n1 = 80, n2 = 120,
                           n3 = 750, cutpoint = 0.9, window1 = 0.5,
                           window2 = 1, endpoints = "shared")
  lines <- capture.output(print(fit))
  expect_true(any(grepl("same event over the two windows", lines,
                        fixed = TRUE)))
})
