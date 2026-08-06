small_design <- function(...) {
  adaptive2in1count_exact(rate1 = 2, rate2 = 0.5, dispersion1 = 0.6,
                    dispersion2 = 1.2, shared = 0.4,
                    n1 = 20, n2 = 30, n3 = 60, ...)
}

test_that("the rejection probability is the sum of its two parts", {
  fit <- small_design(cutpoint = c(0.4, 0.7, 1, 2))
  expect_equal(fit$results$rejection,
               fit$results$rejection_small + fit$results$rejection_large,
               tolerance = 1e-12)
})

test_that("the probability of expansion increases with the cutpoint", {
  fit <- small_design(cutpoint = c(0.2, 0.3, 0.5, 0.8, 1.2, 2, 5))
  expect_true(all(diff(fit$results$go_probability) >= 0))
  expect_true(all(diff(fit$results$expected_size) >= 0))
})

test_that("the expected sample size lies between the two trial sizes", {
  fit <- small_design(cutpoint = c(0.4, 0.8, 1.5))
  expect_true(all(fit$results$expected_size >= 60))
  expect_true(all(fit$results$expected_size <= 120))
})

test_that("evaluating a vector of cutpoints matches one call per cutpoint", {
  cutpoint <- c(0.3, 0.6, 1.1)
  joint <- small_design(cutpoint = cutpoint)
  for (i in seq_along(cutpoint)) {
    single <- small_design(cutpoint = cutpoint[i])
    expect_equal(single$results$rejection, joint$results$rejection[i],
                 tolerance = 1e-12)
    expect_equal(single$results$go_probability, joint$results$go_probability[i],
                 tolerance = 1e-12)
  }
})

test_that("the type I error rate never exceeds the nominal level", {
  grid <- expand.grid(dispersion1 = c(0.4, 1), dispersion2 = c(0.8, 2),
                      rate1 = c(1, 3), rate2 = c(0.3, 0.8))
  for (i in seq_len(nrow(grid))) {
    fit <- adaptive2in1count_exact(rate1 = grid$rate1[i], rate2 = grid$rate2[i],
                             dispersion1 = grid$dispersion1[i],
                             dispersion2 = grid$dispersion2[i],
                             shared = min(grid$dispersion1[i],
                                          grid$dispersion2[i]) / 2,
                             n1 = 15, n2 = 25, n3 = 45,
                             cutpoint = c(0.4, 0.7, 1, 1.6, 4))
    expect_true(all(fit$results$rejection <= 0.025))
  }
})

test_that("the type I error rate of the Wald version stays near the level", {
  grid <- expand.grid(dispersion1 = c(0.4, 1), dispersion2 = c(0.8, 2))
  for (i in seq_len(nrow(grid))) {
    fit <- adaptive2in1count_exact(rate1 = 2, rate2 = 0.6,
                             dispersion1 = grid$dispersion1[i],
                             dispersion2 = grid$dispersion2[i],
                             shared = min(grid$dispersion1[i],
                                          grid$dispersion2[i]) / 2,
                             n1 = 15, n2 = 25, n3 = 45,
                             cutpoint = c(0.4, 1, 4), test = "wald")
    expect_true(all(fit$results$rejection <= 0.035))
  }
})

test_that("a cutpoint that never expands reproduces the fixed sample small trial", {
  smax <- nb_support_max(30 * 0.6, 30 * 2)
  region <- exact_conditional_region(30 * 0.6, 30 * 0.6, smax, 0.025)
  pmf <- nb_count_pmf(30 * 0.6, 30 * 2, smax)
  fit <- small_design(cutpoint = 1e-6)
  expect_equal(fit$results$go_probability, 0, tolerance = 1e-12)
  expect_equal(fit$results$rejection_large, 0, tolerance = 1e-12)
  expect_equal(fit$results$rejection, sum(outer(pmf, pmf) * region),
               tolerance = 1e-6)
  expect_equal(fit$results$expected_size, 60, tolerance = 1e-8)
})

test_that("a cutpoint that always expands reproduces the fixed sample large trial", {
  tmax <- nb_support_max(60 * 1.2, 60 * 0.5)
  region <- exact_conditional_region(60 * 1.2, 60 * 1.2, tmax, 0.025)
  pmf <- nb_count_pmf(60 * 1.2, 60 * 0.5, tmax)
  fit <- small_design(cutpoint = 1e6)
  expect_equal(fit$results$go_probability, 1, tolerance = 1e-8)
  expect_equal(fit$results$rejection_small, 0, tolerance = 1e-12)
  expect_equal(fit$results$rejection, sum(outer(pmf, pmf) * region),
               tolerance = 1e-6)
  expect_equal(fit$results$expected_size, 120, tolerance = 1e-6)
})

test_that("a stronger treatment effect gives a larger rejection probability", {
  null <- small_design(cutpoint = 0.6)
  medium <- small_design(cutpoint = 0.6, ratio1 = 0.6, ratio2 = 0.8)
  strong <- small_design(cutpoint = 0.6, ratio1 = 0.4, ratio2 = 0.6)
  expect_lt(null$results$rejection, medium$results$rejection)
  expect_lt(medium$results$rejection, strong$results$rejection)
  expect_lt(null$results$go_probability, medium$results$go_probability)
  expect_lt(medium$results$go_probability, strong$results$go_probability)
})

test_that("the small trial part does not depend on the shared frailty shape", {
  values <- vapply(c(0, 0.2, 0.4, 0.6), function(s) {
    adaptive2in1count_exact(rate1 = 2, rate2 = 0.5, dispersion1 = 0.6,
                      dispersion2 = 1.2, shared = s,
                      n1 = 20, n2 = 30, n3 = 60,
                      cutpoint = 0.6)$results$rejection_small
  }, numeric(1))
  expect_equal(diff(range(values)), 0, tolerance = 1e-12)
})

test_that("the large trial part grows with the shared frailty shape", {
  values <- vapply(c(0, 0.2, 0.4, 0.6), function(s) {
    adaptive2in1count_exact(rate1 = 2, rate2 = 0.5, dispersion1 = 0.6,
                      dispersion2 = 1.2, shared = s,
                      n1 = 20, n2 = 30, n3 = 60,
                      cutpoint = 0.6)$results$rejection_large
  }, numeric(1))
  expect_true(all(diff(values) > 0))
  expect_equal(values, c(0.002665341, 0.003357119, 0.004141920, 0.005022382),
               tolerance = 1e-6)
})

test_that("the probability of expansion does not depend on the Phase 3 endpoint", {
  first <- small_design(cutpoint = 0.6)
  second <- adaptive2in1count_exact(rate1 = 2, rate2 = 1.4, dispersion1 = 0.6,
                              dispersion2 = 0.9, shared = 0.5,
                              n1 = 20, n2 = 30, n3 = 60, cutpoint = 0.6)
  expect_equal(first$results$go_probability, second$results$go_probability,
               tolerance = 1e-10)
})

test_that("the reference multiple sclerosis design reproduces its stored values", {
  fit <- adaptive2in1count_exact(rate1 = 3, rate2 = 0.6, dispersion1 = 0.5,
                           dispersion2 = 1.5, shared = 0.25,
                           n1 = 50, n2 = 80, n3 = 250,
                           cutpoint = c(0.5, 0.8), alpha = 0.025)
  expect_equal(fit$results$go_probability, c(0.012877470, 0.233621003),
               tolerance = 1e-6)
  expect_equal(fit$results$rejection, c(0.017735431, 0.007297679),
               tolerance = 1e-6)
  expect_equal(fit$results$rejection_small, c(0.017324112, 0.001016789),
               tolerance = 1e-6)
  expect_equal(fit$results$rejection_large, c(0.000411319, 0.006280891),
               tolerance = 1e-6)
  expect_equal(fit$results$expected_size, c(164.378340, 239.431141),
               tolerance = 1e-6)
})

test_that("the reference design under the alternative is reproduced", {
  fit <- adaptive2in1count_exact(rate1 = 3, rate2 = 0.6, dispersion1 = 0.5,
                           dispersion2 = 1.5, shared = 0.25,
                           n1 = 50, n2 = 80, n3 = 250, cutpoint = 0.5,
                           ratio1 = 0.4, ratio2 = 0.7)
  expect_equal(fit$results$go_probability, 0.754010248, tolerance = 1e-6)
  expect_equal(fit$results$rejection, 0.707337474, tolerance = 1e-6)
  expect_equal(fit$results$rejection_small, 0.198354411, tolerance = 1e-6)
  expect_equal(fit$results$rejection_large, 0.508983063, tolerance = 1e-6)
})

test_that("the small reference designs reproduce their stored values", {
  fit <- small_design(cutpoint = c(0.5, 0.6, 1))
  expect_equal(fit$results$go_probability,
               c(0.072199759, 0.140854650, 0.511067872), tolerance = 1e-6)
  expect_equal(fit$results$rejection,
               c(0.007143308, 0.005909000, 0.012018324), tolerance = 1e-6)
  expect_equal(fit$results$rejection_small,
               c(0.004837009, 0.001767080, 0.000039033), tolerance = 1e-5)
  expect_equal(fit$results$rejection_large,
               c(0.002306299, 0.004141920, 0.011979291), tolerance = 1e-6)
})

test_that("the small reference design under the alternative reproduces its values", {
  fit <- small_design(cutpoint = 0.6, ratio1 = 0.4, ratio2 = 0.6)
  expect_equal(fit$results$go_probability, 0.788817945, tolerance = 1e-6)
  expect_equal(fit$results$rejection, 0.247392234, tolerance = 1e-6)
  expect_equal(fit$results$rejection_small, 0.019248050, tolerance = 1e-6)
  expect_equal(fit$results$rejection_large, 0.228144184, tolerance = 1e-6)
})

test_that("the Wald version reproduces its stored values", {
  fit <- small_design(cutpoint = 0.6, test = "wald")
  expect_equal(fit$results$rejection, 0.007298578, tolerance = 1e-6)
  expect_equal(fit$results$rejection_small, 0.002252748, tolerance = 1e-6)
  expect_equal(fit$results$rejection_large, 0.005045831, tolerance = 1e-6)
})

test_that("the extreme shared frailty shapes reproduce their stored values", {
  none <- adaptive2in1count_exact(rate1 = 2, rate2 = 0.5, dispersion1 = 0.6,
                            dispersion2 = 1.2, shared = 0,
                            n1 = 20, n2 = 30, n3 = 60, cutpoint = 0.6)
  full <- adaptive2in1count_exact(rate1 = 2, rate2 = 0.5, dispersion1 = 0.6,
                            dispersion2 = 1.2, shared = 0.6,
                            n1 = 20, n2 = 30, n3 = 60, cutpoint = 0.6)
  expect_equal(none$results$rejection, 0.004432421, tolerance = 1e-6)
  expect_equal(full$results$rejection, 0.006789462, tolerance = 1e-6)
})

test_that("designs without added patients reproduce their stored values", {
  same_interim <- adaptive2in1count_exact(rate1 = 2, rate2 = 0.5, dispersion1 = 0.6,
                                    dispersion2 = 1.2, shared = 0.4,
                                    n1 = 20, n2 = 20, n3 = 60, cutpoint = 0.6)
  same_final <- adaptive2in1count_exact(rate1 = 2, rate2 = 0.5, dispersion1 = 0.6,
                                  dispersion2 = 1.2, shared = 0.4,
                                  n1 = 20, n2 = 40, n3 = 40, cutpoint = 0.6)
  expect_equal(same_interim$results$rejection, 0.004141920, tolerance = 1e-6)
  expect_lt(same_interim$results$rejection_small, 1e-9)
  expect_equal(same_final$results$rejection, 0.008997179, tolerance = 1e-6)
  expect_equal(same_final$results$expected_size, 80, tolerance = 1e-8)
})

test_that("unequal allocation is accepted and controls the type I error rate", {
  for (allocation in list(c(24, 16), c(16, 24), c(30, 10))) {
    fit <- adaptive2in1count_exact(rate1 = 2, rate2 = 0.5, dispersion1 = 0.6,
                             dispersion2 = 1.2, shared = 0.4,
                             n1 = allocation, n2 = allocation * 2,
                             n3 = allocation * 3,
                             cutpoint = c(0.5, 1, 2))
    expect_equal(fit$n1, allocation)
    expect_true(all(fit$results$rejection <= 0.025))
  }
})

test_that("the returned object carries its inputs and class", {
  fit <- small_design(cutpoint = 0.6, alpha = 0.05, ratio1 = 0.5)
  expect_s3_class(fit, "adaptive2in1count_exact")
  expect_equal(fit$alpha, 0.05)
  expect_equal(fit$ratio1, 0.5)
  expect_equal(fit$test, "exact")
  expect_named(fit$support, c("interim", "small", "large"))
  expect_equal(nrow(fit$results), 1L)
})

test_that("the small trial may analyse the Phase 3 endpoint", {
  fit <- small_design(cutpoint = c(0.4, 0.7, 1.5), small_endpoint = "phase3")
  expect_equal(fit$small_endpoint, "phase3")
  expect_true(all(fit$results$rejection <= 0.025))
  expect_equal(fit$results$rejection,
               fit$results$rejection_small + fit$results$rejection_large,
               tolerance = 1e-12)
})

test_that("a Phase 3 small trial reproduces the fixed sample trial at both ends", {
  dispersion2 <- 1.2
  rate2 <- 0.5
  for (side in c("small", "large")) {
    size <- if (side == "small") 30 else 60
    ymax <- nb_support_max(size * dispersion2, size * rate2)
    region <- exact_conditional_region(size * dispersion2, size * dispersion2,
                                       ymax, 0.025)
    pmf <- nb_count_pmf(size * dispersion2, size * rate2, ymax)
    reference <- sum(outer(pmf, pmf) * region)
    cut <- if (side == "small") 1e-6 else 1e6
    fit <- small_design(cutpoint = cut, small_endpoint = "phase3")
    expect_equal(fit$results$rejection, reference, tolerance = 1e-6)
  }
})

test_that("both small trial endpoints agree when the trial never expands", {
  common <- list(rate1 = 1.5, rate2 = 1.5, dispersion1 = 0.9,
                 dispersion2 = 0.9, shared = 0.9, n1 = 20, n2 = 30, n3 = 60,
                 cutpoint = 1e-6)
  first <- do.call(adaptive2in1count, c(common, list(small_endpoint = "phase2")))
  second <- do.call(adaptive2in1count, c(common, list(small_endpoint = "phase3")))
  expect_equal(first$results$go_probability, 0, tolerance = 1e-12)
  expect_equal(first$results$rejection, second$results$rejection,
               tolerance = 1e-6)
})

test_that("the small trial endpoint changes how the interim decision links to it", {
  common <- list(rate1 = 1.5, rate2 = 1.5, dispersion1 = 0.9,
                 dispersion2 = 0.9, shared = 0.9, n1 = 20, n2 = 30, n3 = 60,
                 cutpoint = 0.7)
  first <- do.call(adaptive2in1count, c(common, list(small_endpoint = "phase2")))
  second <- do.call(adaptive2in1count, c(common, list(small_endpoint = "phase3")))
  expect_equal(first$results$go_probability, second$results$go_probability,
               tolerance = 1e-12)
  expect_false(isTRUE(all.equal(first$results$rejection_small,
                                second$results$rejection_small)))
  expect_lte(first$results$rejection, 0.025)
  expect_lte(second$results$rejection, 0.025)
})

test_that("invalid arguments are rejected", {
  expect_error(adaptive2in1count_exact(-1, 0.5, 0.6, 1.2, 0.4, 20, 30, 60, 0.6),
               "positive")
  expect_error(adaptive2in1count_exact(2, 0.5, 0.6, 1.2, 2, 20, 30, 60, 0.6),
               "smaller of the two dispersions")
  expect_error(adaptive2in1count_exact(2, 0.5, 0.6, 1.2, 0.4, 40, 30, 60, 0.6),
               "non-decreasing")
  expect_error(adaptive2in1count_exact(2, 0.5, 0.6, 1.2, 0.4, 20, 30, 60, -1),
               "positive values")
  expect_error(adaptive2in1count_exact(2, 0.5, 0.6, 1.2, 0.4, 20, 30, 60,
                                       numeric(0)), "positive values")
  expect_error(adaptive2in1count_exact(2, 0.5, 0.6, 1.2, 0.4, 20, 30, 60, 0.6,
                                       alpha = 1.2), "between 0 and 1")
  expect_error(adaptive2in1count_exact(2, 0.5, 0.6, 1.2, 0.4, 20, 30, 60, 0.6,
                                       test = "score"), "should be one of")
  expect_error(adaptive2in1count_exact(2, 0.5, 0.6, 1.2, 0.4, 20, 30, 60, 0.6,
                                       small_endpoint = "phase1"),
               "should be one of")
})
