# The closed forms for the correlations are checked against data generated
# directly from the model. The generator below repeats the model definition
# rather than calling the package, so that a sign or a scale error in the
# closed forms cannot be masked.

draw_follow_up <- function(size, window, dropout, shape) {
  if (shape == "complete" || dropout == 0) {
    return(rep(window, size))
  }
  if (shape == "uniform") {
    early <- stats::runif(size) < dropout
    return(ifelse(early, stats::runif(size, 0, window), window))
  }
  pmin(stats::rexp(size, rate = -log(1 - dropout) / window), window)
}

draw_gamma <- function(size, shape_parameter) {
  if (shape_parameter <= 0) {
    return(rep(0, size))
  }
  stats::rgamma(size, shape = shape_parameter, rate = 1)
}

draw_patients <- function(size, rate1, rate2, dispersion1, dispersion2, shared,
                          window1, window2, dropout, shape) {
  follow_up <- draw_follow_up(size, window2, dropout, shape)
  exposure1 <- pmin(follow_up, window1)
  base <- draw_gamma(size, shared)
  count1 <- stats::rpois(size, rate1 * exposure1 *
                           (base + draw_gamma(size, dispersion1 - shared)) /
                           dispersion1)
  count2 <- stats::rpois(size, rate2 * follow_up *
                           (base + draw_gamma(size, dispersion2 - shared)) /
                           dispersion2)
  list(count1 = count1, count2 = count2, exposure1 = exposure1,
       exposure2 = follow_up)
}

test_that("the within-patient endpoint correlation matches the closed form", {
  set.seed(20260807)
  settings <- list(
    list(shared = 0.7, dropout = 0, shape = "complete"),
    list(shared = 0.35, dropout = 0, shape = "complete"),
    list(shared = 0.35, dropout = 0.4, shape = "exponential"),
    list(shared = 0.7, dropout = 0.4, shape = "uniform"))
  discrepancy <- numeric(length(settings))
  for (i in seq_along(settings)) {
    s <- settings[[i]]
    patients <- draw_patients(2e5, 2, 1.3, 0.7, 0.8, s$shared, 0.5, 1,
                              s$dropout, s$shape)
    observed <- stats::cor(patients$count1, patients$count2)
    expected <- endpoint_correlation(2, 1.3, 0.7, 0.8, s$shared, 0.5, 1,
                                     s$dropout, s$shape)
    expect_lt(abs(observed - expected), 0.01)
  }
})

test_that("dropout alone correlates the two endpoints", {
  set.seed(11)
  for (shape in c("uniform", "exponential")) {
    patients <- draw_patients(3e5, 2, 1.3, 0.7, 0.8, shared = 0, 0.5, 1,
                              dropout = 0.5, shape = shape)
    observed <- stats::cor(patients$count1, patients$count2)
    expected <- endpoint_correlation(2, 1.3, 0.7, 0.8, 0, 0.5, 1, 0.5, shape)
    expect_gt(observed, 0.01)
    expect_gt(expected, 0.01)
    expect_lt(abs(observed - expected), 0.01)
  }
})

test_that("complete follow-up without a shared frailty leaves no correlation", {
  set.seed(12)
  patients <- draw_patients(2e5, 2, 1.3, 0.7, 0.8, shared = 0, 0.5, 1,
                            dropout = 0, shape = "complete")
  expect_lt(abs(stats::cor(patients$count1, patients$count2)), 0.01)
})

test_that("the correlation of the statistics matches the covariance factor", {
  set.seed(20261012)
  # Each statistic is the log rate ratio of the negative binomial model with an
  # offset for the exposure, which is the analysis the closed forms describe,
  # with the dispersion held at its design value. Both statistics read the same
  # patients, so their correlation is the correlation of the two estimated log
  # rate ratios. The settings with dropout and no shared frailty, and with a
  # high event rate, separate the covariance factor from one that would let the
  # covariance of the two exposures through.
  replications <- 4000
  size <- 50
  settings <- list(
    list(rate1 = 2, rate2 = 1.3, dispersion1 = 0.7, dispersion2 = 0.8,
         shared = 0.4, window1 = 0.5, dropout = 0, shape = "complete"),
    list(rate1 = 2, rate2 = 1.3, dispersion1 = 0.7, dispersion2 = 0.8,
         shared = 0.4, window1 = 0.5, dropout = 0.4, shape = "exponential"),
    list(rate1 = 2, rate2 = 1.3, dispersion1 = 0.7, dispersion2 = 0.8,
         shared = 0, window1 = 1, dropout = 0.6, shape = "exponential"),
    list(rate1 = 10, rate2 = 10, dispersion1 = 1, dispersion2 = 1,
         shared = 1, window1 = 1, dropout = 0.4, shape = "uniform"))
  for (s in settings) {
    groups <- lapply(seq_len(2), function(g) {
      patients <- draw_patients(replications * size, s$rate1, s$rate2,
                                s$dispersion1, s$dispersion2, s$shared,
                                s$window1, 1, s$dropout, s$shape)
      lapply(patients, function(z) matrix(z, nrow = replications))
    })
    estimate <- function(endpoint, dispersion) {
      count <- paste0("count", endpoint)
      exposure <- paste0("exposure", endpoint)
      negative_binomial_contrast(groups[[1]][[count]], groups[[1]][[exposure]],
                                 groups[[2]][[count]], groups[[2]][[exposure]],
                                 dispersion = dispersion)$estimate
    }
    statistic <- cbind(estimate(1L, s$dispersion1),
                       estimate(2L, s$dispersion2))
    usable <- apply(is.finite(statistic), 1, all)
    observed <- stats::cor(statistic[usable, ])[1, 2]
    expected <- design_components(
      s$rate1, s$rate2, s$dispersion1, s$dispersion2, s$shared, size, size,
      size, s$window1, 1, s$dropout, s$shape)$correlation[["interim_large"]]
    # Monte Carlo error of a sample correlation at this number of replications.
    error <- (1 - expected^2) / sqrt(sum(usable))
    expect_lt(abs(observed - expected), 4 * error)
  }
})

test_that("statistics from disjoint patients are uncorrelated", {
  # A control for the test above: if the analyses did not share patients the
  # correlation would vanish, so the agreement found there is driven by the
  # nesting and not by anything else in the generator.
  set.seed(31)
  replications <- 4000
  size <- 50
  draw_statistic <- function() {
    rate <- vapply(seq_len(2), function(g) {
      patients <- draw_patients(replications * size, 2, 1.3, 0.7, 0.8, 0.35,
                                0.5, 1, 0, "complete")
      count <- matrix(patients$count1, nrow = replications)
      exposure <- matrix(patients$exposure1, nrow = replications)
      rowSums(count) / rowSums(exposure)
    }, numeric(replications))
    log(rate[, 1]) - log(rate[, 2])
  }
  first <- draw_statistic()
  second <- draw_statistic()
  expect_lt(abs(stats::cor(first, second)), 0.04)
})
