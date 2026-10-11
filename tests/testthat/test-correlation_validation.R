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

test_that("the correlations of the three statistics match the closed forms", {
  set.seed(20260808)
  # The tolerances below are set against the Monte Carlo error of a sample
  # correlation at this number of replications, (1 - rho^2) / sqrt(4000),
  # which is about 0.005 for the interim and small trial pair when both read
  # the Phase 2 endpoint and about 0.015 for the other pairs.
  replications <- 4000
  sizes <- c(40, 60, 150)
  # The closed forms apply the delta method to the log of a ratio of totals and
  # treat the total exposure as fixed at its expectation, which is what the
  # negative binomial model does when it conditions on the offset. Under
  # complete follow-up the total exposure is fixed, so the closed forms are
  # exact up to the delta method itself and a tight tolerance applies. When the
  # follow-up period is random the total exposure varies and is itself
  # correlated with the counts and across the two endpoints, so the closed forms
  # carry an additional approximation error. Its size is checked here rather
  # than assumed to be negligible.
  settings <- list(
    list(dropout = 0, shape = "complete", endpoint = "phase2",
         tolerance = 0.035),
    list(dropout = 0.4, shape = "exponential", endpoint = "phase2",
         tolerance = 0.08),
    list(dropout = 0.4, shape = "exponential", endpoint = "phase3",
         tolerance = 0.08),
    list(dropout = 0.4, shape = "uniform", endpoint = "phase2",
         tolerance = 0.08))
  discrepancy <- numeric(length(settings))
  for (i in seq_along(settings)) {
    s <- settings[[i]]
    # The three analyses read the same patients, the interim taking the first
    # block and the later analyses adding to it. Drawing a fresh data set for
    # each analysis would destroy the very dependence that is being checked.
    groups <- lapply(seq_len(2), function(g) {
      patients <- draw_patients(replications * sizes[3], 2, 1.3, 0.7, 0.8,
                                0.35, 0.5, 1, s$dropout, s$shape)
      lapply(patients, function(z) matrix(z, nrow = replications))
    })
    contrast <- function(size, endpoint) {
      count <- paste0("count", endpoint)
      exposure <- paste0("exposure", endpoint)
      rate <- vapply(seq_len(2), function(g) {
        columns <- seq_len(size)
        rowSums(groups[[g]][[count]][, columns, drop = FALSE]) /
          rowSums(groups[[g]][[exposure]][, columns, drop = FALSE])
      }, numeric(replications))
      log(rate[, 1]) - log(rate[, 2])
    }
    small_endpoint <- if (s$endpoint == "phase2") 1L else 2L
    statistic <- cbind(contrast(sizes[1], 1L),
                       contrast(sizes[2], small_endpoint),
                       contrast(sizes[3], 2L))
    usable <- apply(is.finite(statistic), 1, all)
    observed <- stats::cor(statistic[usable, ])
    expected <- design_components(2, 1.3, 0.7, 0.8, 0.35, sizes[1], sizes[2],
                                  sizes[3], 0.5, 1, s$dropout, s$shape,
                                  small_endpoint = s$endpoint)$correlation
    discrepancy[i] <- max(abs(observed[1, 2] - expected[["interim_small"]]),
                          abs(observed[1, 3] - expected[["interim_large"]]))
    expect_lt(abs(observed[1, 2] - expected[["interim_small"]]), s$tolerance)
    expect_lt(abs(observed[1, 3] - expected[["interim_large"]]), s$tolerance)
    # The ordering that carries type I error control must survive in the
    # simulated statistics as well as in the closed forms.
    expect_gte(expected[["interim_small"]], expected[["interim_large"]])
    expect_gt(observed[1, 2] - observed[1, 3], -s$tolerance)
  }
  # A random follow-up period is the only source of the extra approximation
  # error, so the first setting, which has complete follow-up, must agree
  # better than the worst of the settings that have dropout.
  expect_lt(discrepancy[1], max(discrepancy[-1]))
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
