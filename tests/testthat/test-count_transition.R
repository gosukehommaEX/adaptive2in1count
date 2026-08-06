settings <- list(
  list(rate1 = 2, rate2 = 0.5, dispersion1 = 0.6, dispersion2 = 1.2,
       shared = 0.4, n1 = 20, n2 = 30, n3 = 60),
  list(rate1 = 3, rate2 = 1, dispersion1 = 0.5, dispersion2 = 0.5,
       shared = 0.5, n1 = 15, n2 = 25, n3 = 40),
  list(rate1 = 1, rate2 = 0.4, dispersion1 = 2, dispersion2 = 1,
       shared = 0, n1 = 10, n2 = 20, n3 = 30),
  list(rate1 = 4, rate2 = 0.8, dispersion1 = 1.5, dispersion2 = 3,
       shared = 0.75, n1 = 12, n2 = 12, n3 = 36),
  list(rate1 = 2.5, rate2 = 0.6, dispersion1 = 0.8, dispersion2 = 1.6,
       shared = 0.2, n1 = 18, n2 = 30, n3 = 30)
)

build_case <- function(s) {
  amax <- nb_support_max(s$n1 * s$dispersion1, s$n1 * s$rate1)
  smax <- nb_support_max(s$n2 * s$dispersion1, s$n2 * s$rate1)
  tmax <- nb_support_max(s$n3 * s$dispersion2, s$n3 * s$rate2)
  out <- count_transition(rate1 = s$rate1, rate2 = s$rate2,
                          dispersion1 = s$dispersion1,
                          dispersion2 = s$dispersion2, shared = s$shared,
                          n1 = s$n1, n2 = s$n2, n3 = s$n3,
                          amax = amax, smax = smax, tmax = tmax)
  c(out, list(amax = amax, smax = smax, tmax = tmax))
}

test_that("each conditional distribution of the Phase 3 total is sub-probability", {
  for (s in settings) {
    out <- build_case(s)
    total <- rowSums(out$transition)
    expect_true(all(out$transition >= 0))
    expect_true(all(total <= 1 + 1e-12))
  }
})

test_that("the retained mass of the conditional distributions is complete", {
  for (s in settings) {
    out <- build_case(s)
    expect_equal(sum(out$interim * rowSums(out$transition)), 1,
                 tolerance = 1e-8)
  }
})

test_that("the conditional distributions are proper at plausible interim counts", {
  for (s in settings) {
    out <- build_case(s)
    keep <- out$interim > 1e-6
    expect_equal(rowSums(out$transition)[keep], rep(1, sum(keep)),
                 tolerance = 1e-5)
  }
})

test_that("the retained mass decreases as the interim count grows", {
  out <- build_case(settings[[2]])
  total <- rowSums(out$transition)
  expect_lte(total[out$amax + 1], total[1])
  expect_equal(total[1], 1, tolerance = 1e-9)
})

test_that("the interim total is negative binomial", {
  for (s in settings) {
    out <- build_case(s)
    expect_equal(out$interim,
                 nb_count_pmf(s$n1 * s$dispersion1, s$n1 * s$rate1, out$amax),
                 tolerance = 1e-12)
  }
})

test_that("the count from the later patients is negative binomial", {
  for (s in settings) {
    out <- build_case(s)
    expect_equal(out$extra,
                 nb_count_pmf((s$n2 - s$n1) * s$dispersion1,
                              (s$n2 - s$n1) * s$rate1, out$smax),
                 tolerance = 1e-12)
  }
})

test_that("the marginal Phase 3 total is negative binomial", {
  for (s in settings) {
    out <- build_case(s)
    marginal <- as.numeric(out$interim %*% out$transition)
    expect_equal(marginal,
                 nb_count_pmf(s$n3 * s$dispersion2, s$n3 * s$rate2, out$tmax),
                 tolerance = 1e-6)
  }
})

test_that("the propagated correlation matches the closed form", {
  for (s in settings) {
    out <- build_case(s)
    interim <- seq(0, out$amax)
    final <- seq(0, out$tmax)
    pmf_final <- as.numeric(out$interim %*% out$transition)
    mean_interim <- sum(out$interim * interim)
    mean_final <- sum(pmf_final * final)
    sd_interim <- sqrt(sum(out$interim * interim^2) - mean_interim^2)
    sd_final <- sqrt(sum(pmf_final * final^2) - mean_final^2)
    cross <- sum(out$interim * interim * as.numeric(out$transition %*% final))
    observed <- (cross - mean_interim * mean_final) / (sd_interim * sd_final)
    covariance <- s$rate1 * s$rate2 * s$shared / (s$dispersion1 * s$dispersion2)
    variance1 <- s$rate1 * (1 + s$rate1 / s$dispersion1)
    variance2 <- s$rate2 * (1 + s$rate2 / s$dispersion2)
    expected <- s$n1 * covariance /
      sqrt(s$n1 * variance1 * s$n3 * variance2)
    expect_equal(observed, expected, tolerance = 1e-5)
  }
})

test_that("a shared shape of zero makes the two endpoints independent", {
  out <- count_transition(rate1 = 2, rate2 = 0.5, dispersion1 = 0.6,
                          dispersion2 = 1.2, shared = 0,
                          n1 = 20, n2 = 30, n3 = 60,
                          amax = 100, smax = 150, tmax = 100)
  spread <- apply(out$transition, 2, function(z) diff(range(z)))
  expect_equal(max(spread), 0, tolerance = 1e-12)
})

test_that("no patients are added after the interim analysis", {
  out <- count_transition(rate1 = 2, rate2 = 0.5, dispersion1 = 0.6,
                          dispersion2 = 1.2, shared = 0.4,
                          n1 = 20, n2 = 20, n3 = 60,
                          amax = 120, smax = 150, tmax = 100)
  expect_equal(out$extra, c(1, rep(0, 150)))
})

test_that("a single cohort serves all three analyses", {
  out <- count_transition(rate1 = 2, rate2 = 0.5, dispersion1 = 0.6,
                          dispersion2 = 1.2, shared = 0.4,
                          n1 = 20, n2 = 20, n3 = 20,
                          amax = 120, smax = 150, tmax = 100)
  expect_equal(out$extra, c(1, rep(0, 150)))
  expect_equal(rowSums(out$transition), rep(1, 121), tolerance = 1e-8)
})

test_that("the returned pieces have the requested lengths", {
  out <- count_transition(rate1 = 2, rate2 = 0.5, dispersion1 = 0.6,
                          dispersion2 = 1.2, shared = 0.4,
                          n1 = 20, n2 = 30, n3 = 60,
                          amax = 90, smax = 110, tmax = 70)
  expect_length(out$interim, 91)
  expect_length(out$extra, 111)
  expect_equal(dim(out$transition), c(91L, 71L))
})

test_that("invalid arguments are rejected", {
  expect_error(count_transition(2, 0.5, 0.6, 1.2, 2, 20, 30, 60, 100, 150, 100),
               "smaller of the two dispersions")
  expect_error(count_transition(2, 0.5, 0.6, 1.2, -1, 20, 30, 60, 100, 150, 100),
               "smaller of the two dispersions")
  expect_error(count_transition(2, 0.5, 0.6, 1.2, 0.4, 40, 30, 60, 100, 150, 100),
               "non-decreasing")
})
