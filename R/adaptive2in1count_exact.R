#' Exact operating characteristics under complete follow-up
#'
#' Evaluates a two-group 2-in-1 adaptive design in which the Phase 2 part and the
#' Phase 3 part use different overdispersed count endpoints. The trial starts as
#' a small trial. At the interim analysis the observed rate ratio of the Phase 2
#' endpoint is compared with a cutpoint. If the rate ratio is at or below the
#' cutpoint the trial expands into a large trial and the final analysis uses the
#' Phase 3 endpoint on all patients, otherwise the trial remains small and the
#' final analysis uses the endpoint named by \code{small_endpoint}. Both final
#' analyses are carried out at the full one-sided level.
#'
#' Counts follow an additive gamma frailty model, so each endpoint is marginally
#' negative binomial with its own dispersion while the two endpoints are
#' correlated within a patient. All probabilities are obtained by propagating the
#' joint distribution of the interim total and the final total, so the reported
#' operating characteristics carry no Monte Carlo error and are exact up to the
#' truncation tolerance \code{tol}.
#'
#' Only the structure in which the two endpoints count different events is
#' available here. A design whose two endpoints count the same events over the
#' two windows, which \code{adaptive2in1count()} evaluates with
#' \code{endpoints = "shared"}, is not covered by this function.
#'
#' A lower count indicates greater efficacy, so the investigational group is the
#' first group throughout and the null hypothesis is rejected when its rate is
#' significantly lower.
#'
#' @param rate1 Mean count per patient for the Phase 2 endpoint in the control
#'   group, measured over the Phase 2 observation window.
#' @param rate2 Mean count per patient for the Phase 3 endpoint in the control
#'   group, measured over the Phase 3 observation window.
#' @param dispersion1 Dispersion of the Phase 2 endpoint.
#' @param dispersion2 Dispersion of the Phase 3 endpoint.
#' @param shared Shape of the frailty component shared by the two endpoints,
#'   between 0 and the smaller of the two dispersions.
#' @param n1 Number of patients per group at the interim analysis, either a
#'   single value or a vector of length two giving the investigational and the
#'   control group.
#' @param n2 Number of patients per group in the small trial, in the same format
#'   as \code{n1}.
#' @param n3 Number of patients per group in the large trial, in the same format
#'   as \code{n1}.
#' @param cutpoint Numeric vector of cutpoints for the observed rate ratio of the
#'   Phase 2 endpoint. The trial expands when the observed rate ratio is at or
#'   below the cutpoint.
#' @param alpha One-sided significance level applied at each final analysis.
#' @param ratio1 True rate ratio of the investigational group to the control
#'   group for the Phase 2 endpoint.
#' @param ratio2 True rate ratio of the investigational group to the control
#'   group for the Phase 3 endpoint.
#' @param test Final analysis to evaluate, either \code{"exact"} for the exact
#'   conditional test or \code{"wald"} for the Wald test on the log rate ratio.
#' @param small_endpoint Endpoint analysed at the final analysis of the small
#'   trial, either \code{"phase2"} for the endpoint that drives the interim
#'   decision or \code{"phase3"} for the endpoint of the large trial. The
#'   second choice keeps the confirmatory endpoint the same in both arms of the
#'   decision and leaves the Phase 2 endpoint as an intermediate endpoint used
#'   only to decide whether to expand.
#' @param tol Upper tail mass that may be discarded when truncating each support.
#' @return An object of class \code{adaptive2in1count_exact}, which is a list with the
#'   data frame \code{results} holding one row per cutpoint with the probability
#'   of expansion, the overall rejection probability, its split into the small
#'   trial and the large trial, and the expected total sample size, together with
#'   the design inputs and the retained supports.
#' @export
#' @examples
#' fit <- adaptive2in1count_exact(rate1 = 3, rate2 = 0.6, dispersion1 = 0.5,
#'                          dispersion2 = 1.5, shared = 0.25,
#'                          n1 = 50, n2 = 80, n3 = 250,
#'                          cutpoint = c(0.5, 0.8), alpha = 0.025)
#' fit
#'
#' adaptive2in1count_exact(rate1 = 3, rate2 = 0.6, dispersion1 = 0.5,
#'                   dispersion2 = 1.5, shared = 0.25,
#'                   n1 = 50, n2 = 80, n3 = 250, cutpoint = 0.5,
#'                   small_endpoint = "phase3")
#'
#' adaptive2in1count_exact(rate1 = 3, rate2 = 0.6, dispersion1 = 0.5,
#'                   dispersion2 = 1.5, shared = 0.25,
#'                   n1 = 50, n2 = 80, n3 = 250, cutpoint = 0.5,
#'                   ratio1 = 0.4, ratio2 = 0.7)
adaptive2in1count_exact <- function(rate1, rate2, dispersion1, dispersion2, shared,
                              n1, n2, n3, cutpoint, alpha = 0.025,
                              ratio1 = 1, ratio2 = 1,
                              test = c("exact", "wald"),
                              small_endpoint = c("phase2", "phase3"),
                              tol = 1e-11) {
  test <- match.arg(test)
  small_endpoint <- match.arg(small_endpoint)
  for (nm in c("rate1", "rate2", "dispersion1", "dispersion2",
               "ratio1", "ratio2")) {
    value <- get(nm)
    if (!is.numeric(value) || length(value) != 1L || is.na(value) || value <= 0) {
      stop(sprintf("'%s' must be a single positive numeric value.", nm))
    }
  }
  if (!is.numeric(shared) || length(shared) != 1L || is.na(shared) ||
      shared < 0 || shared > min(dispersion1, dispersion2)) {
    stop("'shared' must lie between 0 and the smaller of the two dispersions.")
  }
  if (!is.numeric(alpha) || length(alpha) != 1L || is.na(alpha) ||
      alpha <= 0 || alpha >= 1) {
    stop("'alpha' must be a single numeric value between 0 and 1.")
  }
  if (!is.numeric(cutpoint) || length(cutpoint) == 0L || anyNA(cutpoint) ||
      any(cutpoint <= 0)) {
    stop("'cutpoint' must be a numeric vector of positive values.")
  }
  n1 <- expand_group_size(n1, "n1")
  n2 <- expand_group_size(n2, "n2")
  n3 <- expand_group_size(n3, "n3")
  if (any(n1 > n2) || any(n2 > n3)) {
    stop("'n1', 'n2' and 'n3' must be non-decreasing within each group.")
  }
  mean1 <- c(rate1 * ratio1, rate1)
  mean2 <- c(rate2 * ratio2, rate2)
  amax <- max(vapply(seq_len(2), function(g) {
    nb_support_max(n1[g] * dispersion1, n1[g] * mean1[g], tol)
  }, integer(1)))
  small_dispersion <- if (small_endpoint == "phase2") {
    dispersion1
  } else {
    dispersion2
  }
  small_mean <- if (small_endpoint == "phase2") mean1 else mean2
  smax <- max(vapply(seq_len(2), function(g) {
    nb_support_max(n2[g] * small_dispersion, n2[g] * small_mean[g], tol)
  }, integer(1)))
  tmax <- max(vapply(seq_len(2), function(g) {
    nb_support_max(n3[g] * dispersion2, n3[g] * mean2[g], tol)
  }, integer(1)))
  extra_max <- if (small_endpoint == "phase2") smax else 0L
  part <- lapply(seq_len(2), function(g) {
    count_transition(rate1 = mean1[g], rate2 = mean2[g],
                     dispersion1 = dispersion1, dispersion2 = dispersion2,
                     shared = shared, n1 = n1[g], n2 = n2[g], n3 = n3[g],
                     amax = amax, smax = extra_max, tmax = tmax)
  })
  if (test == "exact") {
    region_small <- exact_conditional_region(n2[1] * small_dispersion,
                                             n2[2] * small_dispersion, smax,
                                             alpha)
    region_large <- exact_conditional_region(n3[1] * dispersion2,
                                             n3[2] * dispersion2, tmax, alpha)
  } else {
    region_small <- wald_region(n2[1], n2[2], small_dispersion, smax, alpha)
    region_large <- wald_region(n3[1], n3[2], dispersion2, tmax, alpha)
  }
  # When the small trial analyses the Phase 2 endpoint, its total is the interim
  # count plus an independent contribution from the patients enrolled later.
  # When it analyses the Phase 3 endpoint, its total is linked to the interim
  # count through the shared frailty in exactly the way the large trial total
  # is, so the same transition matrix applies with the small trial size in place
  # of the large trial size.
  if (small_endpoint == "phase2") {
    link <- lapply(seq_len(2), function(g) {
      shift_matrix_cpp(part[[g]]$extra, as.integer(amax), as.integer(smax))
    })
  } else {
    link <- lapply(seq_len(2), function(g) {
      count_transition(rate1 = mean1[g], rate2 = mean2[g],
                       dispersion1 = dispersion1, dispersion2 = dispersion2,
                       shared = shared, n1 = n1[g], n2 = n1[g], n3 = n2[g],
                       amax = amax, smax = 0L, tmax = smax)$transition
    })
  }
  power_small <- link[[1]] %*% region_small %*% t(link[[2]])
  power_large <- part[[1]]$transition %*% region_large %*%
    t(part[[2]]$transition)
  joint <- outer(part[[1]]$interim, part[[2]]$interim)
  results <- data.frame(cutpoint = cutpoint, go_probability = NA_real_,
                        rejection = NA_real_, rejection_small = NA_real_,
                        rejection_large = NA_real_, expected_size = NA_real_)
  for (i in seq_along(cutpoint)) {
    go <- expansion_region(amax, n1, cutpoint[i])
    go_probability <- sum(joint * go)
    small <- sum(joint * (!go) * power_small)
    large <- sum(joint * go * power_large)
    results$go_probability[i] <- go_probability
    results$rejection[i] <- small + large
    results$rejection_small[i] <- small
    results$rejection_large[i] <- large
    results$expected_size[i] <- (1 - go_probability) * sum(n2) +
      go_probability * sum(n3)
  }
  structure(list(results = results, rate1 = rate1, rate2 = rate2,
                 dispersion1 = dispersion1, dispersion2 = dispersion2,
                 shared = shared, n1 = n1, n2 = n2, n3 = n3, alpha = alpha,
                 ratio1 = ratio1, ratio2 = ratio2, test = test,
                 small_endpoint = small_endpoint,
                 support = c(interim = amax, small = smax, large = tmax),
                 call = match.call()),
            class = "adaptive2in1count_exact")
}
