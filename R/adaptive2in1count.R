#' Operating characteristics of a 2-in-1 adaptive design for count endpoints
#'
#' Evaluates a two-group 2-in-1 adaptive design in which the Phase 2 part and
#' the Phase 3 part use overdispersed count endpoints and patients may leave the
#' trial before the end of the follow-up period. The trial starts as a small
#' trial. At the interim analysis the observed rate ratio of the Phase 2
#' endpoint is compared with a cutpoint. If the observed rate ratio is at or
#' below the cutpoint the trial expands into a large trial and the final
#' analysis uses the Phase 3 endpoint on all patients, otherwise the trial
#' remains small and the final analysis uses the endpoint named by
#' \code{small_endpoint}. Both final analyses are Wald tests of the log rate
#' ratio from a negative binomial model with an offset for the realised
#' exposure, carried out at the full one-sided level.
#'
#' Counts follow an additive gamma frailty model, so each endpoint is marginally
#' negative binomial with its own dispersion. The two endpoints of a patient are
#' correlated through the frailty they share and through the follow-up period
#' they share. Each analysis conditions on the realised exposure through its
#' offset, so only the frailty correlates the three log rate ratios, and the
#' operating characteristics are obtained from their trivariate normal
#' distribution. Setting
#' \code{endpoints = "shared"} instead treats the two endpoints as the same
#' event counted over the two windows, which covers a design whose three
#' analyses read one endpoint.
#'
#'
#' The negative binomial is written as in Homma and Daimon (2019), so
#' \code{dispersion1} and \code{dispersion2} are the shapes of the gamma
#' frailty and a count with mean \code{mu} has variance
#' \code{mu + mu^2 / dispersion}. This is the \code{size} argument of
#' \code{rnbinom()}. Tang (2015) and Zhu and Lakkis (2014) instead report the
#' variance of the frailty, which is the reciprocal, so a dispersion of 2.4 in
#' their sense is \code{dispersion = 1 / 2.4} here.
#' A lower rate indicates greater efficacy, so the investigational group is the
#' first group throughout and the null hypothesis is rejected when its rate is
#' significantly lower.
#'
#' @param rate1 Event rate per unit of exposure for the Phase 2 endpoint in the
#'   control group.
#' @param rate2 Event rate per unit of exposure for the Phase 3 endpoint in the
#'   control group.
#' @param dispersion1 Dispersion of the Phase 2 endpoint.
#' @param dispersion2 Dispersion of the Phase 3 endpoint.
#' @param shared Shape of the frailty component shared by the two endpoints,
#'   between 0 and the smaller of the two dispersions.
#' @param n1 Number of patients per group at the interim analysis, either a
#'   single value or a vector of length two.
#' @param n2 Number of patients per group in the small trial.
#' @param n3 Number of patients per group in the large trial.
#' @param cutpoint Numeric vector of cutpoints for the observed rate ratio of
#'   the Phase 2 endpoint.
#' @param window1 Observation window of the Phase 2 endpoint.
#' @param window2 Observation window of the Phase 3 endpoint, which is also the
#'   prespecified follow-up period.
#' @param dropout Proportion of patients who drop out before the end of the
#'   follow-up period, either a single value or a vector of length two.
#' @param dropout_shape Distribution of the follow-up period, one of
#'   \code{"complete"}, \code{"uniform"} and \code{"exponential"}.
#' @param alpha One-sided significance level applied at each final analysis.
#' @param ratio1 True rate ratio for the Phase 2 endpoint.
#' @param ratio2 True rate ratio for the Phase 3 endpoint.
#' @param small_endpoint Endpoint analysed at the final analysis of the small
#'   trial, either   \code{"phase2"} or \code{"phase3"}.
#' @param endpoints Relationship between the two endpoints, either
#'   \code{"distinct"} or \code{"shared"}. Under \code{"distinct"} the two
#'   endpoints count different events and are independent given the frailty.
#'   Under \code{"shared"} they count the same events over the two windows, so
#'   the Phase 2 count is part of the Phase 3 count. A single process then
#'   governs both endpoints, and the arguments \code{rate1},
#'   \code{dispersion1}, \code{shared} and \code{ratio1} are taken from
#'   \code{rate2}, \code{dispersion2} and \code{ratio2}.
#' @return An object of class \code{adaptive2in1count} holding the data frame
#'   \code{results} with one row per cutpoint, together with the design inputs
#'   and the components returned by \code{design_components()}.
#' @importFrom stats pnorm qnorm
#' @export
#' @examples
#' adaptive2in1count(rate1 = 2.6, rate2 = 1.3, dispersion1 = 0.8,
#'                   dispersion2 = 0.8, shared = 0.4, n1 = 60, n2 = 90,
#'                   n3 = 300, cutpoint = c(0.6, 0.8, 1),
#'                   window1 = 0.5, window2 = 1, dropout = 0.3,
#'                   dropout_shape = "exponential")
#'
#' adaptive2in1count(rate1 = 2.6, rate2 = 1.3, dispersion1 = 0.8,
#'                   dispersion2 = 0.8, shared = 0.4, n1 = 60, n2 = 90,
#'                   n3 = 300, cutpoint = 0.8, window1 = 0.5, window2 = 1,
#'                   dropout = 0.3, dropout_shape = "exponential",
#'                   ratio1 = 0.65, ratio2 = 0.75)
adaptive2in1count <- function(rate1, rate2, dispersion1, dispersion2, shared,
                              n1, n2, n3, cutpoint, window1 = 1, window2 = 1,
                              dropout = 0,
                              dropout_shape = c("complete", "uniform",
                                                "exponential"),
                              alpha = 0.025, ratio1 = 1, ratio2 = 1,
                              small_endpoint = c("phase2", "phase3"),
                              endpoints = c("distinct", "shared")) {
  dropout_shape <- match.arg(dropout_shape)
  small_endpoint <- match.arg(small_endpoint)
  endpoints <- match.arg(endpoints)
  if (endpoints == "shared") {
    rate1 <- rate2
    dispersion1 <- dispersion2
    shared <- dispersion2
    ratio1 <- ratio2
  }
  for (nm in c("rate1", "rate2", "dispersion1", "dispersion2", "ratio1",
               "ratio2")) {
    value <- get(nm)
    if (!is.numeric(value) || length(value) != 1L || is.na(value) ||
        value <= 0) {
      stop(sprintf("'%s' must be a single positive numeric value.", nm))
    }
  }
  if (!is.numeric(alpha) || length(alpha) != 1L || is.na(alpha) ||
      alpha <= 0 || alpha >= 1) {
    stop("'alpha' must be a single numeric value between 0 and 1.")
  }
  if (!is.numeric(cutpoint) || length(cutpoint) == 0L || anyNA(cutpoint) ||
      any(cutpoint <= 0)) {
    stop("'cutpoint' must be a numeric vector of positive values.")
  }
  parts <- design_components(rate1 = rate1, rate2 = rate2,
                             dispersion1 = dispersion1,
                             dispersion2 = dispersion2, shared = shared,
                             n1 = n1, n2 = n2, n3 = n3, window1 = window1,
                             window2 = window2, dropout = dropout,
                             dropout_shape = dropout_shape, ratio1 = ratio1,
                             ratio2 = ratio2, small_endpoint = small_endpoint,
                             endpoints = endpoints)
  n1 <- expand_group_size(n1, "n1")
  n2 <- expand_group_size(n2, "n2")
  n3 <- expand_group_size(n3, "n3")
  deviation <- sqrt(parts$variance)
  boundary_small <- -qnorm(1 - alpha) - parts$mean[["small"]] /
    deviation[["small"]]
  boundary_large <- -qnorm(1 - alpha) - parts$mean[["large"]] /
    deviation[["large"]]
  reject_small <- pnorm(boundary_small)
  results <- data.frame(cutpoint = cutpoint, go_probability = NA_real_,
                        rejection = NA_real_, rejection_small = NA_real_,
                        rejection_large = NA_real_, expected_size = NA_real_)
  for (i in seq_along(cutpoint)) {
    boundary_interim <- (log(cutpoint[i]) - parts$mean[["interim"]]) /
      deviation[["interim"]]
    go_probability <- pnorm(boundary_interim)
    large <- bivariate_normal_lower(boundary_interim, boundary_large,
                                    parts$correlation[["interim_large"]])
    small <- reject_small -
      bivariate_normal_lower(boundary_interim, boundary_small,
                             parts$correlation[["interim_small"]])
    small <- max(small, 0)
    results$go_probability[i] <- go_probability
    results$rejection[i] <- small + large
    results$rejection_small[i] <- small
    results$rejection_large[i] <- large
    results$expected_size[i] <- (1 - go_probability) * sum(n2) +
      go_probability * sum(n3)
  }
  structure(list(results = results, rate1 = rate1, rate2 = rate2,
                 dispersion1 = dispersion1, dispersion2 = dispersion2,
                 shared = shared, n1 = n1, n2 = n2, n3 = n3,
                 window1 = window1, window2 = window2,
                 dropout = parts$dropout, dropout_shape = dropout_shape,
                 alpha = alpha, ratio1 = ratio1, ratio2 = ratio2,
                 small_endpoint = small_endpoint, endpoints = endpoints,
                 components = parts, call = match.call()),
            class = "adaptive2in1count")
}
