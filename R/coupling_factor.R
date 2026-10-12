#' Covariance factor of two log rates estimated on shared patients
#'
#' Returns, for one group, the factor \code{c} such that the covariance of the
#' log rate of the Phase 2 endpoint estimated from \code{m} patients and the log
#' rate of the Phase 3 endpoint estimated from \code{n} patients of a nested set
#' is approximately \code{min(m, n) / (m * n) * c}.
#'
#' Each log rate is estimated by a negative binomial model with an offset for
#' the realised exposure. The score residual of a patient,
#' \code{count - rate * exposure}, has mean zero given the exposure, so the
#' covariance of the two estimators is that of their influence functions,
#' \code{E(w1 * w2 * Cov(Y1, Y2 | T)) / (kappa1 * kappa2)}, where
#' \code{w_k = dispersion_k / (dispersion_k + rate_k * T_k)} is the weight of
#' the score, \code{kappa_k} is the information per patient returned by
#' \code{information_factor()} and \code{T_k = min(T, window_k)} is the
#' exposure. Only the covariance of the two counts given the exposure enters.
#' The covariance of the two exposures, which a random follow-up period adds to
#' the covariance of the counts themselves, does not reach the estimators.
#'
#' Under \code{endpoints = "distinct"} the covariance given the exposure is
#' \code{rate1 * rate2 * T1 * T2 * shared / (dispersion1 * dispersion2)}, so the
#' factor is zero when no frailty is shared and reduces to
#' \code{shared / (dispersion1 * dispersion2)} under complete follow-up. Under
#' \code{endpoints = "shared"} the Phase 2 count is part of the Phase 3 count,
#' the covariance given the exposure is
#' \code{rate * T1 + rate^2 * T1 * T2 / dispersion}, and the factor equals
#' \code{1 / kappa2}, which gives the canonical covariance of group sequential
#' theory.
#'
#' @param rate1 Event rate per unit of exposure for the Phase 2 endpoint.
#' @param rate2 Event rate per unit of exposure for the Phase 3 endpoint.
#' @param dispersion1 Dispersion of the Phase 2 endpoint.
#' @param dispersion2 Dispersion of the Phase 3 endpoint.
#' @param shared Shape of the frailty component shared by the two endpoints.
#' @param window1 Observation window of the Phase 2 endpoint.
#' @param window2 Observation window of the Phase 3 endpoint, which is also the
#'   follow-up period.
#' @param dropout Proportion of patients who drop out before the end of the
#'   follow-up period.
#' @param dropout_shape Distribution of the follow-up period, one of
#'   \code{"complete"}, \code{"uniform"} and \code{"exponential"}.
#' @param endpoints Relationship between the two endpoints, either
#'   \code{"distinct"} or \code{"shared"}. Under \code{"shared"} the arguments
#'   \code{rate1}, \code{dispersion1} and \code{shared} are taken from
#'   \code{rate2} and \code{dispersion2}.
#' @return A numeric value holding the covariance factor.
#' @export
#' @examples
#' # Complete follow-up: the frailty term alone.
#' coupling_factor(rate1 = 2, rate2 = 1.3, dispersion1 = 0.7,
#'                 dispersion2 = 0.8, shared = 0.4, window1 = 0.5)
#'
#' # No shared frailty: the two estimators are uncorrelated even when patients
#' # drop out.
#' coupling_factor(rate1 = 2, rate2 = 1.3, dispersion1 = 0.7,
#'                 dispersion2 = 0.8, shared = 0, window1 = 1, dropout = 0.6,
#'                 dropout_shape = "exponential")
coupling_factor <- function(rate1, rate2, dispersion1, dispersion2, shared,
                            window1 = 1, window2 = 1, dropout = 0,
                            dropout_shape = c("complete", "uniform",
                                              "exponential"),
                            endpoints = c("distinct", "shared")) {
  dropout_shape <- match.arg(dropout_shape)
  endpoints <- match.arg(endpoints)
  if (endpoints == "shared") {
    rate1 <- rate2
    dispersion1 <- dispersion2
    shared <- dispersion2
  }
  if (!is.numeric(shared) || length(shared) != 1L || is.na(shared) ||
      shared < 0 || shared > min(dispersion1, dispersion2)) {
    stop("'shared' must lie between 0 and the smaller of the two dispersions.")
  }
  if (!is.numeric(window1) || length(window1) != 1L || is.na(window1) ||
      window1 <= 0 || window1 > window2) {
    stop("'window1' must be a single positive value not exceeding 'window2'.")
  }
  information1 <- information_factor(rate1, dispersion1, window1, window2,
                                     dropout, dropout_shape)
  information2 <- information_factor(rate2, dispersion2, window2, window2,
                                     dropout, dropout_shape)
  # Covariance of the two counts given the follow-up period, multiplied by the
  # weights of the two scores.
  integrand <- function(t) {
    exposure1 <- pmin(t, window1)
    exposure2 <- pmin(t, window2)
    weight1 <- dispersion1 / (dispersion1 + rate1 * exposure1)
    weight2 <- dispersion2 / (dispersion2 + rate2 * exposure2)
    conditional <- if (endpoints == "shared") {
      rate2 * exposure1 + rate2^2 * exposure1 * exposure2 / dispersion2
    } else {
      rate1 * rate2 * exposure1 * exposure2 * shared /
        (dispersion1 * dispersion2)
    }
    weight1 * weight2 * conditional
  }
  follow_up_expectation(integrand, window2, dropout, dropout_shape) /
    (information1 * information2)
}
