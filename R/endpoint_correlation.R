#' Correlation between the two count endpoints within a patient
#'
#' Returns the correlation between the Phase 2 endpoint count and the Phase 3
#' endpoint count of the same patient. Two mechanisms contribute when the two
#' endpoints count different events. The frailty component shared by the two
#' endpoints correlates the underlying event rates, and a random follow-up
#' period correlates the exposures because a patient who leaves the trial stops
#' contributing to either count. The second mechanism operates even when the
#' frailty is not shared. When the two endpoints count the same events over the
#' two windows, a third mechanism is added, since the Phase 2 count is then
#' part of the Phase 3 count.
#'
#' @param rate1 Event rate per unit of exposure for the Phase 2 endpoint.
#' @param rate2 Event rate per unit of exposure for the Phase 3 endpoint.
#' @param dispersion1 Dispersion of the Phase 2 endpoint.
#' @param dispersion2 Dispersion of the Phase 3 endpoint.
#' @param shared Shape of the frailty component shared by the two endpoints.
#' @param window1 Observation window of the Phase 2 endpoint.
#' @param window2 Observation window of the Phase 3 endpoint.
#' @param dropout Proportion of patients who drop out before the end of the
#'   follow-up period.
#' @param dropout_shape Distribution of the follow-up period.
#' @param endpoints Relationship between the two endpoints, either
#'   \code{"distinct"} or \code{"shared"}. Under \code{"distinct"} the two
#'   endpoints count different events and are independent given the frailty,
#'   which is the additive gamma frailty model. Under \code{"shared"} they
#'   count the same events over the two windows, so the Phase 2 count is part
#'   of the Phase 3 count, and the covariance gains the variance of the shared
#'   events. In that case a single process governs both endpoints, and the
#'   arguments \code{rate1}, \code{dispersion1} and \code{shared} are taken
#'   from \code{rate2} and \code{dispersion2}.
#' @return A numeric value holding the within-patient correlation.
#' @export
#' @examples
#' endpoint_correlation(rate1 = 2.6, rate2 = 1.3, dispersion1 = 0.8,
#'                      dispersion2 = 0.8, shared = 0.4,
#'                      window1 = 0.5, window2 = 1)
#'
#' endpoint_correlation(rate1 = 2.6, rate2 = 1.3, dispersion1 = 0.8,
#'                      dispersion2 = 0.8, shared = 0,
#'                      window1 = 0.5, window2 = 1, dropout = 0.4,
#'                      dropout_shape = "exponential")
#'
#' # The same event counted over half a year and over a year.
#' endpoint_correlation(rate1 = 1.3, rate2 = 1.3, dispersion1 = 0.8,
#'                      dispersion2 = 0.8, shared = 0.8, window1 = 0.5,
#'                      window2 = 1, endpoints = "shared")
endpoint_correlation <- function(rate1, rate2, dispersion1, dispersion2, shared,
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
  if (shared < 0 || shared > min(dispersion1, dispersion2)) {
    stop("'shared' must lie between 0 and the smaller of the two dispersions.")
  }
  moments <- follow_up_moments(window1, window2, dropout, dropout_shape)
  covariance <- rate1 * rate2 *
    (moments$cross * (1 + shared / (dispersion1 * dispersion2)) -
       moments$first * moments$second)
  if (endpoints == "shared") {
    covariance <- covariance + rate1 * moments$first
  }
  variance1 <- rate1 * moments$first +
    rate1^2 * moments$square1 * (1 + 1 / dispersion1) -
    rate1^2 * moments$first^2
  variance2 <- rate2 * moments$second +
    rate2^2 * moments$square2 * (1 + 1 / dispersion2) -
    rate2^2 * moments$second^2
  covariance / sqrt(variance1 * variance2)
}
