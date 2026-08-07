#' Information per patient for a negative binomial rate
#'
#' Returns the quantity that Tang (2015) denotes by the expectation of
#' \code{dispersion * rate * t / (dispersion + rate * t)} over the distribution
#' of the exposure \code{t}. The asymptotic variance of the log rate estimated
#' from a negative binomial model with an offset for the exposure is the
#' reciprocal of this quantity divided by the number of patients, so it plays
#' the role that the number of events plays for a Poisson rate.
#'
#' Under complete follow-up the factor reduces to
#' \code{dispersion * mean / (dispersion + mean)} with \code{mean} the expected
#' count over the whole window, which gives the familiar variance
#' \code{1 / (n * mean) + 1 / (n * dispersion)}.
#'
#' @param rate Event rate per unit of exposure.
#' @param dispersion Dispersion of the negative binomial distribution.
#' @param window Observation window of the endpoint.
#' @param follow_up Prespecified follow-up period, which must be at least
#'   \code{window}.
#' @param dropout Proportion of patients who drop out before the end of the
#'   follow-up period.
#' @param shape Distribution of the follow-up period, one of
#'   \code{"complete"}, \code{"uniform"} and \code{"exponential"}. The two
#'   random distributions both use \code{dropout} as the proportion who leave
#'   early and both reduce to complete follow-up when it is zero.
#' @return A numeric value holding the information per patient.
#' @export
#' @examples
#' information_factor(rate = 1.3, dispersion = 0.8, window = 1, follow_up = 1)
#' information_factor(rate = 1.3, dispersion = 0.8, window = 1, follow_up = 1,
#'                    dropout = 0.3, shape = "exponential")
information_factor <- function(rate, dispersion, window, follow_up,
                               dropout = 0,
                               shape = c("complete", "uniform",
                                         "exponential")) {
  shape <- match.arg(shape)
  if (!is.numeric(rate) || length(rate) != 1L || is.na(rate) || rate <= 0) {
    stop("'rate' must be a single positive numeric value.")
  }
  if (!is.numeric(dispersion) || length(dispersion) != 1L ||
      is.na(dispersion) || dispersion <= 0) {
    stop("'dispersion' must be a single positive numeric value.")
  }
  if (window > follow_up) {
    stop("'window' must not exceed 'follow_up'.")
  }
  follow_up_expectation(function(t) {
    exposure <- pmin(t, window)
    dispersion * rate * exposure / (dispersion + rate * exposure)
  }, follow_up, dropout, shape)
}
