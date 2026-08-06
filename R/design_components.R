#' Variances and correlations of the three test statistics
#'
#' Assembles the quantities that the normal approximation of the 2-in-1 design
#' requires. The interim statistic is the log rate ratio of the Phase 2 endpoint
#' among the patients followed at the interim analysis, the small trial
#' statistic is the log rate ratio of the endpoint named by
#' \code{small_endpoint} among the patients of the small trial, and the large
#' trial statistic is the log rate ratio of the Phase 3 endpoint among all
#' patients.
#'
#' The variance of a log rate is the reciprocal of the number of patients times
#' the information per patient of Tang (2015). The covariance between the two
#' endpoints within a patient is driven by the frailty they share and by the
#' follow-up period they share, and the second contribution is present whenever
#' the follow-up period is random.
#'
#' @param rate1 Event rate per unit of exposure for the Phase 2 endpoint in the
#'   control group.
#' @param rate2 Event rate per unit of exposure for the Phase 3 endpoint in the
#'   control group.
#' @param dispersion1 Dispersion of the Phase 2 endpoint.
#' @param dispersion2 Dispersion of the Phase 3 endpoint.
#' @param shared Shape of the frailty component shared by the two endpoints.
#' @param n1 Number of patients per group at the interim analysis.
#' @param n2 Number of patients per group in the small trial.
#' @param n3 Number of patients per group in the large trial.
#' @param window1 Observation window of the Phase 2 endpoint.
#' @param window2 Observation window of the Phase 3 endpoint.
#' @param dropout Proportion of patients who drop out, either a single value or
#'   a vector of length two giving the investigational and the control group.
#' @param dropout_shape Distribution of the follow-up period.
#' @param ratio1 True rate ratio for the Phase 2 endpoint.
#' @param ratio2 True rate ratio for the Phase 3 endpoint.
#' @param small_endpoint Endpoint analysed in the small trial.
#' @return A list with the three variances in \code{variance}, the three means
#'   in \code{mean}, the two correlations in \code{correlation}, the information
#'   per patient in \code{information}, and the within-patient covariance factor
#'   of the two endpoints in \code{coupling}.
#' @export
#' @examples
#' design_components(rate1 = 2.6, rate2 = 1.3, dispersion1 = 0.8,
#'                   dispersion2 = 0.8, shared = 0.4, n1 = 60, n2 = 90,
#'                   n3 = 300, window1 = 0.5, window2 = 1)
design_components <- function(rate1, rate2, dispersion1, dispersion2, shared,
                              n1, n2, n3, window1 = 1, window2 = 1,
                              dropout = 0,
                              dropout_shape = c("complete", "uniform",
                                                "exponential"),
                              ratio1 = 1, ratio2 = 1,
                              small_endpoint = c("phase2", "phase3")) {
  dropout_shape <- match.arg(dropout_shape)
  small_endpoint <- match.arg(small_endpoint)
  if (shared < 0 || shared > min(dispersion1, dispersion2)) {
    stop("'shared' must lie between 0 and the smaller of the two dispersions.")
  }
  n1 <- expand_group_size(n1, "n1")
  n2 <- expand_group_size(n2, "n2")
  n3 <- expand_group_size(n3, "n3")
  if (any(n1 > n2) || any(n2 > n3)) {
    stop("'n1', 'n2' and 'n3' must be non-decreasing within each group.")
  }
  if (length(dropout) == 1L) {
    dropout <- rep(dropout, 2L)
  }
  if (length(dropout) != 2L || anyNA(dropout) || any(dropout < 0) ||
      any(dropout >= 1)) {
    stop("'dropout' must be one or two values in the interval from 0 to 1.")
  }
  rate <- cbind(c(rate1 * ratio1, rate1), c(rate2 * ratio2, rate2))
  window <- c(window1, window2)
  information <- matrix(NA_real_, 2L, 2L)
  coupling <- numeric(2L)
  for (g in seq_len(2)) {
    for (k in seq_len(2)) {
      information[g, k] <- information_factor(rate[g, k], c(dispersion1,
                                                            dispersion2)[k],
                                              window[k], window2, dropout[g],
                                              dropout_shape)
    }
    moments <- follow_up_moments(window1, window2, dropout[g], dropout_shape)
    coupling[g] <- moments$inflation *
      (1 + shared / (dispersion1 * dispersion2)) - 1
  }
  small_column <- if (small_endpoint == "phase2") 1L else 2L
  variance <- c(interim = sum(1 / (n1 * information[, 1])),
                small = sum(1 / (n2 * information[, small_column])),
                large = sum(1 / (n3 * information[, 2])))
  covariance_small <- if (small_endpoint == "phase2") {
    variance[["small"]]
  } else {
    sum(coupling / n2)
  }
  covariance_large <- sum(coupling / n3)
  correlation <- c(
    interim_small = covariance_small /
      sqrt(variance[["interim"]] * variance[["small"]]),
    interim_large = covariance_large /
      sqrt(variance[["interim"]] * variance[["large"]]))
  effect <- c(interim = log(ratio1),
              small = if (small_endpoint == "phase2") log(ratio1) else log(ratio2),
              large = log(ratio2))
  list(variance = variance, mean = effect, correlation = correlation,
       information = information, coupling = coupling, dropout = dropout,
       dropout_shape = dropout_shape, small_endpoint = small_endpoint)
}
