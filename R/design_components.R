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
#' the information per patient of Tang (2015). The covariance between two
#' statistics computed on shared patients follows from the covariance factor of
#' \code{coupling_factor()}. Each analysis conditions on the realised exposure
#' through its offset, so only the covariance of the two counts given the
#' exposure contributes. When the two endpoints count different events this is
#' the frailty they share, and the follow-up period they share does not
#' correlate the statistics. When they count the same events over the two
#' windows, the data behind the earlier statistic are part of the data behind
#' the later one, so the covariance equals the variance of the later statistic,
#' which is the canonical joint distribution of group sequential theory.
#'
#'
#' The negative binomial is written as in Homma and Daimon (2019), so
#' \code{dispersion1} and \code{dispersion2} are the shapes of the gamma
#' frailty and a count with mean \code{mu} has variance
#' \code{mu + mu^2 / dispersion}. This is the \code{size} argument of
#' \code{rnbinom()}. Tang (2015) and Zhu and Lakkis (2014) instead report the
#' variance of the frailty, which is the reciprocal, so a dispersion of 2.4 in
#' their sense is \code{dispersion = 1 / 2.4} here.
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
#' @param endpoints Relationship between the two endpoints, either
#'   \code{"distinct"} or \code{"shared"}. Under \code{"distinct"} the two
#'   endpoints count different events and are independent given the frailty.
#'   Under \code{"shared"} they count the same events over the two windows, so
#'   a single process governs both and the arguments \code{rate1},
#'   \code{dispersion1}, \code{shared} and \code{ratio1} are taken from
#'   \code{rate2}, \code{dispersion2} and \code{ratio2}.
#' @return A list with the three variances in \code{variance}, the three means
#'   in \code{mean}, the two correlations in \code{correlation}, the information
#'   per patient in \code{information}, and the covariance factor of
#'   \code{coupling_factor()} for each group in \code{coupling}.
#' @export
#' @examples
#' design_components(rate1 = 2.6, rate2 = 1.3, dispersion1 = 0.8,
#'                   dispersion2 = 0.8, shared = 0.4, n1 = 60, n2 = 90,
#'                   n3 = 300, window1 = 0.5, window2 = 1)
#'
#' # The same event counted over half a year at the interim analysis and over
#' # a year at the final analysis.
#' design_components(rate1 = 1.3, rate2 = 1.3, dispersion1 = 0.8,
#'                   dispersion2 = 0.8, shared = 0.8, n1 = 60, n2 = 90,
#'                   n3 = 300, window1 = 0.5, window2 = 1,
#'                   endpoints = "shared")
design_components <- function(rate1, rate2, dispersion1, dispersion2, shared,
                              n1, n2, n3, window1 = 1, window2 = 1,
                              dropout = 0,
                              dropout_shape = c("complete", "uniform",
                                                "exponential"),
                              ratio1 = 1, ratio2 = 1,
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
    coupling[g] <- coupling_factor(rate[g, 1], rate[g, 2], dispersion1,
                                   dispersion2, shared, window1, window2,
                                   dropout[g], dropout_shape, endpoints)
  }
  small_column <- if (small_endpoint == "phase2") 1L else 2L
  variance <- c(interim = sum(1 / (n1 * information[, 1])),
                small = sum(1 / (n2 * information[, small_column])),
                large = sum(1 / (n3 * information[, 2])))
  # When the two endpoints count the same events, the data behind the interim
  # statistic are part of the data behind either final statistic, so the
  # covariance is the variance of the final statistic, which is also what the
  # covariance factor gives in that case. Otherwise only the patients shared by
  # the two analyses contribute, through the covariance factor.
  covariance_small <- if (endpoints == "shared" ||
                          small_endpoint == "phase2") {
    variance[["small"]]
  } else {
    sum(coupling / n2)
  }
  covariance_large <- if (endpoints == "shared") {
    variance[["large"]]
  } else {
    sum(coupling / n3)
  }
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
       dropout_shape = dropout_shape, small_endpoint = small_endpoint,
       endpoints = endpoints)
}
