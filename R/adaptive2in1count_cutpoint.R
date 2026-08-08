#' Cutpoint that attains a target probability of expansion
#'
#' Returns the cutpoint for the observed rate ratio of the Phase 2 endpoint
#' whose probability of expansion equals \code{target} under the stated rate
#' ratio. Under the normal approximation the interim log rate ratio is normal,
#' so the cutpoint is available in closed form.
#'
#' @inheritParams design_components
#' @param target Target probability of expansion.
#' @param endpoints Relationship between the two endpoints, either
#'   \code{"distinct"} or \code{"shared"}. Under \code{"distinct"} the two
#'   endpoints count different events and are independent given the frailty.
#'   Under \code{"shared"} they count the same events over the two windows, so
#'   the Phase 2 count is part of the Phase 3 count. A single process then
#'   governs both endpoints, and the arguments \code{rate1},
#'   \code{dispersion1}, \code{shared} and \code{ratio1} are taken from
#'   \code{rate2}, \code{dispersion2} and \code{ratio2}.
#' @return A list with the selected \code{cutpoint} and the attained
#'   \code{go_probability}.
#' @importFrom stats qnorm
#' @export
#' @examples
#' adaptive2in1count_cutpoint(rate1 = 2.6, rate2 = 1.3, dispersion1 = 0.8,
#'                            dispersion2 = 0.8, shared = 0.4, n1 = 60,
#'                            n2 = 90, n3 = 300, target = 0.5, ratio1 = 0.65,
#'                            window1 = 0.5, window2 = 1, dropout = 0.3,
#'                            dropout_shape = "exponential")
adaptive2in1count_cutpoint <- function(rate1, rate2, dispersion1, dispersion2,
                                       shared, n1, n2, n3, target,
                                       window1 = 1, window2 = 1, dropout = 0,
                                       dropout_shape = c("complete", "uniform",
                                                         "exponential"),
                                       ratio1 = 1, ratio2 = 1,
                                       small_endpoint = c("phase2",
                                                          "phase3"),
                                       endpoints = c("distinct", "shared")) {
  dropout_shape <- match.arg(dropout_shape)
  small_endpoint <- match.arg(small_endpoint)
  endpoints <- match.arg(endpoints)
  if (!is.numeric(target) || length(target) != 1L || is.na(target) ||
      target <= 0 || target >= 1) {
    stop("'target' must be a single numeric value between 0 and 1.")
  }
  parts <- design_components(rate1 = rate1, rate2 = rate2,
                             dispersion1 = dispersion1,
                             dispersion2 = dispersion2, shared = shared,
                             n1 = n1, n2 = n2, n3 = n3, window1 = window1,
                             window2 = window2, dropout = dropout,
                             dropout_shape = dropout_shape, ratio1 = ratio1,
                             ratio2 = ratio2, small_endpoint = small_endpoint,
                             endpoints = endpoints)
  cutpoint <- exp(parts$mean[["interim"]] +
                    sqrt(parts$variance[["interim"]]) * qnorm(target))
  list(cutpoint = cutpoint, go_probability = target)
}
