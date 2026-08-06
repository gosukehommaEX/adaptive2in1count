#' Moments of the exposure of the two count endpoints
#'
#' A single follow-up period determines the exposure of both endpoints, because
#' a patient who leaves the trial stops contributing to either count. With
#' observation windows \code{window1} and \code{window2} the exposures are
#' \code{min(T, window1)} and \code{min(T, window2)}, where \code{T} is the
#' realised follow-up period. The two exposures are increasing functions of the
#' same variable, so they are positively correlated whenever the follow-up
#' period is random. That correlation is a second source of dependence between
#' the two endpoints, separate from the frailty they share.
#'
#' @param window1 Observation window of the Phase 2 endpoint, which must not
#'   exceed \code{window2}.
#' @param window2 Observation window of the Phase 3 endpoint, which is also the
#'   prespecified follow-up period.
#' @param dropout Proportion of patients who drop out before the end of the
#'   follow-up period.
#' @param shape Distribution of the follow-up period, one of
#'   \code{"complete"}, \code{"uniform"} and \code{"exponential"}.
#' @return A list with the mean exposure of each endpoint in \code{first} and
#'   \code{second}, their second moments in \code{square1} and \code{square2},
#'   their product moment in \code{cross}, and the ratio
#'   \code{cross / (first * second)} in \code{inflation}. The inflation equals
#'   one under complete follow-up and exceeds one otherwise.
#' @export
#' @examples
#' follow_up_moments(window1 = 0.5, window2 = 1, shape = "complete")
#' follow_up_moments(window1 = 0.5, window2 = 1, dropout = 0.3,
#'                   shape = "exponential")
follow_up_moments <- function(window1, window2, dropout = 0,
                              shape = c("complete", "uniform",
                                        "exponential")) {
  shape <- match.arg(shape)
  if (!is.numeric(window1) || length(window1) != 1L || is.na(window1) ||
      window1 <= 0 || window1 > window2) {
    stop("'window1' must be a single positive value not exceeding 'window2'.")
  }
  first <- follow_up_expectation(function(t) pmin(t, window1), window2,
                                 dropout, shape)
  second <- follow_up_expectation(function(t) pmin(t, window2), window2,
                                  dropout, shape)
  cross <- follow_up_expectation(function(t) pmin(t, window1) *
                                   pmin(t, window2), window2, dropout, shape)
  square1 <- follow_up_expectation(function(t) pmin(t, window1)^2, window2,
                                   dropout, shape)
  square2 <- follow_up_expectation(function(t) pmin(t, window2)^2, window2,
                                   dropout, shape)
  list(first = first, second = second, square1 = square1, square2 = square2,
       cross = cross, inflation = cross / (first * second))
}
