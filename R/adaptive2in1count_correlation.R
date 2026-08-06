#' Correlation between the interim and the final test statistics
#'
#' Type I error control of the 2-in-1 design rests on the requirement that the
#' interim statistic be at least as strongly correlated with the final statistic
#' of the small trial as with the final statistic of the large trial. This
#' function evaluates the two correlations under the normal approximation and
#' reports whether the requirement is met.
#'
#' @inheritParams design_components
#' @return An object of class \code{adaptive2in1count_correlation} holding the
#'   within-patient endpoint correlation, the correlation of the interim
#'   statistic with each final statistic, and whether the requirement is met.
#' @export
#' @examples
#' adaptive2in1count_correlation(rate1 = 2.6, rate2 = 1.3, dispersion1 = 0.8,
#'                               dispersion2 = 0.8, shared = 0.4, n1 = 60,
#'                               n2 = 90, n3 = 300, window1 = 0.5, window2 = 1,
#'                               dropout = 0.3,
#'                               dropout_shape = "exponential")
adaptive2in1count_correlation <- function(rate1, rate2, dispersion1,
                                          dispersion2, shared, n1, n2, n3,
                                          window1 = 1, window2 = 1,
                                          dropout = 0,
                                          dropout_shape = c("complete",
                                                            "uniform",
                                                            "exponential"),
                                          ratio1 = 1, ratio2 = 1,
                                          small_endpoint = c("phase2",
                                                             "phase3")) {
  dropout_shape <- match.arg(dropout_shape)
  small_endpoint <- match.arg(small_endpoint)
  parts <- design_components(rate1 = rate1, rate2 = rate2,
                             dispersion1 = dispersion1,
                             dispersion2 = dispersion2, shared = shared,
                             n1 = n1, n2 = n2, n3 = n3, window1 = window1,
                             window2 = window2, dropout = dropout,
                             dropout_shape = dropout_shape, ratio1 = ratio1,
                             ratio2 = ratio2, small_endpoint = small_endpoint)
  endpoint <- endpoint_correlation(rate1, rate2, dispersion1, dispersion2,
                                   shared, window1, window2, parts$dropout[2],
                                   dropout_shape)
  structure(list(endpoint = endpoint,
                 interim_small = parts$correlation[["interim_small"]],
                 interim_large = parts$correlation[["interim_large"]],
                 holds = parts$correlation[["interim_small"]] >=
                   parts$correlation[["interim_large"]],
                 coupling = parts$coupling, rate1 = rate1, rate2 = rate2,
                 dispersion1 = dispersion1, dispersion2 = dispersion2,
                 shared = shared, n1 = expand_group_size(n1, "n1"),
                 n2 = expand_group_size(n2, "n2"),
                 n3 = expand_group_size(n3, "n3"), window1 = window1,
                 window2 = window2, dropout = parts$dropout,
                 dropout_shape = dropout_shape,
                 small_endpoint = small_endpoint, call = match.call()),
            class = "adaptive2in1count_correlation")
}
