#' Print the correlation between the interim and the final test statistics
#'
#' Displays the within-patient endpoint correlation, the correlation of the
#' interim statistic with each final statistic, and whether the requirement for
#' type I error control is met.
#'
#' @param x An object of class \code{adaptive2in1count_correlation}.
#' @param digits Number of digits used when rounding the correlations.
#' @param ... Further arguments, which are ignored.
#' @return The object \code{x}, invisibly.
#' @export
#' @examples
#' rho <- adaptive2in1count_correlation(rate1 = 2.6, rate2 = 1.3,
#'                                      dispersion1 = 0.8, dispersion2 = 0.8,
#'                                      shared = 0.4, n1 = 60, n2 = 90,
#'                                      n3 = 300, window1 = 0.5, window2 = 1)
#' print(rho, digits = 5)
print.adaptive2in1count_correlation <- function(x, digits = 4, ...) {
  cat("Correlation structure of a 2-in-1 count design\n\n")
  cat(sprintf("Rates  Phase 2 %g, Phase 3 %g\n", x$rate1, x$rate2))
  cat(sprintf("Dispersions  Phase 2 %g, Phase 3 %g, shared frailty shape %g\n",
              x$dispersion1, x$dispersion2, x$shared))
  cat(sprintf("Windows  Phase 2 %g, Phase 3 %g\n", x$window1, x$window2))
  cat(sprintf("Follow-up  %s, dropout %s\n", x$dropout_shape,
              paste(format(x$dropout, digits = 3), collapse = "/")))
  cat(sprintf("Patients per group  interim %s, small trial %s, large trial %s\n\n",
              paste(x$n1, collapse = "/"), paste(x$n2, collapse = "/"),
              paste(x$n3, collapse = "/")))
  cat(sprintf("Within patient endpoint correlation      %s\n",
              format(round(x$endpoint, digits), nsmall = digits)))
  cat(sprintf("Interim with small trial statistic       %s\n",
              format(round(x$interim_small, digits), nsmall = digits)))
  cat(sprintf("Interim with large trial statistic       %s\n",
              format(round(x$interim_large, digits), nsmall = digits)))
  cat(sprintf("\nRequirement for type I error control is %s\n",
              if (isTRUE(x$holds)) "met" else "not met"))
  invisible(x)
}
