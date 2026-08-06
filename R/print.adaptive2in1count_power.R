#' Print the power of a 2-in-1 count design
#'
#' Displays the design inputs and a table with one row per combination of true
#' rate ratios and cutpoints, holding the probability of expansion, the power
#' with its split into the small trial and the large trial, and the expected
#' total sample size.
#'
#' @param x An object of class \code{adaptive2in1count_power}.
#' @param digits Number of digits used when rounding the probabilities.
#' @param ... Further arguments, which are ignored.
#' @return The object \code{x}, invisibly.
#' @export
#' @examples
#' curve <- adaptive2in1count_power(rate1 = 3, rate2 = 0.6, dispersion1 = 0.5,
#'                                  dispersion2 = 1.5, shared = 0.25,
#'                                  n1 = 50, n2 = 80, n3 = 250,
#'                                  cutpoint = 0.5, ratio1 = 0.4,
#'                                  ratio2 = c(0.7, 1))
#' print(curve, digits = 3)
print.adaptive2in1count_power <- function(x, digits = 4, ...) {
  cat("Power of a 2-in-1 count design\n\n")
  cat(sprintf("Control rates  Phase 2 %g, Phase 3 %g\n", x$rate1, x$rate2))
  cat(sprintf("Dispersions  Phase 2 %g, Phase 3 %g, shared frailty shape %g\n",
              x$dispersion1, x$dispersion2, x$shared))
  cat(sprintf("Patients per group  interim %s, small trial %s, large trial %s\n",
              paste(x$n1, collapse = "/"), paste(x$n2, collapse = "/"),
              paste(x$n3, collapse = "/")))
  cat(sprintf("Follow-up  %s, dropout %s\n", x$dropout_shape,
              paste(format(x$dropout, digits = 3), collapse = "/")))
  cat(sprintf("Wald test at one-sided alpha %g\n\n", x$alpha))
  out <- x$results
  out$go_probability <- round(out$go_probability, digits)
  out$power <- round(out$power, digits)
  out$power_small <- round(out$power_small, digits)
  out$power_large <- round(out$power_large, digits)
  out$expected_size <- round(out$expected_size, 1)
  names(out) <- c("ratio1", "ratio2", "cutpoint", "P(expand)", "power",
                  "small", "large", "E(N)")
  print(out, row.names = FALSE)
  invisible(x)
}
