#' Print the operating characteristics of a 2-in-1 count design
#'
#' Displays the design inputs, the correlations that drive the interim decision
#' and a table with one row per cutpoint holding the probability of expansion,
#' the rejection probability with its split into the small trial and the large
#' trial, and the expected total sample size.
#'
#' @param x An object of class \code{adaptive2in1count}.
#' @param digits Number of digits used when rounding the probabilities.
#' @param ... Further arguments, which are ignored.
#' @return The object \code{x}, invisibly.
#' @export
#' @examples
#' fit <- adaptive2in1count(rate1 = 2.6, rate2 = 1.3, dispersion1 = 0.8,
#'                          dispersion2 = 0.8, shared = 0.4, n1 = 60, n2 = 90,
#'                          n3 = 300, cutpoint = c(0.6, 0.8),
#'                          window1 = 0.5, window2 = 1, dropout = 0.3,
#'                          dropout_shape = "exponential")
#' print(fit, digits = 5)
print.adaptive2in1count <- function(x, digits = 4, ...) {
  cat("2-in-1 adaptive design for count endpoints\n\n")
  cat(sprintf("Phase 2 endpoint  rate %g, dispersion %g, window %g,",
              x$rate1, x$dispersion1, x$window1),
      sprintf(" rate ratio %g\n", x$ratio1), sep = "")
  cat(sprintf("Phase 3 endpoint  rate %g, dispersion %g, window %g,",
              x$rate2, x$dispersion2, x$window2),
      sprintf(" rate ratio %g\n", x$ratio2), sep = "")
  cat(sprintf("Shared frailty shape %g\n", x$shared))
  cat(sprintf("Endpoints  %s\n",
              if (identical(x$endpoints, "shared")) {
                "same event over the two windows"
              } else {
                "different events"
              }))
  cat(sprintf("Follow-up  %s, dropout %s\n", x$dropout_shape,
              paste(format(x$dropout, digits = 3), collapse = "/")))
  cat(sprintf("Patients per group  interim %s, small trial %s, large trial %s\n",
              paste(x$n1, collapse = "/"), paste(x$n2, collapse = "/"),
              paste(x$n3, collapse = "/")))
  cat(sprintf("Small trial endpoint  %s\n",
              if (x$small_endpoint == "phase3") "Phase 3" else "Phase 2"))
  cat(sprintf("Wald test at one-sided alpha %g\n", x$alpha))
  cat(sprintf("Interim correlation  small trial %s, large trial %s\n\n",
              format(round(x$components$correlation[["interim_small"]], digits),
                     nsmall = digits),
              format(round(x$components$correlation[["interim_large"]], digits),
                     nsmall = digits)))
  out <- x$results
  out$go_probability <- round(out$go_probability, digits)
  out$rejection <- round(out$rejection, digits)
  out$rejection_small <- round(out$rejection_small, digits)
  out$rejection_large <- round(out$rejection_large, digits)
  out$expected_size <- round(out$expected_size, 1)
  names(out) <- c("cutpoint", "P(expand)", "rejection", "small", "large",
                  "E(N)")
  print(out, row.names = FALSE)
  invisible(x)
}
