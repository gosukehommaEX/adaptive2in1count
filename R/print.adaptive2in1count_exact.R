#' Print the exact operating characteristics of a 2-in-1 count design
#'
#' Displays the design inputs and a table with one row per cutpoint holding the
#' probability of expansion, the rejection probability with its split into the
#' small trial and the large trial, and the expected total sample size.
#'
#' @param x An object of class \code{adaptive2in1count_exact}.
#' @param digits Number of digits used when rounding the probabilities.
#' @param ... Further arguments, which are ignored.
#' @return The object \code{x}, invisibly.
#' @export
#' @examples
#' fit <- adaptive2in1count_exact(rate1 = 3, rate2 = 0.6, dispersion1 = 0.5,
#'                          dispersion2 = 1.5, shared = 0.25,
#'                          n1 = 50, n2 = 80, n3 = 250,
#'                          cutpoint = c(0.5, 0.8))
#' print(fit, digits = 5)
print.adaptive2in1count_exact <- function(x, digits = 4, ...) {
  cat("2-in-1 count design, exact evaluation under complete follow-up\n\n")
  cat(sprintf("Phase 2 endpoint  control rate %g, dispersion %g,",
              x$rate1, x$dispersion1),
      sprintf(" true rate ratio %g\n", x$ratio1), sep = "")
  cat(sprintf("Phase 3 endpoint  control rate %g, dispersion %g,",
              x$rate2, x$dispersion2),
      sprintf(" true rate ratio %g\n", x$ratio2), sep = "")
  cat(sprintf("Shared frailty shape %g, endpoint correlation %s\n",
              x$shared,
              format(endpoint_correlation(x$rate1, x$rate2, x$dispersion1,
                                          x$dispersion2, x$shared),
                     digits = digits)))
  cat(sprintf("Patients per group  interim %s, small trial %s, large trial %s\n",
              paste(x$n1, collapse = "/"), paste(x$n2, collapse = "/"),
              paste(x$n3, collapse = "/")))
  cat(sprintf("Final analysis  %s test at one-sided alpha %g\n",
              if (x$test == "exact") "exact conditional" else "Wald", x$alpha))
  cat(sprintf("Small trial endpoint  %s\n\n",
              if (isTRUE(x$small_endpoint == "phase3")) "Phase 3" else "Phase 2"))
  out <- x$results
  out$go_probability <- round(out$go_probability, digits)
  out$rejection <- round(out$rejection, digits)
  out$rejection_small <- round(out$rejection_small, digits)
  out$rejection_large <- round(out$rejection_large, digits)
  out$expected_size <- round(out$expected_size, 1)
  names(out) <- c("cutpoint", "P(expand)", "rejection", "small", "large", "E(N)")
  print(out, row.names = FALSE)
  cat(sprintf("\nRetained supports  interim %d, small trial %d, large trial %d\n",
              x$support[["interim"]], x$support[["small"]], x$support[["large"]]))
  invisible(x)
}
