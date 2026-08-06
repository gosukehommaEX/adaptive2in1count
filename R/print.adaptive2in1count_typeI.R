#' Print the maximum type I error rate of a 2-in-1 count design
#'
#' Displays the design inputs, the largest type I error rate attained at each
#' cutpoint over the grid of nuisance rates, and the overall maximum.
#'
#' @param x An object of class \code{adaptive2in1count_typeI}.
#' @param digits Number of digits used when rounding the probabilities.
#' @param ... Further arguments, which are ignored.
#' @return The object \code{x}, invisibly.
#' @export
#' @examples
#' size <- adaptive2in1count_typeI(rate1 = c(2, 3), rate2 = 0.6,
#'                                 dispersion1 = 0.5, dispersion2 = 1.5,
#'                                 shared = 0.25, n1 = 50, n2 = 80, n3 = 250,
#'                                 cutpoint = c(0.5, 1))
#' print(size, digits = 5)
print.adaptive2in1count_typeI <- function(x, digits = 5, ...) {
  cat("Maximum type I error rate of a 2-in-1 count design\n\n")
  cat(sprintf("Dispersions  Phase 2 %g, Phase 3 %g, shared frailty shape %g\n",
              x$dispersion1, x$dispersion2, x$shared))
  cat(sprintf("Patients per group  interim %s, small trial %s, large trial %s\n",
              paste(x$n1, collapse = "/"), paste(x$n2, collapse = "/"),
              paste(x$n3, collapse = "/")))
  cat(sprintf("Follow-up  %s, dropout %s\n", x$dropout_shape,
              paste(format(x$dropout, digits = 3), collapse = "/")))
  cat(sprintf("Wald test at one-sided alpha %g\n", x$alpha))
  cat(sprintf("Nuisance grid  %d combinations of the two control rates\n\n",
              nrow(x$results) / length(unique(x$results$cutpoint))))
  out <- x$maximum
  out$go_probability <- round(out$go_probability, digits)
  out$type1 <- round(out$type1, digits)
  out$expected_size <- round(out$expected_size, 1)
  names(out) <- c("rate1", "rate2", "cutpoint", "P(expand)", "type I", "E(N)")
  print(out, row.names = FALSE)
  cat(sprintf("\nOverall maximum type I error rate %s at cutpoint %g\n",
              format(round(x$overall$type1, digits), nsmall = digits),
              x$overall$cutpoint))
  invisible(x)
}
