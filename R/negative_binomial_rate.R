#' Maximum likelihood rate and information for a single group
#'
#' Solves the score equation of the negative binomial model with an offset for
#' the exposure and a single rate parameter, for many replications at once. The
#' score for the log rate is
#' \code{sum(count) - sum((count + dispersion) * mean / (dispersion + mean))}
#' with \code{mean = rate * exposure}, and it is solved by Newton iteration from
#' the moment estimator. Convergence is quadratic, so a handful of iterations
#' reach machine precision. When every patient shares the same exposure the
#' moment estimator already solves the equation, whatever the dispersion.
#'
#' The information returned is the Fisher information for the log rate, namely
#' the sum over patients of \code{dispersion * mean / (dispersion + mean)}. Its
#' reciprocal is the asymptotic variance of the estimated log rate.
#'
#' @param count Numeric matrix of counts with one row per replication and one
#'   column per patient.
#' @param exposure Numeric matrix of exposures of the same shape as
#'   \code{count}.
#' @param dispersion Dispersion of the negative binomial distribution, either a
#'   single value or a vector with one element per replication.
#' @param iterations Largest number of Newton iterations.
#' @param tol Convergence tolerance on the relative change of the log rate.
#' @return A list with the estimated \code{rate} and the \code{information} for
#'   the log rate, both vectors with one element per replication.
#' @export
#' @examples
#' count <- matrix(c(2, 0, 3, 1, 1, 4), nrow = 2)
#' exposure <- matrix(1, nrow = 2, ncol = 3)
#' negative_binomial_rate(count, exposure, dispersion = 0.8)
negative_binomial_rate <- function(count, exposure, dispersion,
                                   iterations = 25L, tol = 1e-11) {
  if (!is.matrix(count) || !is.matrix(exposure) ||
      any(dim(count) != dim(exposure))) {
    stop("'count' and 'exposure' must be matrices of the same shape.")
  }
  if (!is.numeric(dispersion) || anyNA(dispersion) || any(dispersion <= 0) ||
      !(length(dispersion) %in% c(1L, nrow(count)))) {
    stop("'dispersion' must be positive, of length one or one per replication.")
  }
  total <- rowSums(count)
  rate <- total / rowSums(exposure)
  positive <- total > 0
  rate[!positive] <- 1
  # The two terms that do not change with the rate are formed once.
  inflated <- count + dispersion
  for (step in seq_len(iterations)) {
    mean_count <- rate * exposure
    ratio <- mean_count / (dispersion + mean_count)
    gradient <- total - rowSums(inflated * ratio)
    curvature <- rowSums(inflated * ratio * (1 - ratio))
    change <- gradient / curvature
    change[!is.finite(change)] <- 0
    change <- pmax(pmin(change, 1), -1)
    rate <- rate * exp(change)
    if (max(abs(change)) < tol) {
      break
    }
  }
  rate[!positive] <- 0
  mean_count <- rate * exposure
  information <- rowSums(dispersion * mean_count / (dispersion + mean_count))
  list(rate = rate, information = information)
}
