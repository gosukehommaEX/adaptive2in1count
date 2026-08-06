#' Estimated log rate ratio of two groups under a negative binomial model
#'
#' Fits the negative binomial model with a common dispersion, a separate rate in
#' each group and an offset for the exposure, and returns the estimated log rate
#' ratio with its standard error. Many replications are fitted at once by
#' Newton iteration, so no loop over replications is needed.
#'
#' When \code{dispersion} is supplied it is held fixed. When it is \code{NULL}
#' the dispersion is estimated jointly with the two rates by alternating a
#' Newton step on each log rate with a Newton step on the log dispersion, which
#' is the maximum likelihood fit that \code{MASS::glm.nb} produces for the same
#' model. The standard error uses the expected information for the log rates,
#' as \code{MASS::glm.nb} does, so the uncertainty of the dispersion is not
#' propagated.
#'
#' @param count1 Numeric matrix of counts in the investigational group, with one
#'   row per replication.
#' @param exposure1 Numeric matrix of exposures in the investigational group.
#' @param count2 Numeric matrix of counts in the control group.
#' @param exposure2 Numeric matrix of exposures in the control group.
#' @param dispersion Common dispersion, or \code{NULL} to estimate it.
#' @param iterations Largest number of alternating iterations when the
#'   dispersion is estimated.
#' @param tol Convergence tolerance on the relative change of the parameters.
#' @return A list with the estimated log rate ratio in \code{estimate}, its
#'   standard error in \code{error}, and the dispersion used or estimated in
#'   \code{dispersion}. Replications in which a group observed no event give
#'   \code{NA}.
#' @export
#' @examples
#' set.seed(1)
#' exposure <- matrix(1, nrow = 5, ncol = 60)
#' count1 <- matrix(rnbinom(300, size = 0.8, mu = 0.9), nrow = 5)
#' count2 <- matrix(rnbinom(300, size = 0.8, mu = 1.3), nrow = 5)
#' negative_binomial_contrast(count1, exposure, count2, exposure)
negative_binomial_contrast <- function(count1, exposure1, count2, exposure2,
                                       dispersion = NULL, iterations = 25L,
                                       tol = 1e-8) {
  replications <- nrow(count1)
  total1 <- rowSums(count1)
  total2 <- rowSums(count2)
  usable <- total1 > 0 & total2 > 0
  if (is.null(dispersion)) {
    value <- rep(1, replications)
    inflated1 <- count1
    inflated2 <- count2
    rate1 <- total1 / rowSums(exposure1)
    rate2 <- total2 / rowSums(exposure2)
    rate1[!usable] <- 1
    rate2[!usable] <- 1
    for (step in seq_len(iterations)) {
      previous <- c(rate1, rate2, value)
      {
        for (side in seq_len(2)) {
          count <- if (side == 1L) count1 else count2
          exposure <- if (side == 1L) exposure1 else exposure2
          total <- if (side == 1L) total1 else total2
          rate <- if (side == 1L) rate1 else rate2
          mean_count <- rate * exposure
          ratio <- mean_count / (value + mean_count)
          gradient <- total - rowSums((count + value) * ratio)
          curvature <- rowSums((count + value) * ratio * (1 - ratio))
          change <- gradient / curvature
          change[!is.finite(change)] <- 0
          rate <- rate * exp(pmax(pmin(change, 1), -1))
          if (side == 1L) rate1 <- rate else rate2 <- rate
        }
      }
      gradient <- 0
      curvature <- 0
      for (side in seq_len(2)) {
        count <- if (side == 1L) count1 else count2
        exposure <- if (side == 1L) exposure1 else exposure2
        rate <- if (side == 1L) rate1 else rate2
        mean_count <- rate * exposure
        sum_term <- value + mean_count
        gradient <- gradient +
          rowSums(digamma(count + value) - digamma(value) + log(value) -
                    log(sum_term) + 1 - (count + value) / sum_term)
        curvature <- curvature +
          rowSums(trigamma(count + value) - trigamma(value) + 1 / value -
                    2 / sum_term + (count + value) / sum_term^2)
      }
      # The Newton step is taken on the log dispersion, where the parameter is
      # unconstrained, and is capped so that a flat likelihood cannot send the
      # dispersion to an extreme in a single step.
      slope <- value * gradient
      bend <- value^2 * curvature + slope
      change <- ifelse(bend < 0, -slope / bend, sign(slope) * 0.5)
      change[!is.finite(change)] <- 0
      value <- value * exp(pmax(pmin(change, 0.5), -0.5))
      value <- pmin(pmax(value, 1e-4), 1e6)
      if (max(abs(c(rate1, rate2, value) - previous) /
                pmax(abs(previous), 1e-8)) < tol) {
        break
      }
    }
    dispersion <- value
  } else {
    fit1 <- negative_binomial_rate(count1, exposure1, dispersion)
    fit2 <- negative_binomial_rate(count2, exposure2, dispersion)
    rate1 <- fit1$rate
    rate2 <- fit2$rate
  }
  mean1 <- rate1 * exposure1
  mean2 <- rate2 * exposure2
  information1 <- rowSums(dispersion * mean1 / (dispersion + mean1))
  information2 <- rowSums(dispersion * mean2 / (dispersion + mean2))
  estimate <- rep(NA_real_, replications)
  error <- rep(NA_real_, replications)
  estimate[usable] <- log(rate1[usable]) - log(rate2[usable])
  error[usable] <- sqrt(1 / information1[usable] + 1 / information2[usable])
  list(estimate = estimate, error = error,
       dispersion = if (length(dispersion) == 1L) {
         rep(dispersion, replications)
       } else {
         dispersion
       })
}
