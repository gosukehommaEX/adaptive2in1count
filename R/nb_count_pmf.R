#' Negative binomial probability mass function on a truncated support
#'
#' Evaluates the negative binomial probability mass function with a given
#' dispersion and mean on the integers 0 to \code{ymax}. The probabilities are
#' accumulated on the log scale, so a dispersion multiplied by a large number of
#' patients does not overflow. A non-positive dispersion or mean returns a point
#' mass at zero, which is the degenerate limit used when the shared frailty
#' component or the number of additional patients is zero.
#'
#' @param dispersion Dispersion of the negative binomial distribution. Smaller
#'   values produce greater overdispersion.
#' @param mean Mean of the negative binomial distribution.
#' @param ymax Largest count retained on the support.
#' @return A numeric vector of length \code{ymax + 1} holding the probability of
#'   each count from 0 to \code{ymax}.
#' @export
#' @examples
#' pmf <- nb_count_pmf(dispersion = 0.5, mean = 3, ymax = 40)
#' sum(pmf)
#' sum(pmf * seq(0, 40))
nb_count_pmf <- function(dispersion, mean, ymax) {
  if (!is.numeric(dispersion) || length(dispersion) != 1L || is.na(dispersion)) {
    stop("'dispersion' must be a single numeric value.")
  }
  if (!is.numeric(mean) || length(mean) != 1L || is.na(mean)) {
    stop("'mean' must be a single numeric value.")
  }
  if (!is.numeric(ymax) || length(ymax) != 1L || is.na(ymax) || ymax < 0) {
    stop("'ymax' must be a single non-negative integer.")
  }
  nb_pmf_cpp(as.numeric(dispersion), as.numeric(mean), as.integer(ymax))
}
