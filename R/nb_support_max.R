#' Smallest support that retains all but a negligible upper tail
#'
#' Finds the smallest count at which the upper tail mass of a negative binomial
#' distribution falls below \code{tol}. Every distribution propagated by the
#' package is truncated at this count, so the reported operating characteristics
#' are exact up to the truncation tolerance and floating point precision.
#'
#' @param dispersion Dispersion of the negative binomial distribution.
#' @param mean Mean of the negative binomial distribution.
#' @param tol Upper tail mass that may be discarded.
#' @param cap Largest support that will be returned, which guards against a
#'   search that does not terminate for extreme parameter values.
#' @return An integer giving the largest count to retain.
#' @export
#' @examples
#' nb_support_max(dispersion = 0.5 * 80, mean = 3 * 80)
#' nb_support_max(dispersion = 1.5 * 250, mean = 0.6 * 250)
nb_support_max <- function(dispersion, mean, tol = 1e-11, cap = 20000L) {
  if (!is.numeric(dispersion) || length(dispersion) != 1L || is.na(dispersion)) {
    stop("'dispersion' must be a single numeric value.")
  }
  if (!is.numeric(mean) || length(mean) != 1L || is.na(mean)) {
    stop("'mean' must be a single numeric value.")
  }
  if (!is.numeric(tol) || length(tol) != 1L || is.na(tol) || tol <= 0) {
    stop("'tol' must be a single positive numeric value.")
  }
  if (dispersion <= 0 || mean <= 0) {
    return(0L)
  }
  cap <- as.integer(cap)
  spread <- sqrt(mean * (1 + mean / dispersion))
  hi <- min(cap, max(40L, as.integer(ceiling(mean + 15 * spread + 40))))
  repeat {
    pmf <- nb_pmf_cpp(as.numeric(dispersion), as.numeric(mean), as.integer(hi))
    upper <- 1 - cumsum(pmf)
    hit <- which(upper < tol)
    if (length(hit) > 0L) {
      return(as.integer(hit[1L] - 1L))
    }
    if (hi >= cap) {
      return(cap)
    }
    hi <- min(cap, as.integer(ceiling(hi * 1.6)))
  }
}
