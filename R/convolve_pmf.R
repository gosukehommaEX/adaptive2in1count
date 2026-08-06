#' Convolution of two probability mass functions on a truncated support
#'
#' Returns the probability mass function of the sum of two independent
#' non-negative integer valued variables, truncated at \code{ymax}. The
#' truncation discards only the upper tail that was already outside the retained
#' support of the two inputs.
#'
#' @param f Numeric vector holding the first probability mass function.
#' @param g Numeric vector holding the second probability mass function.
#' @param ymax Largest count retained on the support of the result.
#' @return A numeric vector of length \code{ymax + 1}.
#' @export
#' @examples
#' f <- nb_count_pmf(dispersion = 10, mean = 5, ymax = 60)
#' g <- nb_count_pmf(dispersion = 20, mean = 8, ymax = 60)
#' h <- convolve_pmf(f, g, ymax = 60)
#' sum(h * seq(0, 60))
convolve_pmf <- function(f, g, ymax) {
  if (!is.numeric(f) || length(f) == 0L) {
    stop("'f' must be a non-empty numeric vector.")
  }
  if (!is.numeric(g) || length(g) == 0L) {
    stop("'g' must be a non-empty numeric vector.")
  }
  if (!is.numeric(ymax) || length(ymax) != 1L || is.na(ymax) || ymax < 0) {
    stop("'ymax' must be a single non-negative integer.")
  }
  out <- convolve_rows_cpp(matrix(as.numeric(f), nrow = 1L), as.numeric(g),
                           as.integer(ymax))
  as.numeric(out[1L, ])
}
