#' Rejection region of the Wald test for two negative binomial rates
#'
#' Evaluates the one-sided Wald test on the log rate ratio for every pair of
#' group totals, treating the dispersion as known. A continuity correction of
#' one half is applied to the totals so that the statistic remains finite when a
#' group total is zero. The region is provided as a comparator for
#' \code{exact_conditional_region()}, whose validity does not rely on the normal
#' approximation.
#'
#' @param size1 Number of patients in the investigational group.
#' @param size2 Number of patients in the control group.
#' @param dispersion Common dispersion of the negative binomial distribution.
#' @param smax Largest count retained on the support of each group total.
#' @param alpha One-sided significance level.
#' @return A numeric matrix of zeros and ones with \code{smax + 1} rows and
#'   \code{smax + 1} columns, indexed by the two group totals.
#' @importFrom stats qnorm
#' @export
#' @examples
#' region <- wald_region(size1 = 80, size2 = 80, dispersion = 0.5, smax = 60,
#'                       alpha = 0.025)
#' dim(region)
#' sum(region)
wald_region <- function(size1, size2, dispersion, smax, alpha = 0.025) {
  if (!is.numeric(size1) || length(size1) != 1L || is.na(size1) || size1 <= 0) {
    stop("'size1' must be a single positive numeric value.")
  }
  if (!is.numeric(size2) || length(size2) != 1L || is.na(size2) || size2 <= 0) {
    stop("'size2' must be a single positive numeric value.")
  }
  if (!is.numeric(dispersion) || length(dispersion) != 1L ||
      is.na(dispersion) || dispersion <= 0) {
    stop("'dispersion' must be a single positive numeric value.")
  }
  if (!is.numeric(smax) || length(smax) != 1L || is.na(smax) || smax < 0) {
    stop("'smax' must be a single non-negative integer.")
  }
  if (!is.numeric(alpha) || length(alpha) != 1L || is.na(alpha) ||
      alpha <= 0 || alpha >= 1) {
    stop("'alpha' must be a single numeric value between 0 and 1.")
  }
  wald_region_cpp(as.numeric(size1), as.numeric(size2), as.numeric(dispersion),
                  as.integer(smax), qnorm(1 - alpha))
}
