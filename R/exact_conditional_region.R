#' Rejection region of the exact conditional test for two negative binomial rates
#'
#' Under the null hypothesis of equal rates, a common exposure per patient and a
#' known common dispersion, the two group totals are negative binomial with a
#' common success probability. Conditioning on their sum removes that
#' probability and leaves a negative hypergeometric distribution, so the test is
#' exact and free of the nuisance rate. This is the count analogue of the
#' Fisher exact test for two binomial proportions.
#'
#' @param shape1 Dispersion multiplied by the number of patients in the
#'   investigational group.
#' @param shape2 Dispersion multiplied by the number of patients in the control
#'   group.
#' @param smax Largest count retained on the support of each group total.
#' @param alpha One-sided significance level.
#' @return A numeric matrix of zeros and ones with \code{smax + 1} rows and
#'   \code{smax + 1} columns. Entry \code{[s1, s2]} equals one when the test
#'   rejects the null hypothesis in favour of the investigational group having
#'   the lower rate, given totals \code{s1 - 1} and \code{s2 - 1}.
#' @export
#' @examples
#' region <- exact_conditional_region(shape1 = 40, shape2 = 40, smax = 60,
#'                                    alpha = 0.025)
#' dim(region)
#' sum(region)
exact_conditional_region <- function(shape1, shape2, smax, alpha = 0.025) {
  if (!is.numeric(shape1) || length(shape1) != 1L || is.na(shape1) || shape1 <= 0) {
    stop("'shape1' must be a single positive numeric value.")
  }
  if (!is.numeric(shape2) || length(shape2) != 1L || is.na(shape2) || shape2 <= 0) {
    stop("'shape2' must be a single positive numeric value.")
  }
  if (!is.numeric(smax) || length(smax) != 1L || is.na(smax) || smax < 0) {
    stop("'smax' must be a single non-negative integer.")
  }
  if (!is.numeric(alpha) || length(alpha) != 1L || is.na(alpha) ||
      alpha <= 0 || alpha >= 1) {
    stop("'alpha' must be a single numeric value between 0 and 1.")
  }
  exact_conditional_region_cpp(as.numeric(shape1), as.numeric(shape2),
                               as.integer(smax), as.numeric(alpha))
}
