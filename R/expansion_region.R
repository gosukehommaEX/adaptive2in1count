#' Region of interim counts that leads to expansion
#'
#' Returns the indicator of the pair of interim counts for which the observed
#' rate ratio of the Phase 2 endpoint is at or below the cutpoint, so that the
#' trial expands into a large trial. The comparison is written with
#' multiplications only and carries a relative slack of \code{1e-12}, so that a
#' pair whose observed rate ratio equals the cutpoint exactly is assigned to the
#' expansion region regardless of the order in which the ratio is formed.
#'
#' @param amax Largest count retained on the support of the interim count.
#' @param n1 Number of patients per group at the interim analysis, either a
#'   single value or a vector of length two giving the investigational and the
#'   control group.
#' @param cutpoint Cutpoint for the observed rate ratio of the Phase 2 endpoint.
#' @return A logical matrix with \code{amax + 1} rows and \code{amax + 1}
#'   columns. Entry \code{[a1, a2]} is \code{TRUE} when interim counts
#'   \code{a1 - 1} and \code{a2 - 1} lead to expansion.
#' @export
#' @examples
#' region <- expansion_region(amax = 10, n1 = 20, cutpoint = 0.6)
#' dim(region)
#' region[2, 3]
#'
#' expansion_region(amax = 4, n1 = c(30, 10), cutpoint = 1)
expansion_region <- function(amax, n1, cutpoint) {
  if (!is.numeric(amax) || length(amax) != 1L || is.na(amax) || amax < 0) {
    stop("'amax' must be a single non-negative integer.")
  }
  if (!is.numeric(cutpoint) || length(cutpoint) != 1L || is.na(cutpoint) ||
      cutpoint <= 0) {
    stop("'cutpoint' must be a single positive numeric value.")
  }
  n1 <- expand_group_size(n1, "n1")
  counts <- seq(0, amax) + 0.5
  outer(counts * n1[2], cutpoint * (1 + 1e-12) * counts * n1[1], "<=")
}
