#' Sample size of the large trial for a target power
#'
#' Solves for the number of patients per group in the large trial at which the
#' 2-in-1 design attains a target power, holding the interim and the small trial
#' sizes fixed. The search uses the normal approximation, so the returned value
#' should be checked with \code{adaptive2in1count()} after rounding.
#'
#' Enlarging the large trial only helps when the trial expands, so the power
#' approaches the probability of expansion plus the contribution of the small
#' trial and cannot be raised beyond it. A target above that limit is reported
#' as an error together with the limit, and a cutpoint that expands more often
#' is then the remedy.
#'
#' @inheritParams adaptive2in1count
#' @param target Target power.
#' @param allocation Allocation of the large trial between the investigational
#'   and the control group, a vector of length two that is scaled to the
#'   returned size.
#' @param interval Interval of per-group sizes searched for the solution. The
#'   search is a bisection over whole numbers, so the returned size is the
#'   smallest one in the interval that attains the target.
#' @return A list with the required per-group sizes in \code{n3}, the attained
#'   \code{power} and the expected total sample size in \code{expected_size}.
#' @export
#' @examples
#' adaptive2in1count_size(rate1 = 2.6, rate2 = 1.3, dispersion1 = 0.8,
#'                        dispersion2 = 0.8, shared = 0.4, n1 = 60, n2 = 90,
#'                        cutpoint = 1.4, target = 0.9, ratio1 = 0.65,
#'                        ratio2 = 0.75, window1 = 0.5, window2 = 1,
#'                        dropout = 0.3, dropout_shape = "exponential")
#'
#' # A cutpoint that rarely expands caps the power below the target, and the
#' # limit is reported instead of a size.
#' try(adaptive2in1count_size(rate1 = 2.6, rate2 = 1.3, dispersion1 = 0.8,
#'                            dispersion2 = 0.8, shared = 0.4, n1 = 60,
#'                            n2 = 90, cutpoint = 0.6, target = 0.9,
#'                            ratio1 = 0.65, ratio2 = 0.75, window1 = 0.5,
#'                            window2 = 1))
adaptive2in1count_size <- function(rate1, rate2, dispersion1, dispersion2,
                                   shared, n1, n2, cutpoint, target,
                                   ratio1, ratio2, window1 = 1, window2 = 1,
                                   dropout = 0,
                                   dropout_shape = c("complete", "uniform",
                                                     "exponential"),
                                   alpha = 0.025,
                                   small_endpoint = c("phase2", "phase3"),
                                   allocation = c(1, 1),
                                   interval = c(10, 20000)) {
  dropout_shape <- match.arg(dropout_shape)
  small_endpoint <- match.arg(small_endpoint)
  if (!is.numeric(target) || length(target) != 1L || is.na(target) ||
      target <= 0 || target >= 1) {
    stop("'target' must be a single numeric value between 0 and 1.")
  }
  if (length(cutpoint) != 1L) {
    stop("'cutpoint' must be a single value.")
  }
  allocation <- allocation / max(allocation)
  floor_size <- expand_group_size(n2, "n2")
  sizes_at <- function(size) pmax(round(size * allocation), floor_size)
  power_at <- function(size) {
    adaptive2in1count(rate1 = rate1, rate2 = rate2, dispersion1 = dispersion1,
                      dispersion2 = dispersion2, shared = shared, n1 = n1,
                      n2 = n2, n3 = sizes_at(size), cutpoint = cutpoint,
                      window1 = window1, window2 = window2, dropout = dropout,
                      dropout_shape = dropout_shape, alpha = alpha,
                      ratio1 = ratio1, ratio2 = ratio2,
                      small_endpoint = small_endpoint)$results
  }
  # Expanding is the only route to the large trial, so the power approaches the
  # probability of expansion plus the contribution of the small trial and
  # cannot be raised beyond it. The attainable limit is reported rather than
  # left to a search that would not terminate.
  lower <- max(interval[1], max(floor_size))
  upper <- max(interval[2], lower + 1)
  largest <- power_at(upper)
  if (largest$rejection < target) {
    stop(sprintf(paste("the target power is not attainable at this cutpoint;",
                       "the power approaches %.4f as the large trial grows,",
                       "because the probability of expansion is %.4f."),
                 largest$rejection, largest$go_probability))
  }
  if (power_at(lower)$rejection < target) {
    while (upper - lower > 1) {
      middle <- floor((lower + upper) / 2)
      if (power_at(middle)$rejection >= target) {
        upper <- middle
      } else {
        lower <- middle
      }
    }
  } else {
    upper <- lower
  }
  size <- sizes_at(upper)
  final <- adaptive2in1count(rate1 = rate1, rate2 = rate2,
                             dispersion1 = dispersion1,
                             dispersion2 = dispersion2, shared = shared,
                             n1 = n1, n2 = n2, n3 = size, cutpoint = cutpoint,
                             window1 = window1, window2 = window2,
                             dropout = dropout, dropout_shape = dropout_shape,
                             alpha = alpha, ratio1 = ratio1, ratio2 = ratio2,
                             small_endpoint = small_endpoint)
  list(n3 = size, power = final$results$rejection,
       expected_size = final$results$expected_size)
}
