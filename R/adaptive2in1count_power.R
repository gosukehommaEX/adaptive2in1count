#' Power of a 2-in-1 adaptive design for count endpoints
#'
#' Evaluates the probability of expansion, the power and the expected total
#' sample size over a grid of true rate ratios and a grid of cutpoints. Because
#' the endpoint tested at the final analysis depends on the interim decision,
#' the power is reported together with its split into the small trial and the
#' large trial.
#'
#' @inheritParams adaptive2in1count
#' @param ratio1 Numeric vector of true rate ratios for the Phase 2 endpoint.
#' @param ratio2 Numeric vector of true rate ratios for the Phase 3 endpoint.
#' @return An object of class \code{adaptive2in1count_power} holding the data
#'   frame \code{results} with one row per combination of rate ratios and
#'   cutpoints.
#' @export
#' @examples
#' adaptive2in1count_power(rate1 = 2.6, rate2 = 1.3, dispersion1 = 0.8,
#'                         dispersion2 = 0.8, shared = 0.4, n1 = 60, n2 = 90,
#'                         n3 = 300, cutpoint = 0.8, ratio1 = c(0.6, 0.7),
#'                         ratio2 = c(0.7, 0.8), window1 = 0.5, window2 = 1,
#'                         dropout = 0.3, dropout_shape = "exponential")
adaptive2in1count_power <- function(rate1, rate2, dispersion1, dispersion2,
                                    shared, n1, n2, n3, cutpoint, ratio1,
                                    ratio2, window1 = 1, window2 = 1,
                                    dropout = 0,
                                    dropout_shape = c("complete", "uniform",
                                                      "exponential"),
                                    alpha = 0.025,
                                    small_endpoint = c("phase2", "phase3")) {
  dropout_shape <- match.arg(dropout_shape)
  small_endpoint <- match.arg(small_endpoint)
  if (!is.numeric(ratio1) || length(ratio1) == 0L || anyNA(ratio1) ||
      any(ratio1 <= 0)) {
    stop("'ratio1' must be a numeric vector of positive values.")
  }
  if (!is.numeric(ratio2) || length(ratio2) == 0L || anyNA(ratio2) ||
      any(ratio2 <= 0)) {
    stop("'ratio2' must be a numeric vector of positive values.")
  }
  grid <- expand.grid(ratio1 = ratio1, ratio2 = ratio2, KEEP.OUT.ATTRS = FALSE)
  pieces <- vector("list", nrow(grid))
  for (i in seq_len(nrow(grid))) {
    fit <- adaptive2in1count(rate1 = rate1, rate2 = rate2,
                             dispersion1 = dispersion1,
                             dispersion2 = dispersion2, shared = shared,
                             n1 = n1, n2 = n2, n3 = n3, cutpoint = cutpoint,
                             window1 = window1, window2 = window2,
                             dropout = dropout, dropout_shape = dropout_shape,
                             alpha = alpha, ratio1 = grid$ratio1[i],
                             ratio2 = grid$ratio2[i],
                             small_endpoint = small_endpoint)
    pieces[[i]] <- data.frame(ratio1 = grid$ratio1[i], ratio2 = grid$ratio2[i],
                              cutpoint = fit$results$cutpoint,
                              go_probability = fit$results$go_probability,
                              power = fit$results$rejection,
                              power_small = fit$results$rejection_small,
                              power_large = fit$results$rejection_large,
                              expected_size = fit$results$expected_size)
  }
  results <- do.call(rbind, pieces)
  ordering <- order(results$cutpoint, results$ratio1, results$ratio2)
  results <- results[ordering, , drop = FALSE]
  rownames(results) <- NULL
  structure(list(results = results, rate1 = rate1, rate2 = rate2,
                 dispersion1 = dispersion1, dispersion2 = dispersion2,
                 shared = shared, n1 = expand_group_size(n1, "n1"),
                 n2 = expand_group_size(n2, "n2"),
                 n3 = expand_group_size(n3, "n3"), window1 = window1,
                 window2 = window2, dropout = dropout,
                 dropout_shape = dropout_shape, alpha = alpha,
                 small_endpoint = small_endpoint, call = match.call()),
            class = "adaptive2in1count_power")
}
