#' Maximum type I error rate of a 2-in-1 adaptive design for count endpoints
#'
#' Evaluates the type I error rate under the global null hypothesis, in which
#' the two groups share the same rate for both endpoints, and maximises it over
#' a grid of plausible nuisance rates. The rates remain nuisance parameters
#' because they enter the information per patient and therefore the correlation
#' between the interim and the final statistics.
#'
#' @inheritParams adaptive2in1count
#' @param rate1 Numeric vector of plausible rates for the Phase 2 endpoint.
#' @param rate2 Numeric vector of plausible rates for the Phase 3 endpoint.
#' @return An object of class \code{adaptive2in1count_typeI} holding the data
#'   frame \code{results} with every combination of nuisance rates and
#'   cutpoints, the data frame \code{maximum} with the largest type I error rate
#'   for each cutpoint, and the overall maximum in \code{overall}.
#' @export
#' @examples
#' adaptive2in1count_typeI(rate1 = c(2, 2.6), rate2 = c(1, 1.3),
#'                         dispersion1 = 0.8, dispersion2 = 0.8, shared = 0.4,
#'                         n1 = 60, n2 = 90, n3 = 300, cutpoint = c(0.6, 1),
#'                         window1 = 0.5, window2 = 1, dropout = 0.3,
#'                         dropout_shape = "exponential")
adaptive2in1count_typeI <- function(rate1, rate2, dispersion1, dispersion2,
                                    shared, n1, n2, n3, cutpoint,
                                    window1 = 1, window2 = 1, dropout = 0,
                                    dropout_shape = c("complete", "uniform",
                                                      "exponential"),
                                    alpha = 0.025,
                                    small_endpoint = c("phase2", "phase3")) {
  dropout_shape <- match.arg(dropout_shape)
  small_endpoint <- match.arg(small_endpoint)
  if (!is.numeric(rate1) || length(rate1) == 0L || anyNA(rate1) ||
      any(rate1 <= 0)) {
    stop("'rate1' must be a numeric vector of positive values.")
  }
  if (!is.numeric(rate2) || length(rate2) == 0L || anyNA(rate2) ||
      any(rate2 <= 0)) {
    stop("'rate2' must be a numeric vector of positive values.")
  }
  grid <- expand.grid(rate1 = rate1, rate2 = rate2, KEEP.OUT.ATTRS = FALSE)
  pieces <- vector("list", nrow(grid))
  for (i in seq_len(nrow(grid))) {
    fit <- adaptive2in1count(rate1 = grid$rate1[i], rate2 = grid$rate2[i],
                             dispersion1 = dispersion1,
                             dispersion2 = dispersion2, shared = shared,
                             n1 = n1, n2 = n2, n3 = n3, cutpoint = cutpoint,
                             window1 = window1, window2 = window2,
                             dropout = dropout, dropout_shape = dropout_shape,
                             alpha = alpha, ratio1 = 1, ratio2 = 1,
                             small_endpoint = small_endpoint)
    pieces[[i]] <- data.frame(rate1 = grid$rate1[i], rate2 = grid$rate2[i],
                              cutpoint = fit$results$cutpoint,
                              go_probability = fit$results$go_probability,
                              type1 = fit$results$rejection,
                              expected_size = fit$results$expected_size)
  }
  results <- do.call(rbind, pieces)
  blocks <- split(results, results$cutpoint)
  maximum <- do.call(rbind, lapply(blocks, function(block) {
    block[which.max(block$type1), , drop = FALSE]
  }))
  maximum <- maximum[order(maximum$cutpoint), , drop = FALSE]
  rownames(maximum) <- NULL
  structure(list(results = results, maximum = maximum,
                 overall = maximum[which.max(maximum$type1), , drop = FALSE],
                 dispersion1 = dispersion1, dispersion2 = dispersion2,
                 shared = shared, n1 = expand_group_size(n1, "n1"),
                 n2 = expand_group_size(n2, "n2"),
                 n3 = expand_group_size(n3, "n3"), window1 = window1,
                 window2 = window2, dropout = dropout,
                 dropout_shape = dropout_shape, alpha = alpha,
                 small_endpoint = small_endpoint, call = match.call()),
            class = "adaptive2in1count_typeI")
}
