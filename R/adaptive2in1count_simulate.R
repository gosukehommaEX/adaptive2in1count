#' Simulated operating characteristics of a 2-in-1 adaptive design for counts
#'
#' Generates trials from the model with a random follow-up period and analyses
#' each of them by fitting a negative binomial model with an offset for the
#' realised exposure and applying a Wald test to the log rate ratio. The
#' function checks the normal approximation used by
#' \code{adaptive2in1count()}. With \code{endpoints = "distinct"} the counts
#' come from the additive gamma frailty model. With \code{endpoints =
#' "shared"} a single mixed Poisson process is generated per patient, the Phase
#' 2 count is the number of events in the shorter window and the Phase 3 count
#' adds the events of the remaining period.
#'
#' Three analyses are available. With \code{method = "estimated"} the two rates
#' and the common dispersion are estimated by maximum likelihood, which is the
#' analysis a trial would report. With \code{method = "known"} the dispersion is
#' held at the value assumed for the design, which is the situation the normal
#' approximation describes. Both fit every replication at once by Newton
#' iteration on matrices, so a large number of replications is inexpensive.
#' With \code{method = "glm"} the same model is fitted one replication at a time
#' by \code{MASS::glm.nb}; it gives the same answer far more slowly and is
#' provided so that the vectorised fit can be checked against a reference
#' implementation. Comparing \code{"known"} with \code{"estimated"} separates
#' the error of the normal approximation from the effect of estimating the
#' dispersion.
#'
#' @inheritParams adaptive2in1count
#' @param endpoints Relationship between the two endpoints, either
#'   \code{"distinct"} or \code{"shared"}. Under \code{"distinct"} the two
#'   endpoints count different events and are independent given the frailty.
#'   Under \code{"shared"} they count the same events over the two windows, so
#'   the Phase 2 count is part of the Phase 3 count. A single process then
#'   governs both endpoints, and the arguments \code{rate1},
#'   \code{dispersion1}, \code{shared} and \code{ratio1} are taken from
#'   \code{rate2}, \code{dispersion2} and \code{ratio2}.
#' @param method Analysis applied to each simulated trial, one of
#'   \code{"estimated"} for the maximum likelihood fit of the two rates and the
#'   dispersion, \code{"known"} for the same fit with the dispersion held at the
#'   design value, and \code{"glm"} for \code{MASS::glm.nb}.
#' @param nsim Number of simulated trials.
#' @param batch Number of replications generated at a time by the vectorised
#'   analyses, which bounds the memory used.
#' @param seed Optional seed passed to \code{set.seed()}.
#' @return A data frame with one row per cutpoint holding the simulated
#'   probability of expansion and rejection probability, their Monte Carlo
#'   standard errors, the number of usable replications and the number in which
#'   the analysis failed.
#' @importFrom stats qnorm rexp rgamma rpois runif vcov
#' @export
#' @examples
#' adaptive2in1count_simulate(rate1 = 2.6, rate2 = 1.3, dispersion1 = 0.8,
#'                            dispersion2 = 0.8, shared = 0.4, n1 = 60,
#'                            n2 = 90, n3 = 150, cutpoint = 0.8,
#'                            window1 = 0.5, window2 = 1, dropout = 0.3,
#'                            dropout_shape = "exponential", nsim = 500,
#'                            seed = 1)
adaptive2in1count_simulate <- function(rate1, rate2, dispersion1, dispersion2,
                                       shared, n1, n2, n3, cutpoint,
                                       window1 = 1, window2 = 1, dropout = 0,
                                       dropout_shape = c("complete", "uniform",
                                                         "exponential"),
                                       alpha = 0.025, ratio1 = 1, ratio2 = 1,
                                       small_endpoint = c("phase2", "phase3"),
                                       endpoints = c("distinct", "shared"),
                                       method = c("estimated", "known",
                                                  "glm"),
                                       nsim = 2000, batch = 2000,
                                       seed = NULL) {
  dropout_shape <- match.arg(dropout_shape)
  small_endpoint <- match.arg(small_endpoint)
  endpoints <- match.arg(endpoints)
  method <- match.arg(method)
  if (endpoints == "shared") {
    rate1 <- rate2
    dispersion1 <- dispersion2
    shared <- dispersion2
    ratio1 <- ratio2
  }
  if (!is.null(seed)) {
    set.seed(seed)
  }
  n1 <- expand_group_size(n1, "n1")
  n2 <- expand_group_size(n2, "n2")
  n3 <- expand_group_size(n3, "n3")
  if (any(n1 > n2) || any(n2 > n3)) {
    stop("'n1', 'n2' and 'n3' must be non-decreasing within each group.")
  }
  if (shared < 0 || shared > min(dispersion1, dispersion2)) {
    stop("'shared' must lie between 0 and the smaller of the two dispersions.")
  }
  if (length(dropout) == 1L) {
    dropout <- rep(dropout, 2L)
  }
  rate <- cbind(c(rate1 * ratio1, rate1), c(rate2 * ratio2, rate2))
  dispersion <- c(dispersion1, dispersion2)
  floor_exposure <- window2 * 1e-6
  critical <- qnorm(1 - alpha)
  small_column <- if (small_endpoint == "phase2") 1L else 2L
  draw_follow_up <- function(size, group) {
    if (dropout_shape == "complete") {
      return(rep(window2, size))
    }
    if (dropout_shape == "uniform") {
      # Two point mixture, as in follow_up_expectation(). A proportion
      # 'dropout' leave at a time uniform on the follow-up period and the rest
      # are followed to the end.
      early <- runif(size) < dropout[group]
      period <- ifelse(early, runif(size, 0, window2), window2)
      return(pmax(period, floor_exposure))
    }
    lambda <- -log(1 - dropout[group]) / window2
    pmax(pmin(rexp(size, rate = lambda), window2), floor_exposure)
  }
  # Patient level data for one group, with one row per replication.
  draw_group <- function(size, group) {
    follow_up <- matrix(draw_follow_up(size * n3[group], group), nrow = size)
    exposure1 <- pmin(follow_up, window1)
    if (endpoints == "shared") {
      # One process per patient. The Phase 2 count is the number of events in
      # the shorter window and the Phase 3 count adds an independent increment
      # over the remaining period, so the Phase 2 count is part of it.
      frailty <- matrix(rgamma(size * n3[group], shape = dispersion2,
                               rate = dispersion2), nrow = size)
      count1 <- matrix(rpois(size * n3[group], rate[group, 2] * exposure1 *
                               frailty), nrow = size)
      increment <- matrix(rpois(size * n3[group], rate[group, 2] *
                                 (follow_up - exposure1) * frailty),
                          nrow = size)
      count2 <- count1 + increment
      return(list(count1 = count1, count2 = count2, exposure1 = exposure1,
                  exposure2 = follow_up))
    }
    base <- matrix(rgamma(size * n3[group], shape = shared, rate = 1),
                   nrow = size)
    specific1 <- matrix(rgamma(size * n3[group], shape = dispersion1 - shared,
                               rate = 1), nrow = size)
    specific2 <- matrix(rgamma(size * n3[group], shape = dispersion2 - shared,
                               rate = 1), nrow = size)
    count1 <- matrix(rpois(size * n3[group], rate[group, 1] * exposure1 *
                             (base + specific1) / dispersion1), nrow = size)
    count2 <- matrix(rpois(size * n3[group], rate[group, 2] * follow_up *
                             (base + specific2) / dispersion2), nrow = size)
    list(count1 = count1, count2 = count2, exposure1 = exposure1,
         exposure2 = follow_up)
  }
  fit_glm <- function(count, exposure, group) {
    frame <- data.frame(count = count, offset_term = log(exposure),
                        arm = factor(group, levels = c(2, 1)))
    fit <- try(suppressWarnings(
      MASS::glm.nb(count ~ arm + offset(offset_term), data = frame)),
      silent = TRUE)
    if (inherits(fit, "try-error") ||
        !("arm1" %in% names(fit$coefficients))) {
      return(c(NA_real_, NA_real_))
    }
    variance <- try(vcov(fit)["arm1", "arm1"], silent = TRUE)
    if (inherits(variance, "try-error") || is.na(variance) ||
        variance <= 0) {
      return(c(NA_real_, NA_real_))
    }
    c(unname(fit$coefficients[["arm1"]]), sqrt(variance))
  }
  interim_estimate <- numeric(0)
  statistic_small <- numeric(0)
  statistic_large <- numeric(0)
  done <- 0L
  while (done < nsim) {
    size <- if (method == "glm") 1L else min(batch, nsim - done)
    parts <- lapply(seq_len(2), function(g) draw_group(size, g))
    if (method != "glm") {
      contrast <- function(sizes, endpoint) {
        count_field <- c("count1", "count2")[endpoint]
        exposure_field <- c("exposure1", "exposure2")[endpoint]
        piece <- function(g, field) {
          parts[[g]][[field]][, seq_len(sizes[g]), drop = FALSE]
        }
        negative_binomial_contrast(
          piece(1L, count_field), piece(1L, exposure_field),
          piece(2L, count_field), piece(2L, exposure_field),
          dispersion = if (method == "known") dispersion[endpoint] else NULL)
      }
      interim <- contrast(n1, 1L)
      small <- contrast(n2, small_column)
      large <- contrast(n3, 2L)
      interim_estimate <- c(interim_estimate, interim$estimate)
      statistic_small <- c(statistic_small, small$estimate / small$error)
      statistic_large <- c(statistic_large, large$estimate / large$error)
    } else {
      group <- rep(c(1, 2), times = n3)
      gather <- function(field) {
        unlist(lapply(seq_len(2), function(g) as.numeric(parts[[g]][[field]])))
      }
      count <- cbind(gather("count1"), gather("count2"))
      exposure <- cbind(gather("exposure1"), gather("exposure2"))
      index <- split(seq_along(group), group)
      take <- function(sizes) {
        unlist(lapply(seq_len(2), function(g) index[[g]][seq_len(sizes[g])]))
      }
      first <- take(n1)
      second <- take(n2)
      interim <- fit_glm(count[first, 1], exposure[first, 1], group[first])
      small <- fit_glm(count[second, small_column],
                       exposure[second, small_column], group[second])
      large <- fit_glm(count[, 2], exposure[, 2], group)
      interim_estimate <- c(interim_estimate, interim[1])
      statistic_small <- c(statistic_small, small[1] / small[2])
      statistic_large <- c(statistic_large, large[1] / large[2])
    }
    done <- done + size
  }
  usable <- is.finite(interim_estimate) & is.finite(statistic_small) &
    is.finite(statistic_large)
  failures <- sum(!usable)
  effective <- sum(usable)
  reject_small <- statistic_small[usable] <= -critical
  reject_large <- statistic_large[usable] <= -critical
  interim_value <- interim_estimate[usable]
  results <- data.frame(cutpoint = cutpoint, go_probability = NA_real_,
                        rejection = NA_real_, se_go = NA_real_,
                        se_rejection = NA_real_, nsim = effective,
                        failures = failures)
  for (i in seq_along(cutpoint)) {
    go <- interim_value <= log(cutpoint[i])
    reject <- ifelse(go, reject_large, reject_small)
    results$go_probability[i] <- mean(go)
    results$rejection[i] <- mean(reject)
    results$se_go[i] <- sqrt(mean(go) * (1 - mean(go)) / effective)
    results$se_rejection[i] <- sqrt(mean(reject) * (1 - mean(reject)) /
                                      effective)
  }
  results
}
