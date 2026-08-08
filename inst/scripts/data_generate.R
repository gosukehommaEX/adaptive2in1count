# Generates the operating characteristics shown by the figures and tables of the
# main text. Run from the package root. Results are written as .rds files to
# inst/scripts/data. The exhibits of the Supporting Information are produced by
# the pair of scripts in inst/supplementary.

library(adaptive2in1count)

data_dir <- file.path("inst", "scripts", "data")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

# Reference design in chronic obstructive pulmonary disease. The Phase 2
# endpoint is the number of moderate or severe exacerbations over half a year
# and the Phase 3 endpoint is the annualized exacerbation rate over a year, with
# rates and dispersions in the range reported for the TRISTAN study by Keene et
# al. (2007). That study lost 39 per cent of the placebo group before the end of
# the trial, so a dropout proportion of 0.3 is taken as the reference.
rate1 <- 2
rate2 <- 1.3
dispersion1 <- 0.7
dispersion2 <- 0.8
shared <- 0.4
window1 <- 0.5
window2 <- 1
n1 <- 80
n2 <- 120
n3 <- 750
alpha <- 0.025
reference_dropout <- 0.3

# Design alternative used whenever a single alternative is required.
alt_ratio1 <- 0.65
alt_ratio2 <- 0.75
target_power <- 0.85

# The exponential follow-up period is not defined at a dropout proportion of
# zero, where the period is complete by definition.
shape_for <- function(dropout) {
  if (dropout <= 0) "complete" else "exponential"
}

cutpoints <- seq(0.4, 2, by = 0.05)

## (1) Maximum type I error rate against the cutpoint.
message("(1/9) maximum type I error rate against the cutpoint")
nuisance_rate1 <- c(1.2, 2, 3.2)
nuisance_rate2 <- c(0.8, 1.3, 2)
panel1 <- expand.grid(dispersion2 = c(0.5, 0.8, 1.5),
                      dropout = c(0, 0.2, 0.4),
                      endpoint = c("phase2", "phase3"),
                      stringsAsFactors = FALSE)
typeI_curve <- do.call(rbind, lapply(seq_len(nrow(panel1)), function(i) {
  fit <- adaptive2in1count_typeI(
    rate1 = nuisance_rate1, rate2 = nuisance_rate2,
    dispersion1 = dispersion1, dispersion2 = panel1$dispersion2[i],
    shared = min(shared, dispersion1, panel1$dispersion2[i]),
    n1 = n1, n2 = n2, n3 = n3, cutpoint = cutpoints, window1 = window1,
    window2 = window2, dropout = panel1$dropout[i],
    dropout_shape = shape_for(panel1$dropout[i]), alpha = alpha,
    small_endpoint = panel1$endpoint[i])
  data.frame(dispersion2 = panel1$dispersion2[i], dropout = panel1$dropout[i],
             endpoint = panel1$endpoint[i], fit$maximum,
             stringsAsFactors = FALSE)
}))
saveRDS(typeI_curve, file.path(data_dir, "typeI_curve.rds"))

## (2) Power against the cutpoint, by treatment effect and dropout.
message("(2/9) power against the cutpoint")
panel2 <- expand.grid(ratio2 = c(0.7, 0.8), ratio1 = c(0.55, 0.65, 0.75),
                      dropout = c(0, 0.2, 0.4))
power_curve <- do.call(rbind, lapply(seq_len(nrow(panel2)), function(i) {
  fit <- adaptive2in1count_power(
    rate1 = rate1, rate2 = rate2, dispersion1 = dispersion1,
    dispersion2 = dispersion2, shared = shared, n1 = n1, n2 = n2, n3 = n3,
    cutpoint = cutpoints, ratio1 = panel2$ratio1[i], ratio2 = panel2$ratio2[i],
    window1 = window1, window2 = window2, dropout = panel2$dropout[i],
    dropout_shape = shape_for(panel2$dropout[i]), alpha = alpha)
  data.frame(dropout = panel2$dropout[i], fit$results)
}))
saveRDS(power_curve, file.path(data_dir, "power_curve.rds"))

## (3) Within-patient correlation of the two endpoints, by its two sources.
message("(3/9) within-patient correlation of the two endpoints")
panel4 <- expand.grid(dropout = seq(0, 0.6, by = 0.05),
                      shared = c(0, 0.23, 0.47, 0.7),
                      window1 = c(0.25, 0.5, 1),
                      dispersion2 = c(0.8, 1.5))
correlation_endpoint <- do.call(rbind, lapply(seq_len(nrow(panel4)), function(i) {
  data.frame(dropout = panel4$dropout[i], shared = panel4$shared[i],
             window1 = panel4$window1[i], dispersion2 = panel4$dispersion2[i],
             correlation = endpoint_correlation(
               rate1 = rate1, rate2 = rate2, dispersion1 = dispersion1,
               dispersion2 = panel4$dispersion2[i], shared = panel4$shared[i],
               window1 = panel4$window1[i], window2 = window2,
               dropout = panel4$dropout[i],
               dropout_shape = shape_for(panel4$dropout[i])))
}))
saveRDS(correlation_endpoint, file.path(data_dir, "correlation_endpoint.rds"))

## (4) Correlation of the interim statistic with the two final statistics when
## the two endpoints count the same events over the two windows.
message("(4/9) correlation between the statistics for a shared endpoint")
panel9 <- expand.grid(n3 = seq(120, 900, by = 20), dropout = c(0, 0.2, 0.4),
                      window1 = c(0.5, 1))
correlation_shared <- do.call(rbind, lapply(seq_len(nrow(panel9)), function(i) {
  common <- list(rate1 = rate2, rate2 = rate2, dispersion1 = dispersion2,
                 dispersion2 = dispersion2, shared = dispersion2, n1 = n1,
                 n2 = n2, n3 = panel9$n3[i], window1 = panel9$window1[i],
                 window2 = window2, dropout = panel9$dropout[i],
                 dropout_shape = shape_for(panel9$dropout[i]))
  rho <- do.call(adaptive2in1count_correlation,
                 c(common, list(endpoints = "shared")))
  # The distinct structure at the same rates and windows, so that the two
  # structures can be read off the same panel.
  apart <- adaptive2in1count_correlation(
    rate1 = rate1, rate2 = rate2, dispersion1 = dispersion1,
    dispersion2 = dispersion2, shared = shared, n1 = n1, n2 = n2,
    n3 = panel9$n3[i], window1 = panel9$window1[i], window2 = window2,
    dropout = panel9$dropout[i],
    dropout_shape = shape_for(panel9$dropout[i]))
  data.frame(n3 = panel9$n3[i], dropout = panel9$dropout[i],
             window1 = panel9$window1[i],
             shared_small = rho$interim_small,
             shared_large = rho$interim_large,
             distinct_small = apart$interim_small,
             distinct_large = apart$interim_large)
}))
saveRDS(correlation_shared, file.path(data_dir, "correlation_shared.rds"))

## (5) Operating characteristics at cutpoints calibrated to a target
## probability of expansion.
message("(5/9) operating characteristics at calibrated cutpoints")
targets <- seq(0.3, 0.95, by = 0.05)
calibrated <- vapply(targets, function(target) {
  adaptive2in1count_cutpoint(
    rate1 = rate1, rate2 = rate2, dispersion1 = dispersion1,
    dispersion2 = dispersion2, shared = shared, n1 = n1, n2 = n2, n3 = n3,
    target = target, window1 = window1, window2 = window2,
    dropout = reference_dropout, dropout_shape = "exponential",
    ratio1 = alt_ratio1)$cutpoint
}, numeric(1))
common <- list(rate1 = rate1, rate2 = rate2, dispersion1 = dispersion1,
               dispersion2 = dispersion2, shared = shared, n1 = n1, n2 = n2,
               n3 = n3, cutpoint = calibrated, window1 = window1,
               window2 = window2, dropout = reference_dropout,
               dropout_shape = "exponential", alpha = alpha)
under_null <- do.call(adaptive2in1count, common)$results
under_alternative <- do.call(adaptive2in1count,
                             c(common, list(ratio1 = alt_ratio1,
                                            ratio2 = alt_ratio2)))$results
operating_table <- data.frame(
  target = targets, cutpoint = calibrated,
  go_null = under_null$go_probability, type1 = under_null$rejection,
  size_null = under_null$expected_size,
  go_alternative = under_alternative$go_probability,
  power = under_alternative$rejection,
  power_small = under_alternative$rejection_small,
  power_large = under_alternative$rejection_large,
  size_alternative = under_alternative$expected_size)
saveRDS(operating_table, file.path(data_dir, "operating_table.rds"))

## (6) The ceiling that the cutpoint places on the power.
message("(6/9) the ceiling that the cutpoint places on the power")
ceiling_targets <- seq(0.3, 0.95, by = 0.05)
ceiling_size <- 200000
ceiling_table <- do.call(rbind, lapply(ceiling_targets, function(target) {
  cutpoint <- adaptive2in1count_cutpoint(
    rate1 = rate1, rate2 = rate2, dispersion1 = dispersion1,
    dispersion2 = dispersion2, shared = shared, n1 = n1, n2 = n2, n3 = n3,
    target = target, window1 = window1, window2 = window2,
    dropout = reference_dropout, dropout_shape = "exponential",
    ratio1 = alt_ratio1)$cutpoint
  arguments <- list(rate1 = rate1, rate2 = rate2, dispersion1 = dispersion1,
                    dispersion2 = dispersion2, shared = shared, n1 = n1,
                    n2 = n2, cutpoint = cutpoint, window1 = window1,
                    window2 = window2, dropout = reference_dropout,
                    dropout_shape = "exponential", alpha = alpha,
                    ratio1 = alt_ratio1, ratio2 = alt_ratio2)
  attained <- do.call(adaptive2in1count, c(arguments, list(n3 = n3)))$results
  # Enlarging the large trial cannot raise the power beyond the probability of
  # expansion plus the contribution of the small trial, so the limit is read off
  # at a size far larger than any trial would use.
  limit <- do.call(adaptive2in1count,
                   c(arguments, list(n3 = ceiling_size)))$results
  solved <- try(do.call(adaptive2in1count_size,
                        c(arguments, list(target = target_power))),
                silent = TRUE)
  data.frame(target = target, cutpoint = cutpoint,
             go_probability = attained$go_probability,
             power = attained$rejection, power_small = attained$rejection_small,
             ceiling = limit$rejection,
             required = if (inherits(solved, "try-error")) {
               NA_real_
             } else {
               sum(solved$n3)
             })
}))
saveRDS(ceiling_table, file.path(data_dir, "ceiling_table.rds"))

## (7) The normal approximation against simulation.
message("(7/9) the normal approximation against simulation")
validation_cutpoints <- c(0.8, 1.2, 1.6)
validation_scenarios <- list(
  list(label = "Reference design, global null", ratio1 = 1, ratio2 = 1,
       dropout = 0, sizes = c(n1, n2, n3), nsim = 10000),
  list(label = "Reference design, global null, dropout 0.3", ratio1 = 1,
       ratio2 = 1, dropout = reference_dropout, sizes = c(n1, n2, n3),
       nsim = 10000),
  list(label = "Reference design, both endpoints, dropout 0.3",
       ratio1 = alt_ratio1, ratio2 = alt_ratio2, dropout = reference_dropout,
       sizes = c(n1, n2, n3), nsim = 10000),
  list(label = "Reference design, both endpoints, dropout 0.5",
       ratio1 = alt_ratio1, ratio2 = alt_ratio2, dropout = 0.5,
       sizes = c(n1, n2, n3), nsim = 10000),
  list(label = "Small design, global null, dropout 0.3", ratio1 = 1,
       ratio2 = 1, dropout = reference_dropout, sizes = c(30, 45, 200),
       nsim = 50000),
  list(label = "Small design, both endpoints, dropout 0.3",
       ratio1 = alt_ratio1, ratio2 = alt_ratio2, dropout = reference_dropout,
       sizes = c(30, 45, 200), nsim = 50000))
validation_table <- do.call(rbind, lapply(validation_scenarios, function(z) {
  arguments <- list(rate1 = rate1, rate2 = rate2, dispersion1 = dispersion1,
                    dispersion2 = dispersion2, shared = shared,
                    n1 = z$sizes[1], n2 = z$sizes[2], n3 = z$sizes[3],
                    cutpoint = validation_cutpoints, window1 = window1,
                    window2 = window2, dropout = z$dropout,
                    dropout_shape = shape_for(z$dropout), alpha = alpha,
                    ratio1 = z$ratio1, ratio2 = z$ratio2)
  approximate <- do.call(adaptive2in1count, arguments)$results
  known <- do.call(adaptive2in1count_simulate,
                   c(arguments, list(method = "known", nsim = z$nsim,
                                     seed = 20260806)))
  estimated <- do.call(adaptive2in1count_simulate,
                       c(arguments, list(method = "estimated", nsim = z$nsim,
                                         seed = 20260806)))
  data.frame(label = z$label, cutpoint = validation_cutpoints,
             approximate = approximate$rejection, known = known$rejection,
             estimated = estimated$rejection, se_known = known$se_rejection,
             se_estimated = estimated$se_rejection, nsim = z$nsim,
             stringsAsFactors = FALSE)
}))
saveRDS(validation_table, file.path(data_dir, "validation_table.rds"))

## (8) The normal approximation against simulation for a shared endpoint.
message("(8/9) the normal approximation for a shared endpoint")
shared_scenarios <- list(
  list(label = "Shared endpoint, global null, dropout 0.3", ratio = 1,
       dropout = reference_dropout, window1 = 0.5),
  list(label = "Shared endpoint, both windows, dropout 0.3", ratio = 0.75,
       dropout = reference_dropout, window1 = 0.5),
  list(label = "Shared endpoint, one window, dropout 0.3", ratio = 0.75,
       dropout = reference_dropout, window1 = window2))
validation_shared <- do.call(rbind, lapply(shared_scenarios, function(z) {
  arguments <- list(rate1 = rate2, rate2 = rate2, dispersion1 = dispersion2,
                    dispersion2 = dispersion2, shared = dispersion2, n1 = n1,
                    n2 = n2, n3 = n3, cutpoint = validation_cutpoints,
                    window1 = z$window1, window2 = window2,
                    dropout = z$dropout, dropout_shape = shape_for(z$dropout),
                    alpha = alpha, ratio1 = z$ratio, ratio2 = z$ratio,
                    endpoints = "shared")
  approximate <- do.call(adaptive2in1count, arguments)$results
  known <- do.call(adaptive2in1count_simulate,
                   c(arguments, list(method = "known", nsim = 10000,
                                     seed = 20260806)))
  estimated <- do.call(adaptive2in1count_simulate,
                       c(arguments, list(method = "estimated", nsim = 10000,
                                         seed = 20260806)))
  data.frame(label = z$label, cutpoint = validation_cutpoints,
             approximate = approximate$rejection, known = known$rejection,
             estimated = estimated$rejection, se_known = known$se_rejection,
             se_estimated = estimated$se_rejection, nsim = 10000,
             stringsAsFactors = FALSE)
}))
validation_table <- rbind(validation_table, validation_shared)
saveRDS(validation_table, file.path(data_dir, "validation_table.rds"))

## (9) The tezepelumab development programme in severe asthma, evaluated as a
## single 2-in-1 trial. Every parameter is a design assumption taken from the
## two trials rather than a result of either of them. The protocols state the
## dispersion as the variance of the frailty, in the sense of Tang (2015) and
## Zhu and Lakkis (2014), so it enters here as its reciprocal. The interim
## analysis reads the endpoint on the enrolled patients before their year is
## complete and both final analyses read it at the year, so the small trial and
## the large trial share an endpoint and the interim does not.
message("(9/9) the tezepelumab design")
real_settings <- list(
  list(label = "Phase 3 assumptions", rate = 0.9, dispersion = 1 / 2.4,
       ratio = 0.5),
  list(label = "Phase 2 assumptions", rate = 0.7, dispersion = 1 / 0.7,
       ratio = 0.6))
real_n1 <- 138
real_n2 <- 138
real_n3 <- 530
real_dropout <- 0.1
real_alpha <- 0.025
real_power <- 0.9
# The interim analysis reads the endpoint before the final analysis does, so
# the interim window stays strictly below the year.
real_windows <- c(0.25, 0.5, 0.75)
real_targets <- seq(0.3, 0.95, by = 0.05)

real_arguments <- function(z, window1, dropout, n3, ratio = z$ratio) {
  list(rate1 = z$rate, rate2 = z$rate, dispersion1 = z$dispersion,
       dispersion2 = z$dispersion, shared = z$dispersion, n1 = real_n1,
       n2 = real_n2, n3 = n3, window1 = window1, window2 = 1,
       dropout = dropout, dropout_shape = "uniform", alpha = real_alpha,
       ratio1 = ratio, ratio2 = ratio, small_endpoint = "phase3",
       endpoints = "shared")
}

## (9a) Power against the cutpoint.
panel_real <- expand.grid(dropout = c(0, 0.1, 0.2), window1 = real_windows,
                          label = vapply(real_settings, function(z) z$label,
                                         character(1)),
                          stringsAsFactors = FALSE)
real_power_curve <- do.call(rbind, lapply(seq_len(nrow(panel_real)),
                                          function(i) {
  z <- real_settings[[match(panel_real$label[i],
                            vapply(real_settings, function(w) w$label,
                                   character(1)))]]
  arguments <- real_arguments(z, panel_real$window1[i], panel_real$dropout[i],
                              real_n3)
  fit <- do.call(adaptive2in1count, c(arguments, list(cutpoint = cutpoints)))
  data.frame(label = panel_real$label[i], dropout = panel_real$dropout[i],
             window1 = panel_real$window1[i], fit$results,
             stringsAsFactors = FALSE)
}))
saveRDS(real_power_curve, file.path(data_dir, "real_power_curve.rds"))

## (9b) Operating characteristics at calibrated cutpoints, with the ceiling the
## cutpoint places on the power and the size of the large trial the target power
## requires.
real_table <- do.call(rbind, lapply(real_settings, function(z) {
  arguments <- real_arguments(z, 0.5, real_dropout, real_n3)
  do.call(rbind, lapply(real_targets, function(target) {
    cutpoint <- do.call(adaptive2in1count_cutpoint,
                        c(arguments[setdiff(names(arguments), "alpha")],
                          list(target = target)))$cutpoint
    alternative <- do.call(adaptive2in1count,
                           c(arguments, list(cutpoint = cutpoint)))$results
    null <- do.call(adaptive2in1count,
                    c(real_arguments(z, 0.5, real_dropout, real_n3, ratio = 1),
                      list(cutpoint = cutpoint)))$results
    ceiling <- do.call(adaptive2in1count,
                       c(real_arguments(z, 0.5, real_dropout, 200000),
                         list(cutpoint = cutpoint)))$results
    solved <- try(do.call(adaptive2in1count_size,
                          c(arguments[setdiff(names(arguments), "n3")],
                            list(cutpoint = cutpoint, target = real_power))),
                  silent = TRUE)
    data.frame(label = z$label, target = target, cutpoint = cutpoint,
               type1 = null$rejection, size_null = null$expected_size,
               power = alternative$rejection,
               power_small = alternative$rejection_small,
               ceiling = ceiling$rejection,
               required = if (inherits(solved, "try-error")) {
                 NA_real_
               } else {
                 sum(solved$n3)
               },
               size_alternative = alternative$expected_size,
               stringsAsFactors = FALSE)
  }))
}))
saveRDS(real_table, file.path(data_dir, "real_table.rds"))

message("data generation complete")
