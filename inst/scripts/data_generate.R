# Generates the operating-characteristic data used by the manuscript figures and
# tables. Run from the package root. Results are written as .rds files to
# inst/scripts/data.

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

# The exponential follow-up period is not defined at a dropout proportion of
# zero, where the period is complete by definition.
shape_for <- function(dropout) {
  if (dropout <= 0) "complete" else "exponential"
}

cutpoints <- seq(0.4, 2, by = 0.05)

## (1) Maximum type I error rate against the cutpoint.
message("(1/8) maximum type I error rate against the cutpoint")
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
message("(2/8) power against the cutpoint")
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

## (3) Required size of the large trial and expected sample size against the
## dropout proportion.
message("(3/8) required size against the dropout proportion")
dropout_grid <- seq(0, 0.5, by = 0.05)
panel3 <- expand.grid(dropout = dropout_grid, ratio1 = c(0.55, 0.65, 0.75),
                      ratio2 = c(0.7, 0.75, 0.8))
target_power <- 0.85
size_curve <- do.call(rbind, lapply(seq_len(nrow(panel3)), function(i) {
  shape <- shape_for(panel3$dropout[i])
  cutpoint <- adaptive2in1count_cutpoint(
    rate1 = rate1, rate2 = rate2, dispersion1 = dispersion1,
    dispersion2 = dispersion2, shared = shared, n1 = n1, n2 = n2, n3 = n3,
    target = 0.9, window1 = window1, window2 = window2,
    dropout = panel3$dropout[i], dropout_shape = shape,
    ratio1 = panel3$ratio1[i])$cutpoint
  solved <- try(adaptive2in1count_size(
    rate1 = rate1, rate2 = rate2, dispersion1 = dispersion1,
    dispersion2 = dispersion2, shared = shared, n1 = n1, n2 = n2,
    cutpoint = cutpoint, target = target_power, ratio1 = panel3$ratio1[i],
    ratio2 = panel3$ratio2[i], window1 = window1, window2 = window2,
    dropout = panel3$dropout[i], dropout_shape = shape, alpha = alpha),
    silent = TRUE)
  if (inherits(solved, "try-error")) {
    return(data.frame(dropout = panel3$dropout[i], ratio1 = panel3$ratio1[i],
                      ratio2 = panel3$ratio2[i], cutpoint = cutpoint,
                      required = NA_real_, expected_size = NA_real_,
                      power = NA_real_))
  }
  data.frame(dropout = panel3$dropout[i], ratio1 = panel3$ratio1[i],
             ratio2 = panel3$ratio2[i], cutpoint = cutpoint,
             required = sum(solved$n3), expected_size = solved$expected_size,
             power = solved$power)
}))
saveRDS(size_curve, file.path(data_dir, "size_curve.rds"))

## (4) Within-patient correlation of the two endpoints, by its two sources.
message("(4/8) within-patient correlation of the two endpoints")
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

## (5) Correlation of the interim statistic with the two final statistics.
message("(5/8) correlation between the interim and the final statistics")
panel5 <- expand.grid(n3 = seq(120, 900, by = 20), dropout = c(0, 0.2, 0.4),
                      endpoint = c("phase2", "phase3"),
                      stringsAsFactors = FALSE)
correlation_statistic <- do.call(rbind, lapply(seq_len(nrow(panel5)), function(i) {
  rho <- adaptive2in1count_correlation(
    rate1 = rate1, rate2 = rate2, dispersion1 = dispersion1,
    dispersion2 = dispersion2, shared = shared, n1 = n1, n2 = n2,
    n3 = panel5$n3[i], window1 = window1, window2 = window2,
    dropout = panel5$dropout[i],
    dropout_shape = shape_for(panel5$dropout[i]),
    small_endpoint = panel5$endpoint[i])
  data.frame(n3 = panel5$n3[i], dropout = panel5$dropout[i],
             endpoint = panel5$endpoint[i], endpoint_correlation = rho$endpoint,
             interim_small = rho$interim_small,
             interim_large = rho$interim_large, stringsAsFactors = FALSE)
}))
saveRDS(correlation_statistic, file.path(data_dir, "correlation_statistic.rds"))

## (6) Operating characteristics at cutpoints calibrated to a target
## probability of expansion.
message("(6/8) operating characteristics at calibrated cutpoints")
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

## (7) The ceiling that the cutpoint places on the power.
message("(7/8) the ceiling that the cutpoint places on the power")
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

## (8) The normal approximation against simulation.
message("(8/8) the normal approximation against simulation")
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

message("data generation complete")
