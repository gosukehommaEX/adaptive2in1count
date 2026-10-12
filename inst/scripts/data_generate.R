# Generates the operating characteristics shown by the figures and tables of the
# main text. Run from the package root. Results are written as .rds files to
# inst/scripts/data. The exhibits of the Supporting Information are produced by
# the pair of scripts in inst/supplementary.

library(adaptive2in1count)

data_dir <- file.path("inst", "scripts", "data")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

# Reference design in chronic obstructive pulmonary disease. The Phase 2
# endpoint is the number of moderate or severe exacerbations over half a year
# and the Phase 3 endpoint is the annualized exacerbation rate over a year. The
# Phase 3 rate lies within the negative binomial estimates that Keene et al.
# (2007) report for the TRISTAN study, whose dispersion estimate of 0.46 for the
# variance of the frailty corresponds to a dispersion of about 2.2 here, so the
# dispersions below describe a more heterogeneous population. That study lost
# 39 per cent of the placebo group before the end of the trial, so a dropout
# proportion of 0.3 is taken as the reference.
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
message("(1/10) maximum type I error rate against the cutpoint")
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
message("(2/10) power against the cutpoint")
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

## (3) Within-patient correlation of the two counts, by its two sources, and
## the correlation of the two estimated log rate ratios on the same patients.
message("(3/10) within-patient correlation of the two endpoints")
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
               dropout_shape = shape_for(panel4$dropout[i])),
             # With the same patients behind the interim and the large trial
             # statistics, their correlation is rho_0 of the article. The
             # number of patients cancels, so any common value serves.
             statistic = design_components(
               rate1 = rate1, rate2 = rate2, dispersion1 = dispersion1,
               dispersion2 = panel4$dispersion2[i], shared = panel4$shared[i],
               n1 = 100, n2 = 100, n3 = 100, window1 = panel4$window1[i],
               window2 = window2, dropout = panel4$dropout[i],
               dropout_shape = shape_for(panel4$dropout[i])
             )$correlation[["interim_large"]])
}))
saveRDS(correlation_endpoint, file.path(data_dir, "correlation_endpoint.rds"))

## (4) Correlation of the interim statistic with the two final statistics when
## the two endpoints count the same events over the two windows.
message("(4/10) correlation between the statistics for a shared endpoint")
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
message("(5/10) operating characteristics at calibrated cutpoints")
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
message("(6/10) the ceiling that the cutpoint places on the power")
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
message("(7/10) the normal approximation against simulation")
validation_cutpoints <- c(0.8, 1.2, 1.6)
# A type I error rate near 0.025 has a binomial standard error of about 0.0005
# at one hundred thousand replications, so the scenarios under the global null
# use that number and the scenarios under the alternative use fewer.
null_nsim <- 100000
validation_scenarios <- list(
  list(label = "Reference design, global null", ratio1 = 1, ratio2 = 1,
       dropout = 0, sizes = c(n1, n2, n3), nsim = null_nsim),
  list(label = "Reference design, global null, dropout 0.3", ratio1 = 1,
       ratio2 = 1, dropout = reference_dropout, sizes = c(n1, n2, n3),
       nsim = null_nsim),
  list(label = "Reference design, both endpoints, dropout 0.3",
       ratio1 = alt_ratio1, ratio2 = alt_ratio2, dropout = reference_dropout,
       sizes = c(n1, n2, n3), nsim = 10000),
  list(label = "Reference design, both endpoints, dropout 0.5",
       ratio1 = alt_ratio1, ratio2 = alt_ratio2, dropout = 0.5,
       sizes = c(n1, n2, n3), nsim = 10000),
  list(label = "Small design, global null, dropout 0.3", ratio1 = 1,
       ratio2 = 1, dropout = reference_dropout, sizes = c(30, 45, 200),
       nsim = null_nsim),
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
message("(8/10) the normal approximation for a shared endpoint")
shared_scenarios <- list(
  list(label = "Shared endpoint, global null, dropout 0.3", ratio = 1,
       dropout = reference_dropout, window1 = 0.5, nsim = null_nsim),
  list(label = "Shared endpoint, both windows, dropout 0.3", ratio = 0.75,
       dropout = reference_dropout, window1 = 0.5, nsim = 10000),
  list(label = "Shared endpoint, one window, dropout 0.3", ratio = 0.75,
       dropout = reference_dropout, window1 = window2, nsim = 10000))
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
message("(9/10) the tezepelumab design")
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

## (10) Every number the main text quotes, written as plain text so that the
## manuscript can be checked against the evaluation that produced it rather
## than against a transcription of it.
message("(10/10) numbers quoted by the main text")
res_dir <- file.path("inst", "scripts", "results")
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)
report <- file.path(res_dir, "manuscript_numbers.txt")
con <- file(report, "w")
say <- function(...) writeLines(paste0(...), con)
fmt <- function(x, digits = 4) formatC(x, format = "f", digits = digits)

say("Numbers quoted by the main text")
say("Written by inst/scripts/data_generate.R on ", format(Sys.Date()))
say("")

say("Section 3.2  Type I error rate  (Figure 1)")
say("  combinations evaluated                     ", nrow(typeI_curve))
say("  largest maximum type I error rate          ", fmt(max(typeI_curve$type1)))
for (ep in c("phase2", "phase3")) {
  block <- typeI_curve[typeI_curve$endpoint == ep, ]
  say("  smallest, small trial on ", ep, "          ", fmt(min(block$type1)))
}
key <- interaction(typeI_curve$dispersion2, typeI_curve$endpoint,
                   round(typeI_curve$cutpoint, 3), drop = TRUE)
relative <- vapply(split(typeI_curve$type1, key), function(v) {
  if (length(v) < 2) NA_real_ else (max(v) - min(v)) / min(v)
}, numeric(1))
say("  relative change across dropout, median     ",
    fmt(median(relative, na.rm = TRUE)))
say("  relative change across dropout, largest    ",
    fmt(max(relative, na.rm = TRUE)))
say("")

say("Section 3.3  Power and the cutpoint  (Figure 2, Tables 1 and 2)")
for (r2 in unique(power_curve$ratio2)) {
  for (r1 in unique(power_curve$ratio1)) {
    for (dd in unique(power_curve$dropout)) {
      block <- power_curve[power_curve$ratio1 == r1 & power_curve$ratio2 == r2 &
                             power_curve$dropout == dd, ]
      block <- block[order(block$cutpoint), ]
      trough <- which.min(block$power)
      if (block$power[1] - block$power[trough] > 1e-4) {
        say("  dip: ratio1 ", r1, " ratio2 ", r2, " dropout ", dd,
            "  from ", fmt(block$power[1], 3), " at C = ",
            fmt(block$cutpoint[1], 2), " to ", fmt(block$power[trough], 3),
            " at C = ", fmt(block$cutpoint[trough], 2),
            ", largest power ", fmt(max(block$power), 3))
      }
    }
  }
}
say("  type I error rate over the calibrated cutpoints, smallest ",
    fmt(min(operating_table$type1)), " largest ", fmt(max(operating_table$type1)))
say("  contribution of the small trial, largest    ",
    fmt(max(operating_table$power_small), 3), " smallest ",
    fmt(min(operating_table$power_small), 3))
for (tg in c(0.8, 0.9)) {
  row <- operating_table[which.min(abs(operating_table$target - tg)), ]
  say("  target ", fmt(tg, 2), ": cutpoint ", fmt(row$cutpoint, 3), ", power ",
      fmt(row$power, 3), ", E(N) under H1 ", fmt(row$size_alternative, 0),
      ", E(N) under H0 ", fmt(row$size_null, 0))
}
# Partial null hypothesis: the Phase 2 endpoint carries its design effect and
# the Phase 3 endpoint none. With the small trial on the Phase 2 endpoint only a
# rejection by the large trial is false, and with the small trial on the Phase
# 3 endpoint every rejection is.
at_target <- which.min(abs(operating_table$target - 0.9))
partial_cutpoint <- operating_table$cutpoint[at_target]
partial <- vapply(c("phase2", "phase3"), function(ep) {
  out <- do.call(adaptive2in1count,
                 c(common[setdiff(names(common), "cutpoint")],
                   list(cutpoint = partial_cutpoint, ratio1 = alt_ratio1,
                        ratio2 = 1, small_endpoint = ep)))$results
  if (ep == "phase2") out$rejection_large else out$rejection
}, numeric(1))
say("  partial null, Phase 2 ratio ", alt_ratio1, " and Phase 3 ratio 1, ",
    "at C = ", fmt(partial_cutpoint, 3), ": false rejection with the small ",
    "trial on ",
    "phase2 ", fmt(partial[["phase2"]]), " and on phase3 ",
    fmt(partial[["phase3"]]))
say("  largest ceiling among targets of 0.80 or less ",
    fmt(max(ceiling_table$ceiling[ceiling_table$target <= 0.8]), 4))
attained <- ceiling_table[!is.na(ceiling_table$required), ]
for (i in seq_len(nrow(attained))) {
  say("  target ", fmt(attained$target[i], 2),
      ": patients in the large trial for a power of ", fmt(target_power, 2),
      " ", format(attained$required[i], big.mark = ","))
}
say("")

say("Section 3.4  The two sources of correlation  (Figures 3 and 4)")
for (w in unique(correlation_endpoint$window1)) {
  block <- correlation_endpoint[correlation_endpoint$window1 == w &
                                  correlation_endpoint$dispersion2 ==
                                  dispersion2, ]
  none <- block[block$shared == 0, ]
  say("  window1 ", w, ": correlation with no shared frailty rises from ",
      fmt(min(none$correlation), 3), " to ", fmt(max(none$correlation), 3))
  top <- block[block$shared == max(block$shared), ]
  say("            with the largest shared frailty, from ",
      fmt(top$correlation[order(top$dropout)][1], 3), " to ",
      fmt(top$correlation[order(top$dropout)][nrow(top)], 3))
  middle <- block[abs(block$shared - 0.23) < 1e-8 & block$dropout == 0, ]
  say("            with a shared frailty of 0.23 and no dropout ",
      fmt(middle$correlation[1], 3))
  say("            statistics with no shared frailty, largest absolute value ",
      fmt(max(abs(none$statistic)), 3))
  for (sh in sort(unique(block$shared[block$shared > 0]))) {
    part <- block[abs(block$shared - sh) < 1e-8, ]
    part <- part[order(part$dropout), ]
    say("            statistics with a shared frailty of ", sh, ", from ",
        fmt(part$statistic[1], 3), " to ", fmt(part$statistic[nrow(part)], 3),
        "; counts from ", fmt(part$correlation[1], 3), " to ",
        fmt(part$correlation[nrow(part)], 3))
  }
}
for (w in unique(correlation_shared$window1)) {
  for (dd in unique(correlation_shared$dropout)) {
    block <- correlation_shared[correlation_shared$window1 == w &
                                  correlation_shared$dropout == dd, ]
    block <- block[order(block$n3), ]
    say("  window1 ", w, " dropout ", dd, ": rho_XY ",
        fmt(block$shared_small[1], 3), " (one endpoint) and ",
        fmt(block$distinct_small[1], 3), " (two endpoints);",
        " rho_XZ from ", fmt(block$shared_large[1], 3), " to ",
        fmt(block$shared_large[nrow(block)], 3), " (one endpoint) and from ",
        fmt(block$distinct_large[1], 3), " to ",
        fmt(block$distinct_large[nrow(block)], 3), " (two endpoints)")
  }
}
say("")

say("Section 3.5  Accuracy of the normal approximation  (Table 3)")
say("  cells                                      ", nrow(validation_table))
say("  largest gap in standard errors, dispersion known     ",
    fmt(max(abs(validation_table$known - validation_table$approximate) /
              validation_table$se_known), 2))
say("  largest gap in standard errors, dispersion estimated ",
    fmt(max(abs(validation_table$estimated - validation_table$approximate) /
              validation_table$se_estimated), 2))
say("  cells in which estimating the dispersion raises the rejection ",
    "probability ", sum(validation_table$estimated > validation_table$known),
    " of ", nrow(validation_table))
small_block <- validation_table[grepl("^Small", validation_table$label), ]
say("  small design, smallest and largest increase ",
    fmt(min(small_block$estimated - small_block$known)), " ",
    fmt(max(small_block$estimated - small_block$known)))
say("")

say("Section 4  The tezepelumab design  (Figure 5, Table 4)")
for (lab in unique(real_table$label)) {
  block <- real_table[real_table$label == lab, ]
  say("  ", lab)
  say("    power with the smallest target          ",
      fmt(block$power[which.min(block$target)], 3))
  # Power of the small trial on its own, which the text compares with the
  # power of the design.
  setting <- real_settings[[match(lab, vapply(real_settings,
                                              function(w) w$label,
                                              character(1)))]]
  design <- real_arguments(setting, 0.5, real_dropout, real_n3)
  parts <- do.call(design_components, design[names(design) != "alpha"])
  small_alone <- pnorm(-qnorm(1 - real_alpha) -
                         parts$mean[["small"]] /
                         sqrt(parts$variance[["small"]]))
  say("    power of the small trial alone          ", fmt(small_alone, 3))
  say("    type I error rate, smallest and largest ",
      fmt(min(block$type1)), " ", fmt(max(block$type1)))
  say("    largest target with no attainable size  ",
      fmt(max(block$target[is.na(block$required)]), 2))
  say("    expected sample size, smallest and largest target ",
      fmt(block$size_alternative[which.min(block$target)], 0), " ",
      fmt(block$size_alternative[which.max(block$target)], 0))
  reached <- block[!is.na(block$required), ]
  for (i in seq_len(nrow(reached))) {
    say("    target ", fmt(reached$target[i], 2), ": power ",
        fmt(reached$power[i], 3), ", patients in the large trial for a power ",
        "of ", fmt(real_power, 2), " ", format(reached$required[i],
                                               big.mark = ","),
        ", E(N) ", fmt(reached$size_alternative[i], 0))
  }
}
say("")
# Information about the log rate ratio per pair of patients, one per group.
pair_information <- function(z, dropout) {
  shape <- if (dropout > 0) "uniform" else "complete"
  1 / sum(1 / vapply(c(z$rate * z$ratio, z$rate), function(r) {
    information_factor(r, z$dispersion, 1, 1, dropout, shape)
  }, numeric(1)))
}
for (z in real_settings) {
  for (dd in c(0.1, 0.2)) {
    say("  ", z$label, ": information lost to a uniform dropout proportion ",
        "of ", dd, " ", fmt(1 - pair_information(z, dd) /
                              pair_information(z, 0), 3))
  }
}
say("")
say("  cutpoint at which the power crosses ", fmt(real_power, 2))
for (lab in unique(real_power_curve$label)) {
  setting <- real_settings[[match(lab, vapply(real_settings,
                                              function(w) w$label,
                                              character(1)))]]
  for (dd in unique(real_power_curve$dropout)) {
    crossing <- vapply(sort(unique(real_power_curve$window1)), function(w) {
      block <- real_power_curve[real_power_curve$label == lab &
                                  real_power_curve$dropout == dd &
                                  real_power_curve$window1 == w, ]
      block <- block[order(block$cutpoint), ]
      above <- which(block$rejection >= real_power)[1]
      if (is.na(above) || above == 1) {
        return(NA_real_)
      }
      # The grid brackets the crossing, and the closed form locates it.
      gap <- function(cutpoint) {
        do.call(adaptive2in1count,
                c(real_arguments(setting, w, dd, real_n3),
                  list(cutpoint = cutpoint)))$results$rejection - real_power
      }
      stats::uniroot(gap, block$cutpoint[c(above - 1L, above)],
                     tol = 1e-10)$root
    }, numeric(1))
    say("    ", lab, ", dropout ", dd, ": ",
        paste(fmt(crossing, 3), collapse = " "), " for window1 ",
        paste(sort(unique(real_power_curve$window1)), collapse = " "))
  }
}
close(con)
message("numbers quoted by the main text written to ", report)

message("data generation complete")
