# Generates the operating characteristics shown by the figures and tables of the
# Supporting Information. Run from the package root. Results are written as .rds
# files to inst/supplementary/data. The exhibits of the main text are produced
# by the pair of scripts in inst/scripts.

library(adaptive2in1count)

data_dir <- file.path("inst", "supplementary", "data")
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

## (1) Required size of the large trial and expected sample size against the
## dropout proportion.
message("(1/4) required size against the dropout proportion")
dropout_grid <- seq(0, 0.5, by = 0.05)
panel3 <- expand.grid(dropout = dropout_grid, ratio1 = c(0.55, 0.65, 0.75),
                      ratio2 = c(0.7, 0.75, 0.8))
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

## (2) Correlation of the interim statistic with the two final statistics.
message("(2/4) correlation between the interim and the final statistics")
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

## (3) Discrepancy between the closed-form correlations of the three
## statistics and the correlations of statistics computed from generated
## patients. The generator repeats the model rather than calling the closed
## forms, so the measured discrepancy is an independent check.
message("(3/4) measured discrepancy of the closed-form correlations")
set.seed(20260808)
# Fifteen thousand replications leave a Monte Carlo error of about 0.008 on
# each sample correlation, which is well below the discrepancy being measured.
replications <- 15000
check_sizes <- c(40, 60, 150)
draw_follow_up_check <- function(size, dropout, shape) {
  if (shape == "complete" || dropout == 0) {
    return(rep(window2, size))
  }
  if (shape == "uniform") {
    early <- stats::runif(size) < dropout
    return(ifelse(early, stats::runif(size, 0, window2), window2))
  }
  pmin(stats::rexp(size, rate = -log(1 - dropout) / window2), window2)
}
draw_check <- function(size, w1, dropout, shape, structure) {
  follow_up <- draw_follow_up_check(size, dropout, shape)
  exposure1 <- pmin(follow_up, w1)
  if (structure == "shared") {
    frailty <- stats::rgamma(size, shape = dispersion2, rate = dispersion2)
    count1 <- stats::rpois(size, rate2 * exposure1 * frailty)
    count2 <- count1 + stats::rpois(size, rate2 * (follow_up - exposure1) *
                                      frailty)
  } else {
    base <- stats::rgamma(size, shape = shared, rate = 1)
    count1 <- stats::rpois(size, rate1 * exposure1 *
                             (base + stats::rgamma(size,
                                                   shape = dispersion1 - shared,
                                                   rate = 1)) / dispersion1)
    count2 <- stats::rpois(size, rate2 * follow_up *
                             (base + stats::rgamma(size,
                                                   shape = dispersion2 - shared,
                                                   rate = 1)) / dispersion2)
  }
  list(count1 = count1, count2 = count2, exposure1 = exposure1,
       exposure2 = follow_up)
}
check_settings <- expand.grid(structure = c("distinct", "shared"),
                              dropout = c(0, 0.4), window1 = c(0.5, 1),
                              stringsAsFactors = FALSE)
correlation_check <- do.call(rbind, lapply(seq_len(nrow(check_settings)),
                                           function(i) {
  z <- check_settings[i, ]
  shape <- shape_for(z$dropout)
  groups <- lapply(seq_len(2), function(g) {
    patients <- draw_check(replications * check_sizes[3], z$window1, z$dropout,
                           shape, z$structure)
    lapply(patients, function(x) matrix(x, nrow = replications))
  })
  contrast <- function(size, endpoint) {
    ratio <- vapply(seq_len(2), function(g) {
      columns <- seq_len(size)
      rowSums(groups[[g]][[paste0("count", endpoint)]][, columns,
                                                       drop = FALSE]) /
        rowSums(groups[[g]][[paste0("exposure", endpoint)]][, columns,
                                                            drop = FALSE])
    }, numeric(replications))
    log(ratio[, 1]) - log(ratio[, 2])
  }
  statistic <- cbind(contrast(check_sizes[1], 1L), contrast(check_sizes[2], 1L),
                     contrast(check_sizes[3], 2L))
  usable <- apply(is.finite(statistic), 1, all)
  observed <- stats::cor(statistic[usable, ])
  expected <- design_components(
    rate1 = if (z$structure == "shared") rate2 else rate1, rate2 = rate2,
    dispersion1 = if (z$structure == "shared") dispersion2 else dispersion1,
    dispersion2 = dispersion2,
    shared = if (z$structure == "shared") dispersion2 else shared,
    n1 = check_sizes[1], n2 = check_sizes[2], n3 = check_sizes[3],
    window1 = z$window1, window2 = window2, dropout = z$dropout,
    dropout_shape = shape, endpoints = z$structure)$correlation
  data.frame(structure = z$structure, dropout = z$dropout,
             window1 = z$window1, replications = sum(usable),
             rate1 = if (z$structure == "shared") rate2 else rate1,
             rate2 = rate2,
             dispersion1 = if (z$structure == "shared") dispersion2 else
               dispersion1,
             dispersion2 = dispersion2,
             shared = if (z$structure == "shared") dispersion2 else shared,
             n1 = check_sizes[1], n2 = check_sizes[2], n3 = check_sizes[3],
             observed_small = observed[1, 2], observed_large = observed[1, 3],
             closed_small = expected[["interim_small"]],
             closed_large = expected[["interim_large"]],
             discrepancy = max(abs(observed[1, 2] -
                                     expected[["interim_small"]]),
                               abs(observed[1, 3] -
                                     expected[["interim_large"]])),
             stringsAsFactors = FALSE)
}))
saveRDS(correlation_check, file.path(data_dir, "correlation_check.rds"))


## (4) Every number the Supporting Information quotes, written as plain text so
## that the document can be checked against the evaluation that produced it.
message("(4/4) numbers quoted by the Supporting Information")
res_dir <- file.path("inst", "supplementary", "results")
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)
report <- file.path(res_dir, "supplementary_numbers.txt")
con <- file(report, "w")
say <- function(...) writeLines(paste0(...), con)
fmt <- function(x, digits = 4) formatC(x, format = "f", digits = digits)

say("Numbers quoted by the Supporting Information")
say("Written by inst/supplementary/data_generate.R on ", format(Sys.Date()))
say("")

say("Figure S1  Required size against the dropout proportion")
for (r1 in unique(size_curve$ratio1)) {
  for (r2 in unique(size_curve$ratio2)) {
    block <- size_curve[size_curve$ratio1 == r1 & size_curve$ratio2 == r2, ]
    block <- block[order(block$dropout), ]
    # The dropout grid is built by seq(), so a value is matched by distance
    # rather than by equality.
    at <- function(value) which.min(abs(block$dropout - value))
    base <- block$required[at(0)]
    say("  ratio1 ", r1, " ratio2 ", r2, ": large trial from ",
        format(base, big.mark = ","), " with no dropout to ",
        format(block$required[at(0.3)], big.mark = ","),
        " at a dropout proportion of 0.3 and ",
        format(block$required[nrow(block)], big.mark = ","), " at ",
        block$dropout[nrow(block)], ", a rise of ",
        fmt(100 * (block$required[at(0.3)] / base - 1), 0),
        " and ", fmt(100 * (block$required[nrow(block)] / base - 1), 0),
        " per cent")
    say("            expected sample size from ",
        format(round(block$expected_size[at(0)]), big.mark = ","), " to ",
        format(round(block$expected_size[at(0.3)]), big.mark = ","), " and ",
        format(round(block$expected_size[nrow(block)]), big.mark = ","))
  }
}
say("")

say("Figure S2  Correlation of the interim with the two final statistics")
for (ep in unique(correlation_statistic$endpoint)) {
  for (dd in unique(correlation_statistic$dropout)) {
    block <- correlation_statistic[correlation_statistic$endpoint == ep &
                                     correlation_statistic$dropout == dd, ]
    block <- block[order(block$n3), ]
    say("  small trial on ", ep, ", dropout ", dd, ": rho_XY ",
        fmt(block$interim_small[1], 3), ", rho_XZ from ",
        fmt(block$interim_large[1], 3), " to ",
        fmt(block$interim_large[nrow(block)], 3),
        ", within-patient correlation ", fmt(block$endpoint_correlation[1], 3))
  }
}
say("")

say("Table S1  Measured gap between the closed forms and the simulation")
for (i in seq_len(nrow(correlation_check))) {
  say("  ", correlation_check$structure[i], ", window1 ",
      correlation_check$window1[i], ", dropout ",
      correlation_check$dropout[i], ": largest gap ",
      fmt(correlation_check$discrepancy[i]))
}
say("  largest gap under complete follow-up       ",
    fmt(max(correlation_check$discrepancy[correlation_check$dropout == 0])))
for (st in unique(correlation_check$structure)) {
  block <- correlation_check[correlation_check$structure == st &
                               correlation_check$dropout > 0, ]
  say("  gaps with dropout, ", st, "                ",
      fmt(min(block$discrepancy)), " to ", fmt(max(block$discrepancy)))
}
close(con)
message("numbers quoted by the Supporting Information written to ", report)

message("data generation complete")
