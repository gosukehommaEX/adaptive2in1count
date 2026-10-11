# Reproduces numbers published by others with the functions of the package.
# Run from the package root. The summary is written by the script to
# inst/reproduce/results/reproduce_summary.txt, one line per published number,
# with the verdict PASS, FAIL, EXPLAINED or INFO.
#
# (1) Tang (2015), Tables 1 and 2: required size per group for the negative
#     binomial Wald test with uniform and exponential losses to follow-up. These
#     check information_factor() and both random follow-up distributions.
# (2) Corren et al. (2017): 138 patients per group in the Phase 2 tezepelumab
#     trial, which identifies the dispersion of 0.7 as the variance of the
#     frailty.
# (3) Menzies-Gow et al. (2021), protocol section 9.2: power above 99 per cent
#     and a minimum detectable rate reduction of 27 per cent for the Phase 3
#     tezepelumab trial, which identify the dispersion of 2.4 as the variance of
#     the frailty.

library(adaptive2in1count)

out_dir <- file.path("inst", "reproduce", "results")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

rows <- list()
record <- function(source, item, published, reproduced, verdict, reason = "") {
  rows[[length(rows) + 1L]] <<- data.frame(
    source = source, item = item, published = published,
    reproduced = reproduced, verdict = verdict, reason = reason,
    stringsAsFactors = FALSE)
}

## (1) Tang (2015). The dispersion kappa of Tang is the variance of the frailty,
## so it enters information_factor() as its reciprocal. Two-sided level 0.05
## and power 0.9 throughout.
z_tang <- qnorm(0.975) + qnorm(0.9)
tang_size <- function(rate0, ratio, kappa, period, dropout, shape) {
  first <- information_factor(rate0, 1 / kappa, period, period, dropout, shape)
  second <- information_factor(rate0 * ratio, 1 / kappa, period, period,
                               dropout, shape)
  ceiling((1 / first + 1 / second) * z_tang^2 / log(ratio)^2)
}

# Table 2 of Tang (2015): two years of follow-up and a proportion w_c of
# patients who leave early.
tang2 <- data.frame(
  withdraw = rep(c(0.2, 0.3, 0.4), each = 8),
  kappa = rep(c(1, 1.5, 1, 1.5, 1, 1.5, 1, 1.5), times = 3),
  rate0 = rep(c(0.4, 0.8, 0.4, 0.8, 0.4, 0.8, 0.4, 0.8), times = 3),
  ratio = rep(c(0.5, 0.5, 0.5, 0.5, 0.7, 0.7, 0.7, 0.7), times = 3),
  shape = rep(c("uniform", "uniform", "exponential", "exponential",
                "uniform", "uniform", "exponential", "exponential"),
              times = 3),
  published = c(138, 115, 139, 115, 454, 400, 456, 401,
                145, 119, 146, 120, 476, 414, 481, 418,
                152, 124, 155, 126, 499, 429, 510, 437),
  stringsAsFactors = FALSE)
for (i in seq_len(nrow(tang2))) {
  size <- tang_size(tang2$rate0[i], tang2$ratio[i], tang2$kappa[i], 2,
                    tang2$withdraw[i], tang2$shape[i])
  record("Tang (2015), Table 2",
         sprintf("w_c %.1f, kappa %.1f, rate %.1f, ratio %.1f, %s",
                 tang2$withdraw[i], tang2$kappa[i], tang2$rate0[i],
                 tang2$ratio[i], tang2$shape[i]),
         tang2$published[i], size,
         if (size == tang2$published[i]) "PASS" else "FAIL")
}

# Table 1 of Tang (2015): exponential loss to follow-up with mean 6.952 years,
# placebo rate 0.32 and a 40 per cent reduction, over one or two years.
tang1 <- data.frame(kappa = c(1, 1, 1.25, 1.25), period = c(2, 1, 2, 1),
                    published = c(281, 445, 303, 466))
for (i in seq_len(nrow(tang1))) {
  dropout <- 1 - exp(-tang1$period[i] / 6.952)
  size <- tang_size(0.32, 0.6, tang1$kappa[i], tang1$period[i], dropout,
                    "exponential")
  record("Tang (2015), Table 1",
         sprintf("kappa %.2f, %d year(s)", tang1$kappa[i], tang1$period[i]),
         tang1$published[i], size,
         if (size == tang1$published[i]) "PASS" else "FAIL")
}

## (2) Corren et al. (2017): 80 per cent power, two-sided level 0.1, placebo
## rate 0.7 per year over a year, a 40 per cent reduction, dispersion 0.7 and a
## 10 per cent loss of information to dropout.
z_phase2 <- qnorm(0.95) + qnorm(0.8)
phase2_size <- function(dispersion) {
  first <- information_factor(0.7, dispersion, 1, 1)
  second <- information_factor(0.7 * 0.6, dispersion, 1, 1)
  ceiling((1 / first + 1 / second) * z_phase2^2 / log(0.6)^2 / 0.9)
}
size <- phase2_size(1 / 0.7)
record("Corren et al. (2017)", "patients per group, 0.7 read as the variance",
       138, size, if (size == 138) "PASS" else "FAIL")
record("Corren et al. (2017)", "patients per group, 0.7 read as the shape",
       138, phase2_size(0.7), "INFO",
       "the other reading of the dispersion does not give 138")

## (3) Menzies-Gow et al. (2021), protocol: 530 per group, placebo rate 0.9 per
## year, dispersion 2.4, 10 per cent dropout uniform over the year, two-sided
## level 0.01.
z_phase3 <- qnorm(0.995)
phase3_error <- function(ratio, dispersion) {
  first <- information_factor(0.9, dispersion, 1, 1, 0.1, "uniform")
  second <- information_factor(0.9 * ratio, dispersion, 1, 1, 0.1, "uniform")
  sqrt((1 / first + 1 / second) / 530)
}
detectable <- function(dispersion) {
  ratio <- 0.7
  for (step in seq_len(100)) {
    ratio <- exp(-z_phase3 * phase3_error(ratio, dispersion))
  }
  1 - ratio
}
power <- pnorm(log(2) / phase3_error(0.5, 1 / 2.4) - z_phase3)
record("Menzies-Gow et al. (2021), protocol",
       "power against a 50 per cent reduction, 2.4 read as the variance",
       "> 0.99", round(power, 4), if (power > 0.99) "PASS" else "FAIL")
reduction <- detectable(1 / 2.4)
record("Menzies-Gow et al. (2021), protocol",
       "minimum detectable reduction, 2.4 read as the variance",
       0.27, round(reduction, 4),
       if (abs(round(reduction, 2) - 0.27) < 1e-9) "PASS" else "FAIL")
record("Menzies-Gow et al. (2021), protocol",
       "minimum detectable reduction, 2.4 read as the shape",
       0.27, round(detectable(2.4), 4), "INFO",
       "the other reading of the dispersion does not give 27 per cent")

summary_table <- do.call(rbind, rows)
saveRDS(summary_table, file.path(out_dir, "reproduce_summary.rds"))
report <- file.path(out_dir, "reproduce_summary.txt")
con <- file(report, "w")
writeLines(c("Reproduction of published numbers",
             paste0("Written by inst/reproduce/reproduce_published.R on ",
                    format(Sys.Date())),
             "",
             paste(names(summary_table), collapse = " | ")), con)
writeLines(do.call(paste, c(summary_table, sep = " | ")), con)
counts <- table(factor(summary_table$verdict,
                       levels = c("PASS", "EXPLAINED", "INFO", "FAIL")))
writeLines(c("", paste(names(counts), counts, sep = " ", collapse = ", ")), con)
close(con)
message("reproduction summary written to ", report)
