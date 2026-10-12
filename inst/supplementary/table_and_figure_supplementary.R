# Builds the exhibits of the Supporting Information from the data written by
# data_generate.R in this directory. Run from the package root. Figures are
# written as .pdf and .eps and tables as .tex to inst/supplementary/results.
# The exhibits of the main text are built by the pair of scripts in
# inst/scripts.

library(ggplot2)

data_dir <- file.path("inst", "supplementary", "data")
res_dir <- file.path("inst", "supplementary", "results")
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)

pal2 <- c("#1B9E77", "#D95F02")
pal3 <- c("#1B9E77", "#D95F02", "#7570B3")
pal4 <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A")
lty2 <- c("solid", "dashed")
lty3 <- c("solid", "dashed", "dotdash")
lty4 <- c("solid", "dashed", "dotdash", "dotted")
shape3 <- c(16, 17, 15)
shape4 <- c(16, 17, 15, 18)

# The journal prints a figure across a text width of about six inches, so the
# canvas is kept close to that width and the type is set relative to it. A
# figure saved twelve inches wide at twenty point therefore carries about ten
# point type once it is scaled to the text width, which matches the body text.
theme_ms <- theme_bw(base_size = 20) +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.key.size = grid::unit(1.2, "cm"),
        legend.text = element_text(size = rel(0.95)),
        legend.margin = margin(0, 0, 0, 0),
        axis.title = element_text(size = rel(1.0)),
        axis.text = element_text(size = rel(0.9)),
        strip.background = element_rect(fill = "grey92"),
        strip.text = element_text(size = rel(0.85)),
        panel.grid.minor = element_blank(),
        plot.margin = margin(6, 12, 6, 6))

# Graphs are supplied as .eps at 600 dpi or better, which is the minimum the
# journal accepts for line art, and as .pdf for the working manuscript.
save_fig <- function(plot, name, width = 12, height = 9) {
  ggsave(file.path(res_dir, paste0(name, ".pdf")), plot, width = width,
         height = height, units = "in", device = cairo_pdf)
  ggsave(file.path(res_dir, paste0(name, ".eps")), plot, width = width,
         height = height, units = "in", device = cairo_ps,
         fallback_resolution = 600)
}

reference_line <- function(y) {
  geom_hline(yintercept = y, colour = "grey45", linewidth = 0.7,
             linetype = "twodash")
}

scale_cutpoint <- scale_x_continuous(limits = c(0.4, 2),
                                     breaks = seq(0.5, 2, by = 0.5))

# Writes a table environment that can be uploaded to the manuscript source
# without further editing. Only standard LaTeX commands are used, so no
# additional package is required.
write_tex_table <- function(mat, headers, aligns, file, caption, label,
                            clines = character(0), size = "\\small",
                            colsep = "4pt") {
  con <- file(file, "w")
  on.exit(close(con))
  caption <- sub("[.[:space:]]*$", "", caption)
  writeLines("\\begin{table}[htbp]", con)
  writeLines("\\centering", con)
  writeLines(paste0("\\caption{", caption, ".}"), con)
  writeLines(paste0("\\label{", label, "}"), con)
  writeLines("\\vspace{0.6em}", con)
  writeLines(size, con)
  writeLines(paste0("\\setlength{\\tabcolsep}{", colsep, "}"), con)
  writeLines(paste0("\\begin{tabular}{", aligns, "}"), con)
  writeLines("\\hline", con)
  for (k in seq_along(headers)) {
    writeLines(paste0(paste(headers[[k]], collapse = " & "), " \\\\"), con)
    if (k <= length(clines) && nzchar(clines[k])) {
      writeLines(clines[k], con)
    }
  }
  writeLines("\\hline", con)
  for (i in seq_len(nrow(mat))) {
    writeLines(paste0(paste(mat[i, ], collapse = " & "), " \\\\"), con)
  }
  writeLines("\\hline", con)
  writeLines("\\end{tabular}", con)
  writeLines("\\end{table}", con)
}

fmt <- function(x, digits) formatC(x, format = "f", digits = digits)

## Figure S1. Required size of the large trial and expected sample size against
## the dropout proportion.
d3 <- readRDS(file.path(data_dir, "size_curve.rds"))
keys3 <- d3[c("dropout", "ratio1", "ratio2")]
d3_long <- rbind(
  data.frame(keys3, metric = "Large~trial~size",
             value = d3$required),
  data.frame(keys3, metric = "Expected~sample~size",
             value = d3$expected_size))
d3_long$metric <- factor(d3_long$metric,
                         levels = c("Large~trial~size",
                                    "Expected~sample~size"))
d3_long$effect <- factor(d3_long$ratio1, levels = c(0.55, 0.65, 0.75))
d3_long$col <- factor(d3_long$ratio2, levels = c(0.7, 0.75, 0.8),
                      labels = c("Phase~3~ratio == 0.70",
                                 "Phase~3~ratio == 0.75",
                                 "Phase~3~ratio == 0.80"))
lab_effect <- parse(text = c("Phase~2~rate~ratio == 0.55",
                             "Phase~2~rate~ratio == 0.65",
                             "Phase~2~rate~ratio == 0.75"))
fig3 <- ggplot(d3_long, aes(dropout, value, colour = effect, linetype = effect,
                            shape = effect)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3.2) +
  facet_grid(metric ~ col, scales = "free_y", labeller = label_parsed) +
  scale_colour_manual(values = pal3, labels = lab_effect) +
  scale_linetype_manual(values = lty3, labels = lab_effect) +
  scale_shape_manual(values = shape3, labels = lab_effect) +
  labs(x = expression("Dropout proportion " * italic(d)), y = NULL) +
  theme_ms
save_fig(fig3, "figureS1_required_size")

## Figure S2. Correlation of the interim statistic with the two final
## statistics, for two endpoints that count different events.
d5 <- readRDS(file.path(data_dir, "correlation_statistic.rds"))
keys5 <- d5[c("n3", "dropout", "endpoint")]
d5_long <- rbind(
  data.frame(keys5, statistic = "small", value = d5$interim_small),
  data.frame(keys5, statistic = "large", value = d5$interim_large))
d5_long$statistic <- factor(d5_long$statistic, levels = c("small", "large"))
d5_long$row <- factor(d5_long$endpoint, levels = c("phase2", "phase3"),
                      labels = c("Small~trial~on~Phase~2",
                                 "Small~trial~on~Phase~3"))
d5_long$col <- factor(d5_long$dropout, levels = c(0, 0.2, 0.4),
                      labels = c("No~dropout", "italic(d) == 0.2",
                                 "italic(d) == 0.4"))
lab_stat <- expression(rho[XY], rho[XZ])
fig5 <- ggplot(d5_long, aes(n3, value, colour = statistic,
                            linetype = statistic)) +
  geom_line(linewidth = 1.2) +
  facet_grid(row ~ col, scales = "free_y", labeller = label_parsed) +
  scale_colour_manual(values = pal2, labels = lab_stat) +
  scale_linetype_manual(values = lty2, labels = lab_stat) +
  expand_limits(y = 0) +
  labs(x = expression("Patients per group in the large trial, " * italic(n)[3]),
       y = "Correlation with the interim statistic") +
  theme_ms
save_fig(fig5, "figureS2_statistic_correlation")

## Table S1. Measured discrepancy of the closed-form correlations.
s1 <- readRDS(file.path(data_dir, "correlation_check.rds"))
s1 <- s1[order(s1$structure, s1$window1, s1$dropout), , drop = FALSE]
structure_label <- ifelse(s1$structure == "shared", "One endpoint",
                          "Two endpoints")
matS1 <- cbind(ifelse(duplicated(structure_label), "", structure_label),
               fmt(s1$window1, 2), fmt(s1$dropout, 1),
               fmt(s1$closed_small, 4), fmt(s1$observed_small, 4),
               fmt(s1$closed_large, 4), fmt(s1$observed_large, 4),
               fmt(s1$largest_z, 1))
headS1 <- list(
  c("", "", "", "\\multicolumn{2}{c}{$\\rho_{XY}$}",
    "\\multicolumn{2}{c}{$\\rho_{XZ}$}", ""),
  c("Structure", "$\\tau_1$", "$d$", "Closed form", "Simulated",
    "Closed form", "Simulated", "Largest $|z|$"))
# The design of the check is read from the data, so that the caption cannot
# drift from the values that were simulated.
two <- unique(s1[s1$structure == "distinct",
                 c("rate1", "rate2", "dispersion1", "dispersion2", "shared",
                   "n1", "n2", "n3")])
one <- unique(s1[s1$structure == "shared", c("rate2", "dispersion2")])
stopifnot(nrow(two) == 1, nrow(one) == 1)
planned <- unique(s1$planned)
stopifnot(length(planned) == 1)
usable_text <- if (all(s1$replications == planned)) {
  ""
} else {
  paste0(", of which at least ",
         format(min(s1$replications), big.mark = ",", scientific = FALSE),
         " gave all three statistics")
}
capS1 <- paste(
  "Correlations of the interim statistic with the two final statistics",
  "obtained in closed form and from generated patients. The generator repeats",
  "the model rather than calling the closed forms. Each statistic is the log",
  "rate ratio of the negative binomial model with an offset for the exposure",
  "and the dispersion held at its design value, and both groups have the",
  "rates given below. The small trial is analyzed on the Phase 2 endpoint,",
  "or on the interim window of the single endpoint (configuration A). Each",
  "value uses",
  paste0(format(planned, big.mark = ",", scientific = FALSE),
         " replications", usable_text, " of a trial that randomizes ", two$n1,
         " patients per"),
  paste0("group before the interim analysis, ", two$n2, " per group in the",
         " small trial and ", two$n3, " per group in the large trial, at a",
         " rate of ", two$rate1, " events per year and dispersion ",
         two$dispersion1, " for the Phase 2 endpoint and ", two$rate2,
         " events per year and dispersion ", two$dispersion2, " for the",
         " Phase 3 endpoint, with a shared frailty shape of ", two$shared,
         " for two endpoints and a single process of ", one$rate2,
         " events per year and dispersion ", one$dispersion2, " for one"),
  "endpoint. The follow-up period is complete when $d = 0$ and exponential",
  "truncated at the year otherwise. The last column is the larger of the two",
  "absolute differences between the simulation and the closed form, each",
  "divided by the Monte Carlo standard error $(1 - \\rho^2) / \\sqrt{R}$ of a",
  "sample correlation $\\rho$ at $R$ replications")
write_tex_table(matS1, headS1, "lrrrrrrr",
                file.path(res_dir, "tableS1_correlation_check.tex"),
                capS1, "tab:check",
                clines = c("\\cline{4-5}\\cline{6-7}", ""),
                size = "\\footnotesize", colsep = "3pt")

message("exhibits of the Supporting Information written to ", res_dir)
