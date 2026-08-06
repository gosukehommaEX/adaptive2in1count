# Builds the manuscript exhibits (five figures and three tables) from the data
# written by data_generate.R. Run from the package root. Figures are written as
# .pdf and .eps and tables as .tex to inst/scripts/results.

library(ggplot2)

data_dir <- file.path("inst", "scripts", "data")
res_dir <- file.path("inst", "scripts", "results")
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)

pal2 <- c("#1B9E77", "#D95F02")
pal3 <- c("#1B9E77", "#D95F02", "#7570B3")
pal4 <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A")
lty2 <- c("solid", "dashed")
lty3 <- c("solid", "dashed", "dotdash")
lty4 <- c("solid", "dashed", "dotdash", "dotted")
shape3 <- c(16, 17, 15)
shape4 <- c(16, 17, 15, 18)

theme_ms <- theme_bw(base_size = 20) +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.key.size = grid::unit(2, "cm"),
        legend.text = element_text(size = rel(0.95)),
        strip.background = element_rect(fill = "grey92"),
        panel.grid.minor = element_blank(),
        plot.margin = margin(6, 12, 6, 6))

# Graphs are supplied as .eps at 600 dpi or better, which is the minimum the
# journal accepts for line art, and as .pdf for the working manuscript.
save_fig <- function(plot, name, width = 16, height = 12) {
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
                            clines = character(0)) {
  con <- file(file, "w")
  on.exit(close(con))
  caption <- sub("[.[:space:]]*$", "", caption)
  writeLines("\\begin{table}[htbp]", con)
  writeLines("\\centering", con)
  writeLines(paste0("\\caption{", caption, ".}"), con)
  writeLines(paste0("\\label{", label, "}"), con)
  writeLines("\\vspace{0.6em}", con)
  writeLines("\\small", con)
  writeLines("\\setlength{\\tabcolsep}{4pt}", con)
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

## Figure 1. Maximum type I error rate against the cutpoint.
d1 <- readRDS(file.path(data_dir, "typeI_curve.rds"))
d1$dispersion <- factor(d1$dispersion2, levels = c(0.5, 0.8, 1.5))
d1$row <- factor(d1$endpoint, levels = c("phase2", "phase3"),
                 labels = c("Small~trial~on~the~Phase~2~endpoint",
                            "Small~trial~on~the~Phase~3~endpoint"))
d1$col <- factor(d1$dropout, levels = c(0, 0.2, 0.4),
                 labels = c("No~dropout", "italic(d) == 0.2",
                            "italic(d) == 0.4"))
lab_disp2 <- parse(text = c("kappa[2] == 0.5", "kappa[2] == 0.8",
                            "kappa[2] == 1.5"))
fig1 <- ggplot(d1, aes(cutpoint, type1, colour = dispersion,
                       linetype = dispersion)) +
  reference_line(0.025) +
  geom_line(linewidth = 1.2) +
  facet_grid(row ~ col, labeller = label_parsed) +
  scale_colour_manual(values = pal3, labels = lab_disp2) +
  scale_linetype_manual(values = lty3, labels = lab_disp2) +
  scale_cutpoint +
  labs(x = expression("Cutpoint " * italic(C) *
                        " for the observed rate ratio"),
       y = "Maximum type I error rate") +
  theme_ms
save_fig(fig1, "figure1_type_one_error")

## Figure 2. Power against the cutpoint, by treatment effect and dropout.
d2 <- readRDS(file.path(data_dir, "power_curve.rds"))
d2$loss <- factor(d2$dropout, levels = c(0, 0.2, 0.4))
d2$row <- factor(d2$ratio2, levels = c(0.7, 0.8),
                 labels = c("Phase~3~rate~ratio == 0.70",
                            "Phase~3~rate~ratio == 0.80"))
d2$col <- factor(d2$ratio1, levels = c(0.55, 0.65, 0.75),
                 labels = c("Phase~2~rate~ratio == 0.55",
                            "Phase~2~rate~ratio == 0.65",
                            "Phase~2~rate~ratio == 0.75"))
lab_loss <- parse(text = c("No~dropout", "italic(d) == 0.2",
                           "italic(d) == 0.4"))
fig2 <- ggplot(d2, aes(cutpoint, power, colour = loss, linetype = loss)) +
  reference_line(0.85) +
  geom_line(linewidth = 1.2) +
  facet_grid(row ~ col, labeller = label_parsed) +
  scale_colour_manual(values = pal3, labels = lab_loss) +
  scale_linetype_manual(values = lty3, labels = lab_loss) +
  scale_cutpoint +
  scale_y_continuous(breaks = seq(0.2, 1, by = 0.2)) +
  labs(x = expression("Cutpoint " * italic(C) *
                        " for the observed rate ratio"),
       y = "Power") +
  theme_ms
save_fig(fig2, "figure2_power")

## Figure 3. Required size of the large trial and expected sample size against
## the dropout proportion.
d3 <- readRDS(file.path(data_dir, "size_curve.rds"))
keys3 <- d3[c("dropout", "ratio1", "ratio2")]
d3_long <- rbind(
  data.frame(keys3, metric = "Required~size~of~the~large~trial",
             value = d3$required),
  data.frame(keys3, metric = "Expected~total~sample~size",
             value = d3$expected_size))
d3_long$metric <- factor(d3_long$metric,
                         levels = c("Required~size~of~the~large~trial",
                                    "Expected~total~sample~size"))
d3_long$effect <- factor(d3_long$ratio1, levels = c(0.55, 0.65, 0.75))
d3_long$col <- factor(d3_long$ratio2, levels = c(0.7, 0.75, 0.8),
                      labels = c("Phase~3~rate~ratio == 0.70",
                                 "Phase~3~rate~ratio == 0.75",
                                 "Phase~3~rate~ratio == 0.80"))
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
save_fig(fig3, "figure3_required_size")

## Figure 4. Within-patient correlation of the two endpoints, by its two
## sources.
d4 <- readRDS(file.path(data_dir, "correlation_endpoint.rds"))
d4$frailty <- factor(d4$shared, levels = c(0, 0.23, 0.47, 0.7))
d4$row <- factor(d4$dispersion2, levels = c(0.8, 1.5),
                 labels = c("kappa[2] == 0.8", "kappa[2] == 1.5"))
d4$col <- factor(d4$window1, levels = c(0.25, 0.5, 1),
                 labels = c("tau[1] == 0.25", "tau[1] == 0.50",
                            "tau[1] == 1.00"))
lab_frailty <- parse(text = c("nu[0] == 0", "nu[0] == 0.23",
                              "nu[0] == 0.47", "nu[0] == 0.70"))
fig4 <- ggplot(d4, aes(dropout, correlation, colour = frailty,
                       linetype = frailty, shape = frailty)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3, data = subset(d4, round(dropout * 100) %% 10 == 0)) +
  facet_grid(row ~ col, labeller = label_parsed) +
  scale_colour_manual(values = pal4, labels = lab_frailty) +
  scale_linetype_manual(values = lty4, labels = lab_frailty) +
  scale_shape_manual(values = shape4, labels = lab_frailty) +
  labs(x = expression("Dropout proportion " * italic(d)),
       y = "Within-patient correlation of the two endpoints") +
  theme_ms
save_fig(fig4, "figure4_endpoint_correlation")

## Figure 5. Correlation of the interim statistic with the two final
## statistics.
d5 <- readRDS(file.path(data_dir, "correlation_statistic.rds"))
keys5 <- d5[c("n3", "dropout", "endpoint")]
d5_long <- rbind(
  data.frame(keys5, statistic = "small", value = d5$interim_small),
  data.frame(keys5, statistic = "large", value = d5$interim_large))
d5_long$statistic <- factor(d5_long$statistic, levels = c("small", "large"))
d5_long$row <- factor(d5_long$endpoint, levels = c("phase2", "phase3"),
                      labels = c("Small~trial~on~the~Phase~2~endpoint",
                                 "Small~trial~on~the~Phase~3~endpoint"))
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
save_fig(fig5, "figure5_statistic_correlation")

## Table 1. Operating characteristics at calibrated cutpoints.
t1 <- readRDS(file.path(data_dir, "operating_table.rds"))
mat1 <- cbind(fmt(t1$target, 2), fmt(t1$cutpoint, 3), fmt(t1$go_null, 4),
              fmt(t1$type1, 4), fmt(t1$size_null, 1),
              fmt(t1$go_alternative, 4), fmt(t1$power, 4),
              fmt(t1$power_small, 4), fmt(t1$power_large, 4),
              fmt(t1$size_alternative, 1))
head1 <- list(
  c("\\multicolumn{2}{c}{Cutpoint}", "\\multicolumn{3}{c}{Under $H_0$}",
    "\\multicolumn{5}{c}{Under $H_1$}"),
  c("Target", "$C$", "$\\Pr(\\mathrm{expand})$", "Type I", "$E(N)$",
    "$\\Pr(\\mathrm{expand})$", "Power", "Small", "Large", "$E(N)$"))
cap1 <- paste(
  "Operating characteristics of the 2-in-1 adaptive design for count endpoints",
  "in the reference chronic obstructive pulmonary disease setting. The Phase 2",
  "endpoint has rate $\\theta_1 = 2$ events per year, dispersion",
  "$\\kappa_1 = 0.7$ and observation window $\\tau_1 = 0.5$ years, and the",
  "Phase 3 endpoint has rate $\\theta_2 = 1.3$ events per year, dispersion",
  "$\\kappa_2 = 0.8$ and observation window $\\tau_2 = 1$ year. The shape of the",
  "frailty component shared by the two endpoints is $\\nu_0 = 0.4$. The",
  "follow-up period is exponential truncated at $\\tau_2$ with dropout",
  "proportion $d = 0.3$. The trial randomises 80 patients per group before the",
  "interim analysis, 120 per group in the small trial and 750 per group in the",
  "large trial, and the Wald test is applied at a one-sided significance level",
  "of 0.025. The trial expands when the observed rate ratio of the Phase 2",
  "endpoint is at or below the cutpoint $C$, which is chosen so that the",
  "probability of expansion attains the value in the column headed target when",
  "the Phase 2 rate ratio is 0.65. $H_0$ sets both rate ratios to one and $H_1$",
  "sets the Phase 2 rate ratio to 0.65 and the Phase 3 rate ratio to 0.75. The",
  "columns headed small and large split the power into the small trial and the",
  "large trial, and $E(N)$ is the expected total number of patients")
write_tex_table(mat1, head1, "rrrrrrrrrr",
                file.path(res_dir, "table1_operating_characteristics.tex"),
                cap1, "tab:operating",
                clines = c("\\cline{1-2}\\cline{3-5}\\cline{6-10}", ""))

## Table 2. The ceiling that the cutpoint places on the power.
t2 <- readRDS(file.path(data_dir, "ceiling_table.rds"))
mat2 <- cbind(fmt(t2$target, 2), fmt(t2$cutpoint, 3),
              fmt(t2$go_probability, 4), fmt(t2$power, 4),
              fmt(t2$power_small, 4), fmt(t2$ceiling, 4),
              ifelse(is.na(t2$required), "---",
                     formatC(t2$required, format = "d", big.mark = ",")))
head2 <- list(
  c("\\multicolumn{2}{c}{Cutpoint}", "\\multicolumn{4}{c}{Power under $H_1$}",
    ""),
  c("Target", "$C$", "$\\Pr(\\mathrm{expand})$", "Attained", "Small",
    "Limit", "Patients for 0.85"))
cap2 <- paste(
  "Power of the 2-in-1 adaptive design against the cutpoint, with the limit",
  "that the cutpoint places on it. The design is that of",
  "Table~\\ref{tab:operating} and the alternative sets the Phase 2 rate ratio",
  "to 0.65 and the Phase 3 rate ratio to 0.75. The cutpoint $C$ is chosen so",
  "that the probability of expansion attains the value in the column headed",
  "target. The column headed attained is the power at 750 patients per group in",
  "the large trial, the column headed small is the part of it contributed by",
  "the small trial, and the column headed limit is the power that would be",
  "reached if the large trial were made arbitrarily large. The last column is",
  "the total number of patients the large trial needs for a power of 0.85, and",
  "a dash marks a cutpoint whose limit lies below that target")
write_tex_table(mat2, head2, "rrrrrrr",
                file.path(res_dir, "table2_power_limit.tex"),
                cap2, "tab:limit",
                clines = c("\\cline{1-2}\\cline{3-6}", ""))

## Table 3. The normal approximation against simulation.
t3 <- readRDS(file.path(data_dir, "validation_table.rds"))
mat3 <- cbind(ifelse(duplicated(t3$label), "", t3$label), fmt(t3$cutpoint, 1),
              fmt(t3$approximate, 4), fmt(t3$known, 4), fmt(t3$se_known, 4),
              fmt(t3$estimated, 4), fmt(t3$se_estimated, 4))
head3 <- list(
  c("", "", "", "\\multicolumn{2}{c}{Dispersion known}",
    "\\multicolumn{2}{c}{Dispersion estimated}"),
  c("Scenario", "$C$", "Normal", "Simulated", "Standard error", "Simulated",
    "Standard error"))
cap3 <- paste(
  "Rejection probability of the 2-in-1 adaptive design obtained from the normal",
  "approximation and by simulation, for the design of",
  "Table~\\ref{tab:operating} and for a small design that randomises 30",
  "patients per group before the interim analysis, 45 per group in the small",
  "trial and 200 per group in the large trial. The column headed normal is the",
  "value returned by the closed-form evaluation. The simulated columns fit the",
  "negative binomial model to every generated trial, holding the dispersion at",
  "the design value in the first pair of columns and estimating it from the",
  "data in the second pair. The two pairs are computed on the same generated",
  "trials and differ only in the analysis. The global null sets both rate",
  "ratios to one and the remaining scenarios set the Phase 2 rate ratio to 0.65",
  "and the Phase 3 rate ratio to 0.75. Each simulated value uses",
  paste0(format(min(t3$nsim), big.mark = ",", scientific = FALSE),
         " replications for the reference design and"),
  paste0(format(max(t3$nsim), big.mark = ",", scientific = FALSE),
         " for the small design,"),
  "and the standard error is the binomial standard error at that number of",
  "replications")
write_tex_table(mat3, head3, "lrrrrrr",
                file.path(res_dir, "table3_validation.tex"),
                cap3, "tab:validation",
                clines = c("\\cline{4-5}\\cline{6-7}", ""))

message("manuscript exhibits written to ", res_dir)
