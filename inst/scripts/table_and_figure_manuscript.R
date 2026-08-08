# Builds the exhibits of the main text from the data written by
# data_generate.R. Run from the package root. Figures are written as .pdf and
# .eps and tables as .tex to inst/scripts/results. The exhibits of the
# Supporting Information are built by the pair of scripts in
# inst/supplementary.

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

## Figure 1. Maximum type I error rate against the cutpoint.
d1 <- readRDS(file.path(data_dir, "typeI_curve.rds"))
d1$dispersion <- factor(d1$dispersion2, levels = c(0.5, 0.8, 1.5))
d1$row <- factor(d1$endpoint, levels = c("phase2", "phase3"),
                 labels = c("Small~trial~on~Phase~2",
                            "Small~trial~on~Phase~3"))
d1$col <- factor(d1$dropout, levels = c(0, 0.2, 0.4),
                 labels = c("No~dropout", "italic(d) == 0.2",
                            "italic(d) == 0.4"))
lab_disp2 <- parse(text = c("nu[2] == 0.5", "nu[2] == 0.8",
                            "nu[2] == 1.5"))
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
                 labels = c("Phase~3~ratio == 0.70",
                            "Phase~3~ratio == 0.80"))
d2$col <- factor(d2$ratio1, levels = c(0.55, 0.65, 0.75),
                 labels = c("Phase~2~ratio == 0.55",
                            "Phase~2~ratio == 0.65",
                            "Phase~2~ratio == 0.75"))
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

## Figure 3. Within-patient correlation of the two endpoints, by its two
## sources.
d3 <- readRDS(file.path(data_dir, "correlation_endpoint.rds"))
d3$frailty <- factor(d3$shared, levels = c(0, 0.23, 0.47, 0.7))
d3$row <- factor(d3$dispersion2, levels = c(0.8, 1.5),
                 labels = c("nu[2] == 0.8", "nu[2] == 1.5"))
d3$col <- factor(d3$window1, levels = c(0.25, 0.5, 1),
                 labels = c("tau[1] == 0.25", "tau[1] == 0.50",
                            "tau[1] == 1.00"))
lab_frailty <- parse(text = c("gamma == 0", "gamma == 0.23",
                              "gamma == 0.47", "gamma == 0.70"))
fig3 <- ggplot(d3, aes(dropout, correlation, colour = frailty,
                       linetype = frailty, shape = frailty)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3, data = subset(d3, round(dropout * 100) %% 10 == 0)) +
  facet_grid(row ~ col, labeller = label_parsed) +
  scale_colour_manual(values = pal4, labels = lab_frailty) +
  scale_linetype_manual(values = lty4, labels = lab_frailty) +
  scale_shape_manual(values = shape4, labels = lab_frailty) +
  labs(x = expression("Dropout proportion " * italic(d)),
       y = "Within-patient correlation of the two endpoints") +
  theme_ms
save_fig(fig3, "figure3_endpoint_correlation")

## Figure 4. Correlation of the interim statistic with the two final
## statistics under the two endpoint structures.
d4 <- readRDS(file.path(data_dir, "correlation_shared.rds"))
keys4 <- d4[c("n3", "dropout", "window1")]
d4_long <- rbind(
  data.frame(keys4, series = "shared_small", value = d4$shared_small),
  data.frame(keys4, series = "shared_large", value = d4$shared_large),
  data.frame(keys4, series = "distinct_small", value = d4$distinct_small),
  data.frame(keys4, series = "distinct_large", value = d4$distinct_large))
d4_long$series <- factor(d4_long$series,
                         levels = c("shared_small", "shared_large",
                                    "distinct_small", "distinct_large"))
d4_long$row <- factor(d4_long$window1, levels = c(0.5, 1),
                      labels = c("tau[1] == 0.50", "tau[1] == 1.00"))
d4_long$col <- factor(d4_long$dropout, levels = c(0, 0.2, 0.4),
                      labels = c("No~dropout", "italic(d) == 0.2",
                                 "italic(d) == 0.4"))
lab_series <- expression(rho[XY] * " one endpoint", rho[XZ] * " one endpoint",
                         rho[XY] * " two endpoints",
                         rho[XZ] * " two endpoints")
fig4 <- ggplot(d4_long, aes(n3, value, colour = series, linetype = series)) +
  geom_line(linewidth = 1.2) +
  facet_grid(row ~ col, labeller = label_parsed) +
  scale_colour_manual(values = pal4, labels = lab_series) +
  scale_linetype_manual(values = lty4, labels = lab_series) +
  expand_limits(y = 0) +
  guides(colour = guide_legend(nrow = 2), linetype = guide_legend(nrow = 2)) +
  labs(x = expression("Patients per group in the large trial, " * italic(n)[3]),
       y = "Correlation with the interim statistic") +
  theme_ms
save_fig(fig4, "figure4_endpoint_structure")

## Figure 5. Power of the tezepelumab design against the cutpoint.
d5 <- readRDS(file.path(data_dir, "real_power_curve.rds"))
d5$loss <- factor(d5$dropout, levels = c(0, 0.1, 0.2))
d5$row <- factor(d5$label, levels = c("Phase 3 assumptions",
                                      "Phase 2 assumptions"),
                 labels = c("Phase~3~assumptions", "Phase~2~assumptions"))
d5$col <- factor(d5$window1, levels = c(0.25, 0.5, 1),
                 labels = c("tau[1] == 0.25", "tau[1] == 0.50",
                            "tau[1] == 1.00"))
lab_loss5 <- parse(text = c("No~dropout", "italic(d) == 0.1",
                            "italic(d) == 0.2"))
fig5 <- ggplot(d5, aes(cutpoint, rejection, colour = loss, linetype = loss)) +
  reference_line(0.9) +
  geom_line(linewidth = 1.2) +
  facet_grid(row ~ col, labeller = label_parsed) +
  scale_colour_manual(values = pal3, labels = lab_loss5) +
  scale_linetype_manual(values = lty3, labels = lab_loss5) +
  scale_cutpoint +
  labs(x = expression("Cutpoint " * italic(C) *
                        " for the observed rate ratio"),
       y = "Power") +
  theme_ms
save_fig(fig5, "figure5_real_trial")

## Table 1. Operating characteristics at calibrated cutpoints.
t1 <- readRDS(file.path(data_dir, "operating_table.rds"))
mat1 <- cbind(fmt(t1$target, 2), fmt(t1$cutpoint, 3), fmt(t1$go_null, 3),
              fmt(t1$type1, 4), fmt(t1$size_null, 0),
              fmt(t1$go_alternative, 2), fmt(t1$power, 3),
              fmt(t1$power_small, 3), fmt(t1$power_large, 3),
              fmt(t1$size_alternative, 0))
head1 <- list(
  c("\\multicolumn{2}{c}{Cutpoint}", "\\multicolumn{3}{c}{Under $H_0$}",
    "\\multicolumn{5}{c}{Under $H_1$}"),
  c("Target", "$C$", "$\\Pr(\\mathrm{expand})$", "Type I", "$E(N)$",
    "$\\Pr(\\mathrm{expand})$", "Power", "Small", "Large", "$E(N)$"))
cap1 <- paste(
  "Operating characteristics of the 2-in-1 adaptive design for count endpoints",
  "in the reference chronic obstructive pulmonary disease setting. The Phase 2",
  "endpoint has rate $\\theta_1 = 2$ events per year, dispersion",
  "$\\nu_1 = 0.7$ and observation window $\\tau_1 = 0.5$ years, and the",
  "Phase 3 endpoint has rate $\\theta_2 = 1.3$ events per year, dispersion",
  "$\\nu_2 = 0.8$ and observation window $\\tau_2 = 1$ year. The shape of the",
  "frailty component shared by the two endpoints is $\\gamma = 0.4$. The",
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
                clines = c("\\cline{1-2}\\cline{3-5}\\cline{6-10}", ""),
                size = "\\footnotesize", colsep = "3pt")

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
    "Limit", "$N$ for 0.85"))
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
  "the total number of patients $N$ the large trial needs for a power of 0.85,",
  "and a dash marks a cutpoint whose limit lies below that target")
write_tex_table(mat2, head2, "rrrrrrr",
                file.path(res_dir, "table2_power_limit.tex"),
                cap2, "tab:limit",
                clines = c("\\cline{1-2}\\cline{3-6}", ""))

## Table 3. The normal approximation against simulation.
t3 <- readRDS(file.path(data_dir, "validation_table.rds"))
# The scenario labels carry the design, the hypothesis and the dropout
# proportion, which is more compact than a sentence and keeps the table inside
# the text width.
short3 <- c(
  "Reference design, global null" = "Reference, $H_0$, $d = 0$",
  "Reference design, global null, dropout 0.3" = "Reference, $H_0$, $d = 0.3$",
  "Reference design, both endpoints, dropout 0.3" =
    "Reference, $H_1$, $d = 0.3$",
  "Reference design, both endpoints, dropout 0.5" =
    "Reference, $H_1$, $d = 0.5$",
  "Small design, global null, dropout 0.3" = "Small, $H_0$, $d = 0.3$",
  "Small design, both endpoints, dropout 0.3" = "Small, $H_1$, $d = 0.3$",
  "Shared endpoint, global null, dropout 0.3" = "Shared, $H_0$, $d = 0.3$",
  "Shared endpoint, both windows, dropout 0.3" =
    "Shared, $H_1$, $d = 0.3$, $\\tau_1 = 0.5$",
  "Shared endpoint, one window, dropout 0.3" =
    "Shared, $H_1$, $d = 0.3$, $\\tau_1 = 1$")
stopifnot(all(t3$label %in% names(short3)))
label3 <- unname(short3[t3$label])
mat3 <- cbind(ifelse(duplicated(label3), "", label3), fmt(t3$cutpoint, 1),
              fmt(t3$approximate, 4), fmt(t3$known, 4), fmt(t3$se_known, 4),
              fmt(t3$estimated, 4), fmt(t3$se_estimated, 4))
head3 <- list(
  c("", "", "", "\\multicolumn{2}{c}{Dispersion known}",
    "\\multicolumn{2}{c}{Dispersion estimated}"),
  c("Scenario", "$C$", "Normal", "Simulated", "SE", "Simulated", "SE"))
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
  "trials and differ only in the analysis. The rows headed reference and small",
  "have two endpoints that count different events, and the rows headed shared",
  "have one endpoint counted over the interim window $\\tau_1$ and over the",
  "year, with rate 1.3 events per year and dispersion 0.8. $H_0$ sets the rate",
  "ratios to one, $H_1$ sets the Phase 2 rate ratio to 0.65 and the Phase 3",
  "rate ratio to 0.75 for two endpoints and sets the single rate ratio to 0.75",
  "for a shared endpoint. Each simulated value uses",
  paste0(format(min(t3$nsim), big.mark = ",", scientific = FALSE),
         " replications for the reference design and"),
  paste0(format(max(t3$nsim), big.mark = ",", scientific = FALSE),
         " for the small design,"),
  "and the standard error is the binomial standard error at that number of",
  "replications")
write_tex_table(mat3, head3, "lrrrrrr",
                file.path(res_dir, "table3_validation.tex"),
                cap3, "tab:validation",
                clines = c("\\cline{4-5}\\cline{6-7}", ""),
                size = "\\footnotesize", colsep = "3pt")

## Table 4. The tezepelumab design at calibrated cutpoints.
t4 <- readRDS(file.path(data_dir, "real_table.rds"))
mat4 <- cbind(ifelse(duplicated(t4$label), "", t4$label), fmt(t4$target, 2),
              fmt(t4$cutpoint, 3), fmt(t4$type1, 4), fmt(t4$power, 3),
              fmt(t4$power_small, 3), fmt(t4$ceiling, 3),
              ifelse(is.na(t4$required), "---",
                     formatC(t4$required, format = "d", big.mark = ",")),
              fmt(t4$size_alternative, 0))
head4 <- list(
  c("", "\\multicolumn{2}{c}{Cutpoint}", "",
    "\\multicolumn{4}{c}{Under $H_1$}", ""),
  c("Assumptions", "Target", "$C$", "Type I", "Power", "Small", "Limit",
    "$N$ for 0.9", "$E(N)$"))
cap4 <- paste(
  "The tezepelumab programme in severe asthma evaluated as a single 2-in-1",
  "trial. The annualized rate of asthma exacerbations is read at an interim",
  "analysis over half a year and at the final analysis over a year, so one",
  "endpoint serves all three analyses. The Phase 3 assumptions are those of the",
  "protocol of the confirmatory trial, a placebo rate of 0.9 events per year, a",
  "dispersion parameter of 2.4 in the sense of Tang and of Zhu and Lakkis and",
  "therefore $\\nu = 1 / 2.4$ here, a rate ratio of 0.5 and a dropout",
  "proportion of 0.1 spread uniformly over the year. The Phase 2 assumptions",
  "are those of the earlier trial, a placebo rate of 0.7 events per year,",
  "$\\nu = 1 / 0.7$ and a rate ratio of 0.6. The trial randomises 138 patients",
  "per group before the interim analysis and in the small trial and 530 per",
  "group in the large trial, and the Wald test is applied at a one-sided",
  "significance level of 0.025. The cutpoint $C$ is chosen so that the",
  "probability of expansion attains the value in the column headed target. The",
  "column headed limit is the power that an arbitrarily large trial would",
  "reach, the column headed $N$ for 0.9 is the total number of patients the",
  "large trial needs for a power of 0.9, and a dash marks a cutpoint whose",
  "limit lies below that target")
write_tex_table(mat4, head4, "lrrrrrrrr",
                file.path(res_dir, "table4_real_trial.tex"),
                cap4, "tab:real",
                clines = c("\\cline{2-3}\\cline{5-8}", ""),
                size = "\\footnotesize", colsep = "3pt")

message("exhibits of the main text written to ", res_dir)
