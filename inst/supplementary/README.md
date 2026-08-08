# Scripts for the Supporting Information

Run both scripts from the package root, in the order below. The exhibits of the
main text are produced by the pair of scripts in `inst/scripts`, which are
independent of these.

1. `data_generate.R` evaluates the operating characteristics and writes one
   `.rds` file per exhibit to `data/`.
2. `table_and_figure_supplementary.R` reads those files and writes the figures
   as `.pdf` and `.eps` and the table as `.tex` to `results/`.

Step 3 of `data_generate.R` generates fifteen thousand replicate trials for
each of eight settings and takes a few minutes; the rest of the script takes
seconds.

## Exhibits

| File | Content |
| --- | --- |
| `figureS1_required_size` | Patients required in the large trial and expected sample size against the dropout proportion |
| `figureS2_statistic_correlation` | Correlation of the interim statistic with the two final statistics, for two endpoints that count different events |
| `tableS1_correlation_check.tex` | Closed-form correlations of the three statistics against the correlations of statistics computed from generated patients |

Table S1 measures the approximation that the closed form carries when the
follow-up period is random. The closed form treats the total exposure as fixed
at its expectation, which is exact under complete follow-up and an
approximation otherwise, and the table reports the size of the gap rather than
bounding it.

## Conventions

The design parameters at the head of `data_generate.R` are the same as those of
the main text, so the two scripts describe one reference design. The
`dispersion` argument is the shape of the gamma frailty, as in Homma and Daimon
(2019) and as in the `size` argument of `rnbinom()`, and is the reciprocal of
the dispersion parameter reported by Tang (2015) and by Zhu and Lakkis (2014).

Figures are written twelve inches wide with twenty point type, so that
`\includegraphics[width=\textwidth]` scales them to about six inches and leaves
roughly ten point type on the page.
