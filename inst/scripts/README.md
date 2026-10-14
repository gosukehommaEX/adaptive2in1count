# Scripts for the main text

Run both scripts from the package root, in the order below. The exhibits of the
Supporting Information are produced by the pair of scripts in
`inst/supplementary`, which are independent of these.

1. `data_generate.R` evaluates the operating characteristics and writes one
   `.rds` file per exhibit to `data/`.
2. `table_and_figure_manuscript.R` reads those files and writes the figures as
   `.pdf` and `.eps` and the tables as `.tex` to `results/`.

Separating the two steps means that a change to a figure or a table does not
require the operating characteristics to be evaluated again. Steps 7, 8 and 8b
of `data_generate.R` simulate trials, one hundred thousand per scenario under
the global null hypothesis, and take many hours; the rest of the script takes
seconds. Their seeds are fixed, so by default the script reads their saved
results back from `data/validation_table.rds` and `data/shared_rerun.rds`.
Setting `rerun_simulations <- TRUE` at the top of the script simulates them
again.

## Exhibits

| File | Content |
| --- | --- |
| `figure1_type_one_error` | Maximum type I error rate against the cutpoint, by dropout, dispersion and small trial endpoint |
| `figure2_power` | Power against the cutpoint, by treatment effect and dropout |
| `figure3_endpoint_correlation` | Within-patient correlation of the two counts, separated into its frailty and follow-up sources, against the correlation of the two estimated log rate ratios |
| `figure4_endpoint_structure` | Correlation of the interim statistic with the two final statistics under both endpoint structures |
| `figure5_real_trial` | Power of the tezepelumab design against the cutpoint, by interim window, dropout and set of design assumptions |
| `table1_operating_characteristics.tex` | Operating characteristics at cutpoints calibrated to a target probability of expansion |
| `table2_power_limit.tex` | Power against the cutpoint, with the limit the cutpoint places on it |
| `table3_validation.tex` | Normal approximation against simulation, for a reference design, a small design and a shared endpoint |
| `table4_real_trial.tex` | The tezepelumab design at calibrated cutpoints, under the assumptions of each of the two trials |
| `manuscript_numbers.txt` | Every number the main text quotes, written as plain text |

Step 10 writes `manuscript_numbers.txt`, which collects every number the main
text quotes in the order in which the text uses them. Checking a sentence
against that file is the same as checking it against the evaluation, so no
number in the manuscript rests on a transcription. Statements about the shape
of a curve (a rise, a fall, a trough, or one curve lying above another) are
checked there on grids finer than the ones the figures are drawn from:
cutpoints in steps of 0.001 for Figures 2 and 5, and dropout proportions in
steps of 0.0025 for Figure 3.

Each `.tex` file is a complete `table` environment that uses only standard
LaTeX commands, so it can be uploaded to the manuscript source without further
editing.

## Figures and type size

Figures are written twelve inches wide with twenty point type, so that
`\includegraphics[width=\textwidth]` scales them to about six inches and leaves
roughly ten point type on the page, which matches the body text. Changing
`save_fig()` therefore requires `base_size` in `theme_ms` to be changed in the
same proportion, otherwise the printed type no longer matches the body text.

## The dispersion parameter

The package and the manuscript write the negative binomial as in Homma and
Daimon (2019), so `dispersion` is the shape of the gamma frailty and the
variance of a count with mean `mu` is `mu + mu^2 / dispersion`. This is the
`size` argument of `rnbinom()`. Tang (2015) and Zhu and Lakkis (2014) instead
report the variance of the frailty, which is the reciprocal. A protocol that
states a dispersion of 2.4 in their sense therefore enters these scripts as
`dispersion = 1 / 2.4`, which step 9 of `data_generate.R` does for the
tezepelumab design.
