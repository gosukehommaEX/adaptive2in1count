# Manuscript scripts

Run both scripts from the package root, in the order below.

1. `data_generate.R` evaluates the operating characteristics and writes one
   `.rds` file per exhibit to `data/`.
2. `table_and_figure_manuscript.R` reads those files and writes the figures as
   `.pdf` and `.eps` and the tables as `.tex` to `results/`.

Separating the two steps means that a change to a figure or a table does not
require the operating characteristics to be evaluated again. Step 8 of
`data_generate.R` simulates ten thousand trials for each of six scenarios and
takes several minutes; the rest of the script takes seconds.

## Exhibits

| File | Content |
| --- | --- |
| `figure1_type_one_error` | Maximum type I error rate against the cutpoint, by dropout, dispersion and small trial endpoint |
| `figure2_power` | Power against the cutpoint, by treatment effect and dropout |
| `figure3_required_size` | Patients required in the large trial and expected sample size against the dropout proportion |
| `figure4_endpoint_correlation` | Within-patient correlation of the two endpoints, separated into its frailty and follow-up sources |
| `figure5_statistic_correlation` | Correlation of the interim statistic with the two final statistics |
| `table1_operating_characteristics.tex` | Operating characteristics at cutpoints calibrated to a target probability of expansion |
| `table2_power_limit.tex` | Power against the cutpoint, with the limit the cutpoint places on it |
| `table3_validation.tex` | Normal approximation against simulation, for a reference and a small design, with the dispersion held fixed and estimated |

Each `.tex` file is a complete `table` environment that uses only standard
LaTeX commands, so it can be uploaded to the manuscript source without further
editing.
