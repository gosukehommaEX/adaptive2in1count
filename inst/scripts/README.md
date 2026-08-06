# Manuscript scripts

`data_generate.R` will evaluate the operating characteristics and write one
`.rds` file per exhibit to `data/`, and `table_and_figure_manuscript.R` will
read those files and write the figures as `.pdf` and `.eps` and the tables as
`.tex` to `results/`. Separating the two steps means that a change to a figure
or a table does not require the operating characteristics to be evaluated
again.

Both scripts are being rewritten for the design that allows patients to leave
before the end of the follow-up period, and are not yet present.
