# adaptive2in1count 0.1.0

* First release.
* `adaptive2in1count()` evaluates the probability of expansion, the rejection
  probability and the expected sample size of a 2-in-1 adaptive design with
  overdispersed count endpoints in closed form, under complete follow-up or a
  uniform or exponential follow-up period.
* Two endpoint structures through `endpoints`: `"distinct"` for two endpoints
  that count different events (additive gamma frailty) and `"shared"` for one
  endpoint read over two windows.
* `adaptive2in1count_typeI()`, `adaptive2in1count_power()`,
  `adaptive2in1count_size()`, `adaptive2in1count_cutpoint()` and
  `adaptive2in1count_correlation()` for the maximum type I error rate, power
  grids, the size of the large trial, calibrated cutpoints and the correlation
  condition.
* `adaptive2in1count_simulate()` with `method = "estimated"`, `"known"` or
  `"glm"`, and `adaptive2in1count_exact()` for complete follow-up.
* `inst/scripts` and `inst/supplementary` reproduce every figure and table of
  the article and its Supporting Information, and `inst/reproduce` reproduces
  sample sizes and power published by Tang (2015), Corren et al. (2017) and
  Menzies-Gow et al. (2021).
