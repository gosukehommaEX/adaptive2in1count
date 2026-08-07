# adaptive2in1count

<!-- badges: start -->
<!-- badges: end -->

Design and evaluation of 2-in-1 adaptive trials whose Phase 2 and Phase 3 parts
use overdispersed count endpoints, allowing patients to leave before the end of
the follow-up period.

## The design

The trial starts as a small trial. At the interim analysis the observed rate
ratio of the Phase 2 endpoint is compared with a cutpoint. If the observed rate
ratio is at or below the cutpoint the trial expands into a large trial and the
final analysis uses the Phase 3 endpoint on all patients, otherwise the trial
remains small and the final analysis uses the endpoint named by
`small_endpoint`. Both final analyses are Wald tests of the log rate ratio from
a negative binomial model with an offset for the realized exposure, carried out
at the full one-sided level.

## The model

Counts follow an additive gamma frailty model. Patient `i` carries a frailty
component shared by the two endpoints, with shape `shared`, and a component
specific to each endpoint, so each endpoint is marginally negative binomial with
its own dispersion. A single follow-up period determines the exposure of both
endpoints, `min(T, window1)` and `min(T, window2)`. The period is either
complete, a two point mixture in which a proportion `dropout` leave at a time
uniform on the follow-up window, or an exponential variable truncated at that
window with the same `dropout` as the proportion who leave, following Tang
(2015) and Homma and Daimon (2021). Both random forms reduce to complete
follow-up when `dropout` is zero.

The two endpoints of a patient are therefore correlated through two mechanisms.
The shared frailty correlates the underlying event rates, and a random follow-up
period correlates the exposures because a patient who leaves stops contributing
to either count. The second mechanism operates even when the frailty is not
shared, so ignoring dropout understates the correlation that drives the interim
decision.

## Installation

``` r
# install.packages("pak")
pak::pak("gosukehommaEX/adaptive2in1count")
```

## Example

An asthma or chronic obstructive pulmonary disease design in which the Phase 2
part counts exacerbations over half a year and the Phase 3 part is the
annualized exacerbation rate.

``` r
library(adaptive2in1count)

adaptive2in1count(
  rate1 = 2.6, rate2 = 1.3,
  dispersion1 = 0.8, dispersion2 = 0.8, shared = 0.4,
  n1 = 60, n2 = 90, n3 = 300,
  cutpoint = c(0.6, 0.8, 1.0, 1.4),
  window1 = 0.5, window2 = 1,
  dropout = 0.3, dropout_shape = "exponential"
)
```

The type I error rate is maximised over a grid of plausible rates, because the
rates enter the information per patient and therefore the correlation between
the interim and the final statistics.

``` r
adaptive2in1count_typeI(
  rate1 = c(1.5, 2.6, 4), rate2 = c(0.8, 1.3, 2),
  dispersion1 = 0.8, dispersion2 = 0.8, shared = 0.4,
  n1 = 60, n2 = 90, n3 = 300,
  cutpoint = c(0.4, 0.6, 0.8, 1, 1.4, 2),
  window1 = 0.5, window2 = 1, dropout = 0.3,
  dropout_shape = "exponential"
)
```

Power is reported with its split into the small trial and the large trial,
because the endpoint tested at the final analysis depends on the interim
decision.

``` r
adaptive2in1count_power(
  rate1 = 2.6, rate2 = 1.3,
  dispersion1 = 0.8, dispersion2 = 0.8, shared = 0.4,
  n1 = 60, n2 = 90, n3 = 300, cutpoint = 1,
  ratio1 = c(0.5, 0.65, 0.8), ratio2 = c(0.7, 0.85),
  window1 = 0.5, window2 = 1, dropout = 0.3,
  dropout_shape = "exponential"
)
```

Enlarging the large trial only helps when the trial expands, so the power
approaches the probability of expansion plus the contribution of the small
trial. `adaptive2in1count_size()` reports that limit when a target is beyond it.

``` r
adaptive2in1count_size(
  rate1 = 2.6, rate2 = 1.3,
  dispersion1 = 0.8, dispersion2 = 0.8, shared = 0.4,
  n1 = 60, n2 = 90, cutpoint = 1.4, target = 0.9,
  ratio1 = 0.65, ratio2 = 0.75,
  window1 = 0.5, window2 = 1, dropout = 0.3,
  dropout_shape = "exponential"
)
```

Type I error control of the 2-in-1 design requires the interim statistic to be
at least as strongly correlated with the small trial statistic as with the large
trial statistic. Both correlations are available in closed form.

``` r
adaptive2in1count_correlation(
  rate1 = 2.6, rate2 = 1.3,
  dispersion1 = 0.8, dispersion2 = 0.8, shared = 0.4,
  n1 = 60, n2 = 90, n3 = 300,
  window1 = 0.5, window2 = 1, dropout = 0.3,
  dropout_shape = "exponential"
)
```

## Checking the normal approximation

`adaptive2in1count_simulate()` generates trials from the model and analyses each
of them as a trial would be analysed. Three analyses are available.

| `method` | Analysis | Purpose |
| --- | --- | --- |
| `"estimated"` | Maximum likelihood fit of the two rates and the common dispersion, vectorised over replications | The analysis a trial would report |
| `"known"` | The same fit with the dispersion held at the design value | Isolates the error of the normal approximation |
| `"glm"` | `MASS::glm.nb`, one replication at a time | Reference implementation |

Comparing `"known"` with `"estimated"` separates the error of the normal
approximation from the effect of estimating the dispersion.

``` r
adaptive2in1count_simulate(
  rate1 = 2.6, rate2 = 1.3,
  dispersion1 = 0.8, dispersion2 = 0.8, shared = 0.4,
  n1 = 60, n2 = 90, n3 = 300, cutpoint = 1,
  window1 = 0.5, window2 = 1, dropout = 0.3,
  dropout_shape = "exponential",
  ratio1 = 0.65, ratio2 = 0.75, nsim = 10000, seed = 1
)
```

## Main functions

| Function | Purpose |
| --- | --- |
| `adaptive2in1count()` | Probability of expansion, rejection probability and expected sample size for one or more cutpoints |
| `adaptive2in1count_typeI()` | Maximum type I error rate over a grid of nuisance rates |
| `adaptive2in1count_power()` | Power over a grid of true rate ratios |
| `adaptive2in1count_size()` | Size of the large trial for a target power |
| `adaptive2in1count_cutpoint()` | Cutpoint attaining a target probability of expansion |
| `adaptive2in1count_correlation()` | Correlations between the interim and the two final statistics |
| `adaptive2in1count_simulate()` | Simulated operating characteristics |
| `design_components()` | Variances, means and correlations of the three test statistics |
| `follow_up_moments()` | Moments of the exposure of the two endpoints |
| `information_factor()` | Information per patient for a negative binomial rate |
| `endpoint_correlation()` | Within-patient correlation of the two endpoints |
| `negative_binomial_contrast()` | Log rate ratio of two groups, fitted over many replications at once |

## Complete follow-up

`adaptive2in1count_exact()` evaluates the same design without the normal
approximation, by propagating the joint distribution of the interim total and
the final total. It applies only to complete follow-up, and its final analyses
may be an exact conditional test that is free of the nuisance rate or a Wald
test.

## References

Chen C, Anderson K, Mehrotra DV, Rubin EH, Tse A. A 2-in-1 adaptive phase 2/3
design for expedited oncology drug development. *Contemporary Clinical Trials*
2018; 64: 238-242. <https://doi.org/10.1016/j.cct.2017.09.006>

Tang Y. Sample size estimation for negative binomial regression comparing rates
of recurrent events with unequal follow-up time. *Journal of Biopharmaceutical
Statistics* 2015; 25: 1100-1113.
<https://doi.org/10.1080/10543406.2014.971167>

Homma G, Daimon T. Sample size calculation for "gold-standard" noninferiority
trials with fixed margins and negative binomial endpoints. *Statistics in
Biopharmaceutical Research* 2021.
