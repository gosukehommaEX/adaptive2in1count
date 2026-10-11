# Reproduction of published numbers

`reproduce_published.R` recomputes numbers that other authors have published,
using the functions of the package, and writes a one-line verdict for each to
`results/reproduce_summary.txt`. Run it from the package root.

| Source | Numbers | What they check |
| --- | --- | --- |
| Tang (2015), Tables 1 and 2 | 28 sample sizes per group | `information_factor()` and the uniform and exponential follow-up distributions |
| Corren et al. (2017) | 138 patients per group | the dispersion of the Phase 2 tezepelumab trial is the variance of the frailty |
| Menzies-Gow et al. (2021), protocol | power above 0.99, minimum detectable reduction of 27 per cent | the dispersion of the Phase 3 tezepelumab trial is the variance of the frailty |

The verdicts are PASS when the published number is reproduced exactly at the
precision the source reports, FAIL when it is not, and INFO for a line that
reports a contrasting calculation, such as the other reading of a dispersion
parameter. EXPLAINED is reserved for a published number that cannot be
reproduced for a documented reason; none is needed at present.

The summary file is written by the script and is not edited by hand.
