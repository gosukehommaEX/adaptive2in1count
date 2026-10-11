quiet_print <- function(x, ...) {
  lines <- capture.output(returned <- withVisible(print(x, ...)))
  list(lines = lines, visible = returned$visible, value = returned$value)
}

test_that("the power print method reports the design inputs", {
  curve <- adaptive2in1count_power(rate1 = 2.6, rate2 = 1.3, dispersion1 = 0.8,
                                   dispersion2 = 0.8, shared = 0.4, n1 = 60,
                                   n2 = 90, n3 = 300, cutpoint = 0.8,
                                   ratio1 = 0.65, ratio2 = 0.75,
                                   window1 = 0.5, window2 = 1)
  out <- quiet_print(curve)
  expect_false(out$visible)
  expect_true(any(grepl("Power of a 2-in-1 count design", out$lines,
                        fixed = TRUE)))
})
