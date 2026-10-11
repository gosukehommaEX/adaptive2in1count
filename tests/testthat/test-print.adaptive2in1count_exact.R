quiet_print <- function(x, ...) {
  lines <- capture.output(returned <- withVisible(print(x, ...)))
  list(lines = lines, visible = returned$visible, value = returned$value)
}

test_that("the exact print method still works", {
  fit <- adaptive2in1count_exact(rate1 = 2, rate2 = 0.5, dispersion1 = 0.6,
                                 dispersion2 = 1.2, shared = 0.4, n1 = 20,
                                 n2 = 30, n3 = 60, cutpoint = 0.6)
  out <- quiet_print(fit)
  expect_false(out$visible)
  expect_true(any(grepl("exact evaluation", out$lines, fixed = TRUE)))
})
