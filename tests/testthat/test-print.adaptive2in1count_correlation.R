quiet_print <- function(x, ...) {
  lines <- capture.output(returned <- withVisible(print(x, ...)))
  list(lines = lines, visible = returned$visible, value = returned$value)
}

base <- list(rate1 = 2.6, rate2 = 1.3, dispersion1 = 0.8, dispersion2 = 0.8,
             shared = 0.4, n1 = 60, n2 = 90, n3 = 300, window1 = 0.5,
             window2 = 1)

test_that("the correlation print method reports whether the requirement is met", {
  rho <- do.call(adaptive2in1count_correlation, base)
  out <- quiet_print(rho)
  expect_false(out$visible)
  expect_true(any(grepl("Correlation structure", out$lines, fixed = TRUE)))
  expect_true(any(grepl("is met", out$lines, fixed = TRUE)))
})
