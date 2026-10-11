quiet_print <- function(x, ...) {
  lines <- capture.output(returned <- withVisible(print(x, ...)))
  list(lines = lines, visible = returned$visible, value = returned$value)
}

base <- list(rate1 = 2.6, rate2 = 1.3, dispersion1 = 0.8, dispersion2 = 0.8,
             shared = 0.4, n1 = 60, n2 = 90, n3 = 300, window1 = 0.5,
             window2 = 1)

test_that("the design print method returns its input invisibly", {
  fit <- do.call(adaptive2in1count, c(base, list(cutpoint = c(0.6, 1))))
  out <- quiet_print(fit)
  expect_false(out$visible)
  expect_identical(out$value, fit)
  expect_true(any(grepl("2-in-1 adaptive design for count endpoints",
                        out$lines, fixed = TRUE)))
  expect_true(any(grepl("Interim correlation", out$lines, fixed = TRUE)))
  expect_true(any(grepl("P(expand)", out$lines, fixed = TRUE)))
})

test_that("the follow-up assumption is reported", {
  fit <- do.call(adaptive2in1count,
                 c(base, list(cutpoint = 0.8, dropout = 0.3,
                              dropout_shape = "exponential")))
  out <- quiet_print(fit)
  expect_true(any(grepl("exponential", out$lines, fixed = TRUE)))
  expect_true(any(grepl("Small trial endpoint  Phase 2", out$lines,
                        fixed = TRUE)))
})

test_that("the digits argument changes the printed precision", {
  fit <- do.call(adaptive2in1count, c(base, list(cutpoint = 0.8)))
  expect_false(identical(quiet_print(fit, digits = 2)$lines,
                         quiet_print(fit, digits = 7)$lines))
})
