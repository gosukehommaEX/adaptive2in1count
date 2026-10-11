test_that("the correlation print method reports whether the requirement is met", {
  rho <- do.call(adaptive2in1count_correlation, base)
  out <- quiet_print(rho)
  expect_false(out$visible)
  expect_true(any(grepl("Correlation structure", out$lines, fixed = TRUE)))
  expect_true(any(grepl("is met", out$lines, fixed = TRUE)))
})
