test_that("the type I error print method reports the overall maximum", {
  size <- adaptive2in1count_typeI(rate1 = c(2, 2.6), rate2 = 1.3,
                                  dispersion1 = 0.8, dispersion2 = 0.8,
                                  shared = 0.4, n1 = 60, n2 = 90, n3 = 300,
                                  cutpoint = c(0.6, 1), window1 = 0.5,
                                  window2 = 1)
  out <- quiet_print(size)
  expect_false(out$visible)
  expect_true(any(grepl("Maximum type I error rate", out$lines, fixed = TRUE)))
  expect_true(any(grepl("Overall maximum", out$lines, fixed = TRUE)))
})
