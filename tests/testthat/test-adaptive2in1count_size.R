base <- list(rate1 = 2.6, rate2 = 1.3, dispersion1 = 0.8, dispersion2 = 0.8,
             shared = 0.4, n1 = 60, n2 = 90, cutpoint = 1.4, ratio1 = 0.65,
             ratio2 = 0.75, window1 = 0.5, window2 = 1)

test_that("the solved size attains the target power", {
  for (target in c(0.8, 0.9)) {
    out <- do.call(adaptive2in1count_size, c(base, list(target = target)))
    expect_gte(out$power, target)
    expect_true(all(out$n3 >= 90))
  }
})

test_that("a higher target needs a larger trial", {
  low <- do.call(adaptive2in1count_size, c(base, list(target = 0.8)))
  high <- do.call(adaptive2in1count_size, c(base, list(target = 0.95)))
  expect_gt(sum(high$n3), sum(low$n3))
})

test_that("dropout needs a larger trial", {
  complete <- do.call(adaptive2in1count_size, c(base, list(target = 0.9)))
  reduced <- do.call(adaptive2in1count_size,
                     c(base, list(target = 0.9, dropout = 0.4,
                                  dropout_shape = "exponential")))
  expect_gt(sum(reduced$n3), sum(complete$n3))
})

test_that("unequal allocation is returned in the requested ratio", {
  out <- do.call(adaptive2in1count_size,
                 c(base, list(target = 0.9, allocation = c(2, 1))))
  expect_equal(out$n3[1] / out$n3[2], 2, tolerance = 0.05)
  expect_gte(out$power, 0.9)
})

test_that("the expected sample size lies between the two trial sizes", {
  out <- do.call(adaptive2in1count_size, c(base, list(target = 0.9)))
  expect_gte(out$expected_size, 180)
  expect_lte(out$expected_size, sum(out$n3))
})

test_that("an unattainable target is reported with the limit", {
  expect_error(do.call(adaptive2in1count_size,
                       modifyList(base, list(target = 0.99, cutpoint = 0.6))),
               "not attainable")
})

test_that("invalid arguments are rejected", {
  expect_error(do.call(adaptive2in1count_size, c(base, list(target = 0))),
               "between 0 and 1")
  expect_error(do.call(adaptive2in1count_size,
                       modifyList(base, list(target = 0.9,
                                             cutpoint = c(0.5, 1)))),
               "single value")
})
