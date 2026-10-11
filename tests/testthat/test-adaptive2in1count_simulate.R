test_that("the simulation agrees with the normal approximation", {
  arguments <- list(rate1 = 2.6, rate2 = 1.3, dispersion1 = 0.8,
                    dispersion2 = 0.8, shared = 0.4, n1 = 60, n2 = 90,
                    n3 = 200, cutpoint = c(0.6, 1), window1 = 0.5,
                    window2 = 1)
  approximate <- do.call(adaptive2in1count, arguments)$results
  simulated <- do.call(adaptive2in1count_simulate,
                       c(arguments, list(method = "known", nsim = 3000,
                                 seed = 20260806)))
  gap <- abs(simulated$go_probability - approximate$go_probability) /
    pmax(simulated$se_go, 1e-12)
  expect_true(all(gap < 5))
  expect_true(all(simulated$failures == 0L))
})

test_that("the simulation agrees under dropout and an alternative", {
  arguments <- list(rate1 = 2.6, rate2 = 1.3, dispersion1 = 0.8,
                    dispersion2 = 0.8, shared = 0.4, n1 = 60, n2 = 90,
                    n3 = 200, cutpoint = 0.8, window1 = 0.5, window2 = 1,
                    dropout = 0.3, dropout_shape = "exponential",
                    ratio1 = 0.6, ratio2 = 0.7)
  approximate <- do.call(adaptive2in1count, arguments)$results
  simulated <- do.call(adaptive2in1count_simulate,
                       c(arguments, list(method = "estimated", nsim = 1000,
                                 seed = 5)))
  expect_lt(abs(simulated$rejection - approximate$rejection) /
              simulated$se_rejection, 5)
})

test_that("the simulated probability of expansion tracks the cutpoint", {
  out <- adaptive2in1count_simulate(rate1 = 2.6, rate2 = 1.3,
                                    dispersion1 = 0.8, dispersion2 = 0.8,
                                    shared = 0.4, n1 = 60, n2 = 90, n3 = 150,
                                    cutpoint = c(0.5, 0.8, 1.5),
                                    window1 = 0.5, window2 = 1, method = "known",
                                    nsim = 1000, seed = 3)
  expect_equal(nrow(out), 3L)
  expect_true(all(diff(out$go_probability) >= 0))
  expect_true(all(out$se_rejection >= 0))
})

test_that("invalid arguments are rejected", {
  expect_error(adaptive2in1count_simulate(2.6, 1.3, 0.8, 0.8, 0.4, 200, 90,
                                          300, 0.8, nsim = 5),
               "non-decreasing")
  expect_error(adaptive2in1count_simulate(2.6, 1.3, 0.8, 0.8, 2, 60, 90, 300,
                                          0.8, nsim = 5),
               "smaller of the two dispersions")
})
