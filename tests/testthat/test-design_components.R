reference <- list(rate1 = 2.6, rate2 = 1.3, dispersion1 = 0.8,
                  dispersion2 = 0.8, shared = 0.4, n1 = 60, n2 = 90, n3 = 300,
                  window1 = 0.5, window2 = 1)

test_that("the variances match the closed form under complete follow-up", {
  out <- do.call(design_components, reference)
  factor1 <- information_factor(2.6, 0.8, 0.5, 1)
  factor2 <- information_factor(1.3, 0.8, 1, 1)
  expect_equal(out$variance[["interim"]], 2 / (60 * factor1), tolerance = 1e-9)
  expect_equal(out$variance[["small"]], 2 / (90 * factor1), tolerance = 1e-9)
  expect_equal(out$variance[["large"]], 2 / (300 * factor2), tolerance = 1e-9)
})

test_that("a Phase 2 small trial gives the classical nesting correlation", {
  out <- do.call(design_components, reference)
  expect_equal(out$correlation[["interim_small"]], sqrt(60 / 90),
               tolerance = 1e-9)
})

test_that("the coupling reduces to the frailty term under complete follow-up", {
  out <- do.call(design_components, reference)
  expect_equal(out$coupling, rep(0.4 / (0.8 * 0.8), 2), tolerance = 1e-9)
})

test_that("dropout raises the coupling above the frailty term", {
  out <- do.call(design_components, c(reference,
                                      list(dropout = 0.4,
                                           dropout_shape = "exponential")))
  expect_true(all(out$coupling > 0.4 / (0.8 * 0.8)))
})

test_that("the coupling vanishes without a shared frailty under dropout", {
  out <- do.call(design_components,
                 c(modifyList(reference, list(shared = 0)),
                   list(dropout = 0.4, dropout_shape = "exponential")))
  expect_equal(out$coupling, c(0, 0))
  expect_equal(out$correlation[["interim_large"]], 0)
})

test_that("the coupling of one endpoint gives the canonical covariance", {
  out <- do.call(design_components,
                 c(reference, list(dropout = 0.4,
                                   dropout_shape = "exponential",
                                   endpoints = "shared")))
  expect_equal(sum(out$coupling / 300), out$variance[["large"]],
               tolerance = 1e-8)
})

test_that("dropout inflates every variance", {
  complete <- do.call(design_components, reference)
  reduced <- do.call(design_components, c(reference,
                                          list(dropout = 0.4,
                                               dropout_shape = "exponential")))
  expect_true(all(reduced$variance > complete$variance))
})

test_that("the requirement on the two correlations is met across a grid", {
  grid <- expand.grid(dispersion = c(0.3, 0.8, 2), dropout = c(0, 0.2, 0.5),
                      n3 = c(90, 300, 800), endpoint = c("phase2", "phase3"),
                      stringsAsFactors = FALSE)
  for (i in seq_len(nrow(grid))) {
    shape <- if (grid$dropout[i] == 0) "complete" else "exponential"
    out <- design_components(rate1 = 2.6, rate2 = 1.3,
                             dispersion1 = grid$dispersion[i],
                             dispersion2 = grid$dispersion[i],
                             shared = grid$dispersion[i], n1 = 60, n2 = 90,
                             n3 = grid$n3[i], window1 = 0.5, window2 = 1,
                             dropout = grid$dropout[i], dropout_shape = shape,
                             small_endpoint = grid$endpoint[i])
    expect_gte(out$correlation[["interim_small"]],
               out$correlation[["interim_large"]])
  }
})

test_that("a Phase 3 small trial lowers the correlation with the interim", {
  first <- do.call(design_components, reference)
  second <- do.call(design_components,
                    c(reference, list(small_endpoint = "phase3")))
  expect_lt(second$correlation[["interim_small"]],
            first$correlation[["interim_small"]])
})

test_that("the means are the log rate ratios", {
  out <- do.call(design_components, c(reference,
                                      list(ratio1 = 0.65, ratio2 = 0.75)))
  expect_equal(out$mean[["interim"]], log(0.65))
  expect_equal(out$mean[["small"]], log(0.65))
  expect_equal(out$mean[["large"]], log(0.75))
})

test_that("invalid arguments are rejected", {
  expect_error(do.call(design_components, c(reference[-5], list(shared = 2))),
               "smaller of the two dispersions")
  expect_error(do.call(design_components,
                       c(reference, list(dropout = c(0.1, 0.2, 0.3)))),
               "one or two values")
})
