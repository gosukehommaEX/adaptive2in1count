test_that("the region has the requested shape and is logical", {
  for (amax in c(0, 1, 10, 40)) {
    region <- expansion_region(amax, 20, 0.6)
    expect_type(region, "logical")
    expect_equal(dim(region), c(amax + 1L, amax + 1L))
  }
})

test_that("the region grows with the cutpoint", {
  sizes <- vapply(c(0.2, 0.5, 1, 2, 5),
                  function(c) sum(expansion_region(60, 20, c)), numeric(1))
  expect_true(all(diff(sizes) > 0))
  for (c in c(0.2, 0.5, 1, 2)) {
    expect_true(all(expansion_region(60, 20, c) <=
                      expansion_region(60, 20, c * 2)))
  }
})

test_that("a cutpoint of one splits on the two interim counts", {
  region <- expansion_region(30, 20, 1)
  expect_true(all(diag(region)))
  expect_true(region[1, 31])
  expect_false(region[31, 1])
})

test_that("the region agrees with the observed rate ratio", {
  amax <- 25
  n1 <- c(30, 20)
  counts <- seq(0, amax) + 0.5
  observed <- outer(counts / n1[1], counts / n1[2], "/")
  region <- expansion_region(amax, n1, 0.9)
  expect_true(all(region[observed < 0.9 * (1 - 1e-9)]))
  expect_false(any(region[observed > 0.9 * (1 + 1e-9)]))
})

test_that("a boundary tie is assigned to the expansion region", {
  region <- expansion_region(4, 20, 0.6)
  expect_true(region[2, 3])
  expect_true(expansion_region(4, c(10, 20), 1.2)[2, 3])
})

test_that("unequal allocation shifts the region", {
  balanced <- sum(expansion_region(60, c(20, 20), 1))
  favoured <- sum(expansion_region(60, c(40, 20), 1))
  expect_gt(favoured, balanced)
})

test_that("invalid arguments are rejected", {
  expect_error(expansion_region(-1, 20, 0.6), "non-negative")
  expect_error(expansion_region(10, 20, 0), "positive")
  expect_error(expansion_region(10, c(1, 2, 3), 0.6), "one or two")
})
