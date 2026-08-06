test_that("a single value is used for both groups", {
  expect_equal(expand_group_size(50, "n1"), c(50, 50))
  expect_equal(expand_group_size(1, "n1"), c(1, 1))
})

test_that("a vector of length two is returned unchanged", {
  expect_equal(expand_group_size(c(60, 40), "n1"), c(60, 40))
  expect_equal(expand_group_size(c(1, 999), "n1"), c(1, 999))
})

test_that("a numeric vector is returned", {
  expect_type(expand_group_size(50L, "n1"), "double")
  expect_length(expand_group_size(50, "n1"), 2L)
})

test_that("the argument name appears in the error message", {
  expect_error(expand_group_size(0, "n3"), "'n3'")
})

test_that("invalid arguments are rejected", {
  expect_error(expand_group_size(c(1, 2, 3), "n1"), "one or two")
  expect_error(expand_group_size(-5, "n1"), "one or two")
  expect_error(expand_group_size(0, "n1"), "one or two")
  expect_error(expand_group_size(10.5, "n1"), "one or two")
  expect_error(expand_group_size(NA_real_, "n1"), "one or two")
  expect_error(expand_group_size("50", "n1"), "one or two")
  expect_error(expand_group_size(numeric(0), "n1"), "one or two")
})
