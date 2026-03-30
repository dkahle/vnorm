fake_cmdstan_fit <- function(draws_df) {
  list(
    draws = function(format = "df", inc_warmup = FALSE) {
      draws_df
    }
  )
}

test_that("variety_solve validates input class and dimensionality", {
  expect_error(
    variety_solve(123),
    "`polylist` must be an mpoly or mpolyList object."
  )

  expect_error(
    variety_solve(mp(c("x", "y", "x + y")), n = 10),
    "not zero-dimensional"
  )
})

test_that("variety_solve returns rounded posterior means using rvnorm samples", {
  p <- mp("x")
  fake_df <- data.frame(x = c(0.11, 0.09, 0.10), lp__ = c(0, 0, 0))

  testthat::local_mocked_bindings(
    rvnorm = function(...) fake_cmdstan_fit(fake_df),
    .package = "vnorm"
  )

  out <- variety_solve(p, n = 3, sig_digit = 2)

  expect_type(out, "double")
  expect_named(out, "x")
  expect_equal(unname(out), 0.1)
})

test_that("variety_solve can return both stanfit and results", {
  p <- mp("x")
  fake_df <- data.frame(x = c(1.01, 0.99), lp__ = c(0, 0))
  fake_fit <- fake_cmdstan_fit(fake_df)

  testthat::local_mocked_bindings(
    rvnorm = function(...) fake_fit,
    .package = "vnorm"
  )

  out <- variety_solve(p, n = 2, sig_digit = 2, stanfit = TRUE)

  expect_true(is.list(out))
  expect_true(all(c("stanfit", "results") %in% names(out)))
  expect_identical(out$stanfit, fake_fit)
  expect_equal(unname(out$results), 1.0)
})

test_that("variety_solve mpolyList returns correct posterior means", {
  p <- mp(c("x - 1", "y + 2"))
  fake_df <- data.frame(
    x = c(1.001, 0.999, 1.000),
    y = c(-2.001, -1.999, -2.000),
    lp__ = c(0, 0, 0)
  )

  testthat::local_mocked_bindings(
    rvnorm = function(...) fake_cmdstan_fit(fake_df),
    .package = "vnorm"
  )

  out <- variety_solve(p, n = 3, sig_digit = 2)
  expect_type(out, "double")
  expect_named(out, c("x", "y"))
  expect_equal(unname(out), c(1, -2))
})

test_that("variety_solve stanfit=TRUE with mpolyList returns list", {
  p <- mp(c("x - 1", "y + 2"))
  fake_df <- data.frame(
    x = c(1.01, 0.99), y = c(-2.01, -1.99), lp__ = c(0, 0)
  )
  fake_fit <- fake_cmdstan_fit(fake_df)

  testthat::local_mocked_bindings(
    rvnorm = function(...) fake_fit,
    .package = "vnorm"
  )

  out <- variety_solve(p, n = 2, sig_digit = 2, stanfit = TRUE)
  expect_true(is.list(out))
  expect_true(all(c("stanfit", "results") %in% names(out)))
  expect_identical(out$stanfit, fake_fit)
  expect_named(out$results, c("x", "y"))
})

test_that("variety_solve sig_digit parameter controls rounding", {
  p <- mp("x")
  fake_df <- data.frame(x = c(0.12345, 0.12355), lp__ = c(0, 0))

  testthat::local_mocked_bindings(
    rvnorm = function(...) fake_cmdstan_fit(fake_df),
    .package = "vnorm"
  )

  out_2 <- variety_solve(p, n = 2, sig_digit = 2)
  out_4 <- variety_solve(p, n = 2, sig_digit = 4)
  expect_equal(unname(out_2), 0.12)
  expect_equal(unname(out_4), 0.1235)
})
