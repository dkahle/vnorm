test_that("make_multivariety_stan_codes returns homoskedastic code with Jacobian", {
  code <- make_multivariety_stan_codes(
    num_of_vars = c(2, 2),
    totaldeg = c(1, 2),
    num_of_poly = 2,
    homo = TRUE,
    windowed = TRUE
  )

  expect_type(code, "character")
  expect_match(code, "data \\{")
  expect_match(code, "real w;")
  expect_match(code, "vector\\[2\\] g")
  expect_match(code, "matrix\\[2,2\\] J")
  expect_match(code, "normal_lpdf\\(0\\.00 \\| J' \\* \\(\\(J\\*J'\\) \\\\ g\\), si\\)")
})

test_that("make_multivariety_stan_codes heteroskedastic path uses identity Jacobian surrogate", {
  code <- make_multivariety_stan_codes(
    num_of_vars = c(2, 1),
    totaldeg = c(2, 1),
    num_of_poly = 2,
    homo = FALSE,
    windowed = FALSE
  )

  expect_false(grepl("real w;", code, fixed = TRUE))
  expect_match(code, "matrix\\[2,2\\] J")
  expect_match(code, "[ 1, 0 ]", fixed = TRUE)
  expect_match(code, "[ 0, 1 ]", fixed = TRUE)
})
