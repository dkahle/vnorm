test_that("create_stan_code builds single-polynomial Stan code", {
  p <- mp("x^2 + y^2 - 1")

  code_homo <- create_stan_code(p, sd = 0.1, n_eqs = 1L, homo = TRUE)
  code_hetero_w <- create_stan_code(p, sd = 0.1, n_eqs = 1L, w = 5, homo = FALSE)

  expect_type(code_homo, "character")
  expect_length(code_homo, 1)
  expect_match(code_homo, "normal_lpdf\\(")
  expect_match(code_homo, "real g =")
  expect_match(code_homo, "real ndg = sqrt\\(")

  expect_match(code_hetero_w, "real<lower=-5,upper=5> x;")
  expect_match(code_hetero_w, "real<lower=-5,upper=5> y;")
  expect_match(code_hetero_w, "real ndg = 1;")
})

test_that("create_stan_code supports named bounds for single-polynomial case", {
  p <- mp("x^2 + y^2 - 1")
  w <- list(x = c(-2, 1), y = c(-3, 3))

  code <- create_stan_code(p, sd = 0.1, n_eqs = 1L, w = w, homo = TRUE)

  expect_match(code, "real<lower=-2,upper=1> x;")
  expect_match(code, "real<lower=-3,upper=3> y;")
})

test_that("create_stan_code builds multi-polynomial Stan code", {
  p <- mp(c("x^2 + y^2 - 1", "y"))

  code_scalar <- create_stan_code(p, sd = 0.1, n_eqs = 2L, homo = TRUE)
  code_matrix <- create_stan_code(p, sd = diag(c(1, 2)), n_eqs = 2L, homo = FALSE)

  expect_match(code_scalar, "vector\\[2\\] g =")
  expect_match(code_scalar, "matrix\\[2,2\\] J =")
  expect_match(code_scalar, "normal_lpdf\\(")

  expect_match(code_matrix, "cov_matrix\\[2\\] si")
  expect_match(code_matrix, "multi_normal_lpdf\\(")
  expect_match(code_matrix, "matrix\\[2,2\\] J =")
})

test_that("create_stan_code single-variable polynomial", {
  code <- create_stan_code(mp("x^2 - 1"), sd = 0.1, n_eqs = 1L, homo = TRUE)
  expect_match(code, "real x;")
  expect_match(code, "real g =")
  expect_match(code, "normal_lpdf")
  expect_false(grepl("real y;", code))
})

test_that("create_stan_code overdetermined system (n_eqs > n_vars)", {
  p <- mp(c("x", "y", "x + y"))
  code <- create_stan_code(p, sd = 0.1, n_eqs = 3L, homo = TRUE)
  expect_match(code, "vector\\[3\\] g")
  expect_match(code, "matrix\\[3,2\\] J")
  # overdetermined uses (J'*J) \ (J'*g)
  expect_match(code, "J'\\*J", fixed = FALSE)
})

test_that("create_stan_code underdetermined system (n_eqs < n_vars)", {
  p <- mp(c("x + y + z", "x - y"))
  code <- create_stan_code(p, sd = 0.1, n_eqs = 2L, homo = TRUE)
  expect_match(code, "vector\\[2\\] g")
  expect_match(code, "matrix\\[2,3\\] J")
  # underdetermined uses J' * ((J*J') \ g)
  expect_match(code, "J\\*J'", fixed = FALSE)
})

test_that("create_stan_code heteroskedastic multivariate with matrix sd", {
  p <- mp(c("x^2 + y^2 - 1", "y"))
  code <- create_stan_code(p, sd = diag(c(1, 2)), n_eqs = 2L, homo = FALSE)
  expect_match(code, "cov_matrix\\[2\\] si")
  expect_match(code, "multi_normal_lpdf")
})

test_that("create_stan_code w as named list for multivariate case", {
  p <- mp(c("x^2 + y^2 - 1", "y"))
  w <- list(x = c(-2, 2))
  code <- create_stan_code(p, sd = 0.1, n_eqs = 2L, w = w, homo = TRUE)
  expect_match(code, "real<lower=-2,upper=2> x;")
  expect_match(code, "real y;")
})
