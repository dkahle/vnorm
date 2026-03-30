test_that("pdvnorm matches normal density for p(x)=x", {
  p <- mp("x")
  x <- c(-1, 0, 1)

  dens <- pdvnorm(x, p, sigma = 1)
  log_dens <- pdvnorm(x, p, sigma = 1, log = TRUE)

  expect_equal(dens, dnorm(x))
  expect_equal(log_dens, dnorm(x, log = TRUE))
})

test_that("pdvnorm supports matrix and data.frame input for mpolyList", {
  p <- mp(c("x", "y"))
  x <- rbind(c(0, 0), c(1, 2), c(-1, 1))

  out_homo <- pdvnorm(x, p, sigma = c(1, 2), homo = TRUE)
  out_hetero <- pdvnorm(as.data.frame(x), p, sigma = diag(c(1, 4)), homo = FALSE)

  expect_type(out_homo, "double")
  expect_length(out_homo, nrow(x))
  expect_true(all(is.finite(out_homo)))

  expect_type(out_hetero, "double")
  expect_length(out_hetero, nrow(x))
  expect_true(all(is.finite(out_hetero)))
})

test_that("pdvnorm validates sigma shape and positivity", {
  p <- mp(c("x", "y"))
  x <- c(0, 0)

  expect_error(
    pdvnorm(x, p, sigma = c(1, 2, 3), homo = FALSE),
    "length\\(sigma\\) must be m"
  )
  expect_error(
    pdvnorm(x, p, sigma = c(1, -1), homo = TRUE),
    "vector entries must be positive"
  )
  expect_error(
    pdvnorm(x, p, sigma = matrix(c(1, 2, 2, 4), 2, 2), homo = TRUE),
    "positive definite"
  )
})

test_that("pdvnorm errors for invalid poly class", {
  expect_error(pdvnorm(1, poly = 123, sigma = 1), "'poly' must be either")
})

test_that("pdvnorm errors for invalid sigma in single-polynomial case", {
  expect_error(pdvnorm(1, mp("x"), sigma = -1), "single positive numeric")
  expect_error(pdvnorm(1, mp("x"), sigma = c(1, 2)), "single positive numeric")
  expect_error(pdvnorm(1, mp("x"), sigma = "a"), "single positive numeric")
})

test_that("pdvnorm errors for wrong x dimensions (single poly)", {
  expect_error(pdvnorm(c(1, 2, 3), mp("x^2 + y^2 - 1"), sigma = 1), "length n or a matrix")
  expect_error(pdvnorm(matrix(1:9, ncol = 3), mp("x^2 + y^2 - 1"), sigma = 1), "length n or a matrix")
})

test_that("pdvnorm errors for non-finite x in multivariate case", {
  expect_error(pdvnorm(c(1, NaN), mp(c("x", "y")), sigma = 1), "'x' must be finite")
})

test_that("pdvnorm errors for non-finite sigma in multivariate case", {
  expect_error(pdvnorm(c(0, 0), mp(c("x", "y")), sigma = Inf), "'sigma' must be finite")
})

test_that("pdvnorm errors for non-positive scalar sigma in multivariate case", {
  expect_error(pdvnorm(c(0, 0), mp(c("x", "y")), sigma = -1), "'sigma' must be positive")
})

test_that("pdvnorm validates sigma matrix dimensions for homo and hetero", {
  p <- mp(c("x", "y", "x + y"))
  expect_error(pdvnorm(c(0, 0), p, sigma = diag(3), homo = TRUE), "n x n when homo=TRUE")
  expect_error(pdvnorm(c(0, 0), p, sigma = diag(2), homo = FALSE), "m x m when homo=FALSE")
})

test_that("pdvnorm single-mpoly with matrix and data.frame input", {
  p <- mp("x^2 + y^2 - 1")
  x_mat <- matrix(c(1, 0, 0, 1), nrow = 2)
  d_mat <- pdvnorm(x_mat, p, sigma = 0.1)
  expect_length(d_mat, 2)
  expect_true(all(is.finite(d_mat)))

  d_df <- pdvnorm(data.frame(x = c(1, 0), y = c(0, 1)), p, sigma = 0.1)
  expect_length(d_df, 2)
  expect_true(all(is.finite(d_df)))
})

test_that("pdvnorm returns Inf for zero gradient in homoskedastic case", {
  d <- pdvnorm(0, mp("x^2"), sigma = 1, homo = TRUE)
  expect_equal(d, Inf)
})

test_that("pdvnorm log=TRUE works for multivariate case", {
  p <- mp(c("x", "y"))
  x <- c(1, 2)
  d_log <- pdvnorm(x, p, sigma = 1, log = TRUE)
  d_exp <- pdvnorm(x, p, sigma = 1, log = FALSE)
  expect_equal(exp(d_log), d_exp, tolerance = 1e-10)
})

test_that("pdvnorm handles underdetermined and overdetermined systems", {
  # Underdetermined: 1 poly in 2 vars
  p_under <- mp("x + y")
  d1 <- pdvnorm(c(1, 1), p_under, sigma = 1, homo = TRUE)
  expect_true(is.finite(d1))
  expect_true(d1 > 0)

  # Overdetermined: 3 polys in 2 vars
  p_over <- mp(c("x", "y", "x + y"))
  d2 <- pdvnorm(c(0, 0), p_over, sigma = diag(2), homo = TRUE)
  expect_true(is.finite(d2))
  expect_true(d2 > 0)
})
