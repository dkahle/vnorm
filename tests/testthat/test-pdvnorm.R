test_that("pdvnorm matches normal density for p(x)=x", {
  p <- mp("x")
  x <- c(-1, 0, 1)

  dens <- pdvnorm(x, p, sd = 1)
  log_dens <- pdvnorm(x, p, sd = 1, log = TRUE)

  expect_equal(dens, dnorm(x))
  expect_equal(log_dens, dnorm(x, log = TRUE))
})

test_that("pdvnorm uses scalar Sigma as a variance for single polynomials", {
  p <- mp("x")
  x <- c(-1, 0, 1)

  expect_equal(pdvnorm(x, p, Sigma = 4), dnorm(x, sd = 2))
})

test_that("pdvnorm supports Sigma for single multivariable polynomials", {
  p <- mp("x + y")
  x <- c(1, 1)
  Sigma <- diag(c(1, 4))

  normal_distance <- sum(x) / sqrt(2)
  normal_sd <- sqrt(as.numeric(t(c(1, 1) / sqrt(2)) %*% Sigma %*% (c(1, 1) / sqrt(2))))
  expect_equal(
    pdvnorm(x, p, Sigma = Sigma),
    dnorm(normal_distance, sd = normal_sd)
  )
  expect_equal(
    pdvnorm(x, p, Sigma = diag(c(4, 4))),
    pdvnorm(x, p, sd = 2)
  )
  expect_error(
    pdvnorm(x, p, Sigma = Sigma, homo = FALSE),
    "homo = TRUE"
  )
})

test_that("pdvnorm supports matrix and data.frame input for mpolyList", {
  p <- mp(c("x", "y"))
  x <- rbind(c(0, 0), c(1, 2), c(-1, 1))

  out_homo <- pdvnorm(x, p, sd = c(1, 2), homo = TRUE)
  out_hetero <- pdvnorm(as.data.frame(x), p, Sigma = diag(c(1, 4)), homo = FALSE)

  expect_type(out_homo, "double")
  expect_length(out_homo, nrow(x))
  expect_true(all(is.finite(out_homo)))

  expect_type(out_hetero, "double")
  expect_length(out_hetero, nrow(x))
  expect_true(all(is.finite(out_hetero)))
})

test_that("pdvnorm sd vector agrees with diagonal Sigma covariance", {
  p <- mp(c("x", "y"))
  x <- rbind(c(0, 0), c(1, 2), c(-1, 1))

  expect_equal(
    pdvnorm(x, p, sd = c(1, 2)),
    pdvnorm(x, p, Sigma = diag(c(1, 4)))
  )
})

test_that("pdvnorm handles full non-diagonal Sigma", {
  p <- mp(c("x", "y"))
  x <- c(0.5, -0.25)
  Sigma <- matrix(c(1, 0.3, 0.3, 2), 2, 2)

  log_det <- as.numeric(determinant(Sigma, logarithm = TRUE)$modulus)
  quad <- as.numeric(t(x) %*% solve(Sigma, x))
  expected_log <- -0.5 * (2 * log(2 * pi) + log_det + quad)

  expect_equal(pdvnorm(x, p, Sigma = Sigma, log = TRUE), expected_log)
})

test_that("pdvnorm heteroskedastic mpolyList fast path matches manual density", {
  p <- mp(c("x^2 + y", "x - y^2"))
  x <- rbind(c(0.5, -0.25), c(1, 2), c(-1, 0.75))
  Sigma <- matrix(c(1, 0.2, 0.2, 2), 2, 2)
  g_vals <- cbind(x[, 1]^2 + x[, 2], x[, 1] - x[, 2]^2)
  log_det <- as.numeric(determinant(Sigma, logarithm = TRUE)$modulus)
  expected_log <- apply(g_vals, 1, function(g) {
    quad <- as.numeric(t(g) %*% solve(Sigma, g))
    -0.5 * (2 * log(2 * pi) + log_det + quad)
  })

  expect_equal(
    pdvnorm(x, p, Sigma = Sigma, homo = FALSE, log = TRUE),
    expected_log
  )
  expect_equal(
    pdvnorm(x, p, Sigma = Sigma, homo = FALSE),
    exp(expected_log)
  )
})

test_that("pdvnorm constant square homoskedastic fast path matches manual density", {
  p <- mp(c("x", "y"))
  x <- rbind(c(0.5, -0.25), c(1, 2), c(-1, 0.75))
  Sigma <- matrix(c(2, 0.3, 0.3, 1), 2, 2)
  log_det <- as.numeric(determinant(Sigma, logarithm = TRUE)$modulus)
  expected_log <- apply(x, 1, function(v) {
    quad <- as.numeric(t(v) %*% solve(Sigma, v))
    -0.5 * (2 * log(2 * pi) + log_det + quad)
  })

  expect_equal(pdvnorm(x, p, Sigma = Sigma, log = TRUE), expected_log)
  expect_equal(pdvnorm(x, p, Sigma = Sigma), exp(expected_log))
})

test_that("pdvnorm constant overdetermined fast path matches manual density", {
  p <- mp(c("x", "y", "x + y"))
  x <- rbind(c(0.5, -0.25), c(1, 2), c(-1, 0.75))
  Sigma <- matrix(c(2, 0.3, 0.3, 1), 2, 2)
  J <- rbind(c(1, 0), c(0, 1), c(1, 1))
  Jp <- solve(t(J) %*% J) %*% t(J)
  g_vals <- cbind(x[, 1], x[, 2], x[, 1] + x[, 2])
  v_vals <- t(Jp %*% t(g_vals))
  log_det <- as.numeric(determinant(Sigma, logarithm = TRUE)$modulus)
  expected_log <- apply(v_vals, 1, function(v) {
    quad <- as.numeric(t(v) %*% solve(Sigma, v))
    -0.5 * (2 * log(2 * pi) + log_det + quad)
  })

  expect_equal(pdvnorm(x, p, Sigma = Sigma, log = TRUE), expected_log)
  expect_equal(pdvnorm(x, p, Sigma = Sigma), exp(expected_log))
})

test_that("pdvnorm validates sd and Sigma shape and positivity", {
  p <- mp(c("x", "y"))
  x <- c(0, 0)

  expect_error(
    pdvnorm(x, p, sd = c(1, 2, 3), homo = FALSE),
    "length\\(`sd`\\) must be m"
  )
  expect_error(
    pdvnorm(x, p, sd = c(1, -1), homo = TRUE),
    "`sd` vector entries must be positive"
  )
  expect_error(
    pdvnorm(x, p, sd = diag(2), homo = TRUE),
    "Use `Sigma`"
  )
  expect_error(
    pdvnorm(x, p, Sigma = matrix(c(1, 2, 2, 4), 2, 2), homo = TRUE),
    "positive definite"
  )
})

test_that("pdvnorm errors for invalid poly class", {
  expect_error(pdvnorm(1, poly = 123, sd = 1), "'poly' must be either")
})

test_that("pdvnorm errors for invalid sd in single-polynomial case", {
  expect_error(pdvnorm(1, mp("x"), sd = -1), "single positive numeric")
  expect_error(pdvnorm(1, mp("x"), sd = c(1, 2)), "single positive numeric")
  expect_error(pdvnorm(1, mp("x"), sd = "a"), "single positive numeric")
  expect_error(
    pdvnorm(1, mp("x"), Sigma = diag(2)),
    "n x n"
  )
})

test_that("pdvnorm errors for wrong x dimensions (single poly)", {
  expect_error(pdvnorm(c(1, 2, 3), mp("x^2 + y^2 - 1"), sd = 1), "length n or a matrix")
  expect_error(pdvnorm(matrix(1:9, ncol = 3), mp("x^2 + y^2 - 1"), sd = 1), "length n or a matrix")
})

test_that("pdvnorm errors for non-numeric or non-finite x in single-polynomial case", {
  expect_error(pdvnorm(Inf, mp("x"), sd = 1), "'x' must be finite")
  expect_error(
    pdvnorm(data.frame(x = c("a", "b")), mp("x"), sd = 1),
    "'x' must be finite"
  )
})

test_that("pdvnorm errors for non-finite x in multivariate case", {
  expect_error(pdvnorm(c(1, NaN), mp(c("x", "y")), sd = 1), "'x' must be finite")
})

test_that("pdvnorm errors for non-finite Sigma in multivariate case", {
  expect_error(pdvnorm(c(0, 0), mp(c("x", "y")), Sigma = Inf), "`Sigma` must be finite")
})

test_that("pdvnorm errors for non-positive scalar sd in multivariate case", {
  expect_error(pdvnorm(c(0, 0), mp(c("x", "y")), sd = -1), "`sd` must be positive")
})

test_that("pdvnorm validates Sigma matrix dimensions for homo and hetero", {
  p <- mp(c("x", "y", "x + y"))
  expect_error(pdvnorm(c(0, 0), p, Sigma = diag(3), homo = TRUE), "n x n when homo=TRUE")
  expect_error(pdvnorm(c(0, 0), p, Sigma = diag(2), homo = FALSE), "m x m when homo=FALSE")
})

test_that("pdvnorm single-mpoly with matrix and data.frame input", {
  p <- mp("x^2 + y^2 - 1")
  x_mat <- matrix(c(1, 0, 0, 1), nrow = 2)
  d_mat <- pdvnorm(x_mat, p, sd = 0.1)
  expect_length(d_mat, 2)
  expect_true(all(is.finite(d_mat)))

  d_df <- pdvnorm(data.frame(x = c(1, 0), y = c(0, 1)), p, sd = 0.1)
  expect_length(d_df, 2)
  expect_true(all(is.finite(d_df)))
})

test_that("pdvnorm returns Inf for zero gradient in homoskedastic case", {
  d <- pdvnorm(0, mp("x^2"), sd = 1, homo = TRUE)
  expect_equal(d, Inf)
})

test_that("pdvnorm log=TRUE works for multivariate case", {
  p <- mp(c("x", "y"))
  x <- c(1, 2)
  d_log <- pdvnorm(x, p, sd = 1, log = TRUE)
  d_exp <- pdvnorm(x, p, sd = 1, log = FALSE)
  expect_equal(exp(d_log), d_exp, tolerance = 1e-10)
})

test_that("pdvnorm handles underdetermined and overdetermined systems", {
  # underdetermined: 1 poly in 2 vars
  p_under <- mp("x + y")
  d1 <- pdvnorm(c(1, 1), p_under, sd = 1, homo = TRUE)
  expect_true(is.finite(d1))
  expect_true(d1 > 0)

  # overdetermined: 3 polys in 2 vars
  p_over <- mp(c("x", "y", "x + y"))
  d2 <- pdvnorm(c(0, 0), p_over, Sigma = diag(2), homo = TRUE)
  expect_true(is.finite(d2))
  expect_true(d2 > 0)
})

test_that("pdvnorm accepts legacy sigma argument", {
  p <- mp(c("x", "y"))
  x <- c(0.5, 1)

  expect_equal(
    pdvnorm(x, p, sigma = diag(c(1, 4))),
    pdvnorm(x, p, Sigma = diag(c(1, 4)))
  )
  expect_equal(pdvnorm(0, mp("x"), sigma = 2), dnorm(0, sd = 2))
  expect_error(
    pdvnorm(0, mp("x"), sd = 1, sigma = 2),
    "Specify only one"
  )
})
