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

test_that("pdvnorm length-one mpolyList agrees with single mpoly", {
  p <- mp("x^2 + y^2 - 1")
  p_list <- structure(list(p), class = "mpolyList")
  x <- rbind(c(1, 0), c(0.5, 0.5), c(2, -1))

  expect_equal(
    pdvnorm(x, p_list, sd = 0.1, homo = TRUE),
    pdvnorm(x, p, sd = 0.1, homo = TRUE)
  )
  expect_equal(
    pdvnorm(x, p_list, sd = 0.1, homo = FALSE),
    pdvnorm(x, p, sd = 0.1, homo = FALSE)
  )
  expect_equal(
    pdvnorm(x, p_list, sd = 0.1, homo = TRUE, log = TRUE),
    pdvnorm(x, p, sd = 0.1, homo = TRUE, log = TRUE)
  )
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

test_that("pdvnorm is invariant to term order in the polynomial", {
  p_xy <- mp("-1 x^2 (x + 1) + y^2")
  p_yx <- mp("y^2 - x^2 (x + 1)")
  x <- rbind(c(0.5, 0.8), c(-0.3, 0.2), c(1, 1))

  expect_equal(pdvnorm(x, p_yx, sd = .1), pdvnorm(x, p_xy, sd = .1))
  expect_equal(
    pdvnorm(x, p_yx, sd = .1, homo = FALSE),
    pdvnorm(x, p_xy, sd = .1, homo = FALSE)
  )
  expect_equal(
    pdvnorm(x, p_yx, sd = .1, log = TRUE),
    pdvnorm(x, p_xy, sd = .1, log = TRUE)
  )
  expect_equal(pdvnorm(x, p_yx, Sigma = .1^2), pdvnorm(x, p_xy, Sigma = .1^2))

  # the previously-correct x-first evaluation must not change
  expect_equal(
    pdvnorm(c(0.5, 0.8), p_xy, sd = .1), 2.136453,
    tolerance = 1e-6
  )
})

test_that("constant-gradient fast path respects sorted variable order", {
  # y appears first in the terms, so vars() order is (y, x); the gradient with
  # respect to sorted order (x, y) is (-3, 2)
  p <- mp("2 y - 3 x")
  x <- rbind(c(0.5, 0.8), c(0, 0), c(-1, 2))
  g <- 2 * x[, 2] - 3 * x[, 1]
  expect_equal(pdvnorm(x, p, sd = .1), dnorm(g / sqrt(13), sd = .1))
})

test_that("gradient evaluation handles components that drop variables", {
  # gradient of y^2 + x is (1, 2y): the x-component has no variables, which
  # previously made the gradient closure mis-evaluate points
  p <- mp("y^2 + x")
  x <- rbind(c(0.25, 0.5), c(-1, 2))
  g <- x[, 2]^2 + x[, 1]
  ndg <- sqrt(1 + (2 * x[, 2])^2)
  expect_equal(pdvnorm(x, p, sd = .2), dnorm(g / ndg, sd = .2))
})

test_that("pdvnorm mpolyList path is invariant to term order", {
  ps_a <- mp(c("y x - 1", "x^2 + y^2 - 4"))
  ps_b <- mp(c("x y - 1", "x^2 + y^2 - 4"))
  x <- rbind(c(0.5, 2), c(1, 1))

  expect_equal(pdvnorm(x, ps_a, sd = .1), pdvnorm(x, ps_b, sd = .1))
  expect_equal(
    pdvnorm(x, ps_a, sd = c(.1, .2), homo = FALSE),
    pdvnorm(x, ps_b, sd = c(.1, .2), homo = FALSE)
  )
})

test_that("pdvnorm matches x entries to variables by name", {
  p <- mp("x^2 + y^2 - 1")

  expect_equal(
    pdvnorm(cbind(y = 0.8, x = 0.5), p, sd = .1),
    pdvnorm(cbind(x = 0.5, y = 0.8), p, sd = .1)
  )
  expect_equal(
    pdvnorm(cbind(y = 0.8, x = 0.5), p, sd = .1),
    pdvnorm(c(0.5, 0.8), p, sd = .1)
  )
  expect_equal(
    pdvnorm(c(y = 0.8, x = 0.5), p, sd = .1),
    pdvnorm(c(0.5, 0.8), p, sd = .1)
  )

  # mpolyList path
  ps <- mp(c("x^2 + y^2 - 1", "x y - 0.25"))
  expect_equal(
    pdvnorm(cbind(y = 0.3, x = 0.8), ps, sd = .1),
    pdvnorm(c(0.8, 0.3), ps, sd = .1)
  )

  # names disjoint from the variables (e.g. data frame defaults) are positional
  X <- rbind(c(1, 1), c(0, 3))
  expect_equal(pdvnorm(as.data.frame(X), p, sd = 1), pdvnorm(X, p, sd = 1))
})

test_that("pdvnorm errors on partial or extra variable names", {
  p <- mp("x^2 + y^2 - 1")
  expect_error(pdvnorm(cbind(x = 0.5, z = 0.8), p, sd = .1), "`y`")
  expect_error(pdvnorm(cbind(x = 0.5, y = 0.8, z = 1), p, sd = .1), "`z`")
})

test_that("pdvnorm density has the symmetry of the variety, not its reflection", {
  # y^2 = x^2 (x + 1) is symmetric about the x-axis; the historical swapped
  # evaluation was symmetric about the y-axis instead
  p_yx <- mp("y^2 - x^2 (x + 1)")
  pts <- rbind(c(0.5, 0.8), c(0.5, -0.8), c(-0.5, 0.8))
  d <- pdvnorm(pts, p_yx, sd = .1)
  expect_equal(d[1], d[2])
  expect_false(isTRUE(all.equal(d[1], d[3])))
})

test_that("Sigma-matrix single-polynomial branch is term-order invariant", {
  p_xy <- mp("-1 x^2 (x + 1) + y^2")
  p_yx <- mp("y^2 - x^2 (x + 1)")
  S <- matrix(c(.04, .01, .01, .09), nrow = 2)
  x <- rbind(c(0.5, 0.8), c(-0.3, 0.2))

  expect_equal(pdvnorm(x, p_yx, Sigma = S), pdvnorm(x, p_xy, Sigma = S))
  expect_equal(
    pdvnorm(cbind(y = 0.8, x = 0.5), p_yx, Sigma = S),
    pdvnorm(cbind(x = 0.5, y = 0.8), p_yx, Sigma = S)
  )
})

test_that("pdvnorm vanishes at critical points off the variety", {
  # x = 0 is a critical point of x^2 - 1 with g(0) = -1: density tends to 0
  p <- mp("x^2 - 1")
  expect_equal(pdvnorm(0, p, sd = .1), 0)
  expect_equal(pdvnorm(0, p, sd = .1, log = TRUE), -Inf)

  # the origin is a singular point ON V(x y): density diverges
  p_sing <- mp("x y")
  expect_equal(pdvnorm(c(0, 0), p_sing, sd = .1), Inf)
})

test_that("pdvnorm rejects constant polynomials", {
  expect_error(pdvnorm(1, mp("3"), sd = 1), "at least one variable")
})

test_that("pdvnorm evaluates polynomial coefficients at full precision", {
  # mpoly's as.function() prints coefficients with getOption("digits") (7),
  # which evaluated x - 1.23456789012345 as x - 1.234568
  g0 <- -1.23456789012345
  expect_equal(
    pdvnorm(0, mp("x - 1.23456789012345"), sd = 1, homo = FALSE, log = TRUE),
    stats::dnorm(g0, log = TRUE),
    tolerance = 1e-12
  )
  expect_equal(
    pdvnorm(
      c(0, 0), mp(c("x - 1.23456789012345", "y")),
      sd = 1, homo = FALSE, log = TRUE
    ),
    sum(stats::dnorm(c(g0, 0), log = TRUE)),
    tolerance = 1e-12
  )
})

test_that("pdvnorm mpolyList homo density is 0 off the variety where J loses rank", {
  # J = diag(2x, 2y) is singular on the axes; the pseudoinverse dropped the
  # residual component outside its column space, so the off-variety origin
  # scored the maximal density, the same as the solution (1, 1)
  p <- mp(c("x^2 - 1", "y^2 - 1"))
  expect_equal(pdvnorm(c(0, 0), p, sd = 0.1), 0)
  expect_equal(pdvnorm(c(0, 0.5), p, sd = 0.1, log = TRUE), -Inf)
  expect_gt(pdvnorm(c(1, 1), p, sd = 0.1), 15)

  # constant rank-deficient jacobians: an inconsistent system has an empty
  # variety; a dependent but consistent one keeps the pseudoinverse density
  expect_equal(
    pdvnorm(rbind(c(0, 0), c(1, 2)), mp(c("x + y", "x + y - 1")), sd = 1),
    c(0, 0)
  )
  expect_gt(pdvnorm(c(0.1, 0.2), mp(c("x + y", "2 x + 2 y")), sd = 1), 0)

  # overdetermined systems keep the least-squares convention, even at a
  # rank-deficient point (J = 0 at x = 0)
  expect_true(all(
    pdvnorm(cbind(x = c(0, 1)), mp(c("x^2 - 1", "x^2 - 4")), sd = 0.3) > 0
  ))
})
