test_that("project_onto_variety reduces residual for a simple curve", {
  p <- mp("x^2 + y^2 - 1")
  pf <- as.function(p, varorder = c("x", "y"), silent = TRUE)
  x0 <- c(1.4, 0.3)

  x_proj <- project_onto_variety(
    x0, p, al = c(1, 1), adaptive = TRUE, dt = 0.02, tol = 1e-4
  )

  expect_type(x_proj, "double")
  expect_length(x_proj, 2)
  expect_true(all(is.finite(x_proj)))
  expect_lt(abs(pf(x_proj)), abs(pf(x0)))
})

test_that("project_onto_variety supports matrix and data.frame batch input", {
  p <- mp("x^2 + y^2 - 1")
  pts <- rbind(c(1.2, 0.1), c(-1.3, 0.4), c(0.2, -1.4))
  colnames(pts) <- c("x", "y")

  mat_out <- project_onto_variety(
    pts, p, al = c(1, 1), adaptive = FALSE, dt = 0.02, tol = 1e-4
  )
  df_out <- project_onto_variety(
    as.data.frame(pts), p, al = c(1, 1), adaptive = FALSE, dt = 0.02, tol = 1e-4
  )

  expect_true(is.matrix(mat_out))
  expect_equal(dim(mat_out), dim(pts))
  expect_true(all(is.finite(mat_out)))

  expect_true(is.data.frame(df_out))
  expect_equal(dim(df_out), dim(pts))
  expect_true(all(is.finite(as.matrix(df_out))))
})

test_that("project_onto_variety warns and rebuilds helpers when partial helpers are supplied", {
  p <- mp("x^2 + y^2 - 1")
  gfunc <- as.function(p, varorder = c("x", "y"), silent = TRUE)

  expect_warning(
    project_onto_variety(c(1.1, 0.2), p, al = c(1, 1), gfunc = gfunc, tol = 1e-3),
    "computed internally"
  )
})

test_that("alternative projection methods return finite points", {
  p <- mp("x^2 + y^2 - 1")
  x0 <- c(1.3, 0.2)

  x_lag <- project_onto_variety_lagrange(x0, p, tol = 1e-4)
  x_gd <- suppressWarnings(project_onto_variety_gradient_descent(
    x0, p, method = "fixed", ga = 0.01, maxit = 200, tol = 1e-4
  ))
  x_newt <- suppressWarnings(project_onto_variety_newton(
    x0, p, method = "fixed", ga = 0.1, maxit = 50, tol = 1e-4
  ))

  expect_true(all(is.finite(x_lag)))
  expect_true(all(is.finite(x_gd)))
  expect_true(all(is.finite(x_newt)))
  expect_equal(length(x_lag), 2)
  expect_equal(length(x_gd), 2)
  expect_equal(length(x_newt), 2)
})

test_that("project_onto_variety validates required inputs and patch vector", {
  p <- mp("x^2 + y^2 - 1")

  expect_error(project_onto_variety(poly = p), "`x0` must be supplied")
  expect_error(project_onto_variety(c(1, 0)), "`poly` must be supplied")
  expect_error(
    project_onto_variety(c(1, 0), p, al = c(NA_real_, 1)),
    "`al` must be a numeric vector"
  )
  expect_error(
    project_onto_variety(c(1, 0), p, dt = 0),
    "`dt` must be a positive"
  )
  expect_error(
    project_onto_variety(c(1, 0), p, n_correct = 1.5),
    "`n_correct` must be a non-negative integer"
  )
  expect_error(
    project_onto_variety(data.frame(x = "a", y = "b"), p),
    "`x0` must be finite numeric"
  )
  expect_error(
    project_onto_variety(c(1, 0, 2), p),
    "one entry per variable"
  )
})

test_that("project_onto_variety errors instead of hanging at dt_min", {
  p <- mp("x^2 + y^2 - 1")
  expect_error(
    project_onto_variety(
      c(1, 1), p,
      dt = 0.01,
      dt_min = 0.01,
      error_tol = .Machine$double.eps,
      al = c(1, 1)
    ),
    "failed to meet `error_tol`"
  )
})

test_that("project_onto_variety supports fixed-step path with messages", {
  p <- mp("x^2 + y^2 - 1")
  expect_message(
    project_onto_variety(
      c(1.2, 0.1), p, al = c(1, 1), adaptive = FALSE, dt = 0.05, message = TRUE
    )
  )
})

test_that("project_onto_variety fixed-step path reaches t_end", {
  p <- mp("x^2 + y^2 - 1")
  pf <- as.function(p, varorder = c("x", "y"), silent = TRUE)

  expect_warning(
    x_proj <- project_onto_variety(
      c(1, 1), p, adaptive = FALSE, dt = 0.03, al = c(1, 1), tol = 1e-6
    ),
    NA
  )
  expect_lt(abs(pf(x_proj)), 1e-6)
})

test_that("lagrange method supports optim-based branch and tibble batch input", {
  p <- mp("x^2 + y^2 - 1")
  x0 <- c(1.15, 0.2)

  x_lag_bfgs <- suppressWarnings(
    project_onto_variety_lagrange(x0, p, method = "BFGS", maxit = 50)
  )
  expect_true(all(is.finite(x_lag_bfgs)))

  pts_tbl <- tibble::tibble(x = c(1.1, -1.2), y = c(0.2, 0.3))
  out_tbl <- project_onto_variety_lagrange(pts_tbl, p, method = "newton")
  expect_true(tibble::is_tibble(out_tbl))
  expect_equal(dim(out_tbl), dim(pts_tbl))
})

test_that("lagrange optim branch minimizes sum of squared stationarity residuals", {
  p <- mp("x^2 + y^2 - 1")
  pf <- as.function(p, varorder = c("x", "y"), silent = TRUE)

  expect_warning(
    x_lag_bfgs <- project_onto_variety_lagrange(
      c(-3, -3), p, method = "BFGS", maxit = 100, tol = 1e-4
    ),
    NA
  )

  expect_lt(abs(pf(x_lag_bfgs)), 1e-4)
})

test_that("gradient descent and newton methods cover line-search branches", {
  p <- mp("x^2 + y^2 - 1")
  x0 <- c(1.2, 0.25)

  x_gd_line <- suppressWarnings(project_onto_variety_gradient_descent(
    x0, p, method = "line", maxit = 150, tol = 1e-4
  ))
  x_gd_opt <- suppressWarnings(project_onto_variety_gradient_descent(
    x0, p, method = "optimal", maxit = 150, tol = 1e-4
  ))
  x_newt_line <- suppressWarnings(project_onto_variety_newton(
    x0, p, method = "line", maxit = 50, tol = 1e-4
  ))

  expect_true(all(is.finite(x_gd_line)))
  expect_true(all(is.finite(x_gd_opt)))
  expect_true(all(is.finite(x_newt_line)))
})

test_that("gradient descent and newton support data frame batch inputs", {
  p <- mp("x^2 + y^2 - 1")
  pts <- data.frame(x = c(1.1, -1.2), y = c(0.15, 0.4))

  gd_out <- suppressWarnings(project_onto_variety_gradient_descent(
    pts, p, method = "fixed", maxit = 100, tol = 1e-4
  ))
  newt_out <- suppressWarnings(project_onto_variety_newton(
    pts, p, method = "fixed", maxit = 50, tol = 1e-4
  ))

  expect_true(is.data.frame(gd_out))
  expect_true(is.data.frame(newt_out))
  expect_equal(dim(gd_out), dim(pts))
  expect_equal(dim(newt_out), dim(pts))
})

test_that("project_onto_variety n_correct parameter affects convergence", {
  p <- mp("x^2 + y^2 - 1")
  pf <- as.function(p, varorder = c("x", "y"), silent = TRUE)
  x_few <- suppressWarnings(
    project_onto_variety(c(1.4, 0.3), p, al = c(1, 1), n_correct = 0, tol = 1)
  )
  x_many <- project_onto_variety(c(1.4, 0.3), p, al = c(1, 1), n_correct = 5, tol = 1e-6)
  expect_true(abs(pf(x_many)) <= abs(pf(x_few)) + 1e-6)
})

test_that("project_onto_variety works on non-circle polynomial (ellipse)", {
  p <- mp("x^2 + 4 y^2 - 1")
  pf <- as.function(p, varorder = c("x", "y"), silent = TRUE)
  x_proj <- project_onto_variety(c(1.5, 0.5), p, al = c(1, 1), tol = 1e-3)
  expect_length(x_proj, 2)
  expect_lt(abs(pf(x_proj)), 1e-2)
})

test_that("project_onto_variety works in 3D", {
  p <- mp("x^2 + y^2 + z^2 - 1")
  pf <- as.function(p, varorder = c("x", "y", "z"), silent = TRUE)
  x_proj <- project_onto_variety(c(1, 1, 1), p, al = c(1, 1), tol = 1e-3)
  expect_length(x_proj, 3)
  expect_lt(abs(pf(x_proj)), 1e-2)
})

test_that("project_onto_variety tibble input returns tibble", {
  p <- mp("x^2 + y^2 - 1")
  pts <- tibble::tibble(x = c(1.2, -1.1), y = c(0.1, 0.3))
  out <- project_onto_variety(pts, p, al = c(1, 1), tol = 1e-3)
  expect_true(tibble::is_tibble(out))
  expect_equal(nrow(out), 2)
  expect_equal(ncol(out), 2)
})

test_that("project_onto_variety emits warning on convergence failure", {
  p <- mp("x^2 + y^2 - 1")
  expect_warning(
    project_onto_variety(
      c(100, 100), p, al = c(1, 1),
      adaptive = FALSE, dt = 0.5, tol = 1e-12
    ),
    "Tolerance not met"
  )
})

test_that("all projection functions use sorted variable order by default", {
  # vars() appearance order of this polynomial is (y, x); every projector
  # must interpret x0 in sorted order (x, y)
  p <- mp("y^2 - x^3 - x")
  pf <- as.function(p, varorder = c("x", "y"), silent = TRUE)
  x0 <- c(0.5, 0.9)

  # defaults must agree with an explicit sorted varorder
  for (f in list(
    project_onto_variety,
    project_onto_variety_lagrange,
    project_onto_variety_newton,
    project_onto_variety_gradient_descent
  )) {
    out <- suppressWarnings(f(x0, p))
    out_explicit <- suppressWarnings(f(x0, p, varorder = c("x", "y")))
    expect_equal(out, out_explicit)
  }

  # convergent methods land on the variety in sorted coordinates
  for (f in list(
    project_onto_variety,
    project_onto_variety_lagrange,
    project_onto_variety_newton,
    project_onto_variety_gradient_descent
  )) {
    expect_lt(abs(pf(suppressWarnings(f(x0, p)))), 1e-4)
  }
})

test_that("projection functions reject constant polynomials", {
  expect_error(project_onto_variety(c(1, 1), mp("3")), "at least one variable")
  expect_error(
    project_onto_variety_lagrange(c(1, 1), mp("3")),
    "at least one variable"
  )
  expect_error(
    project_onto_variety_gradient_descent(c(1, 1), mp("3")),
    "at least one variable"
  )
  expect_error(
    project_onto_variety_newton(c(1, 1), mp("3")),
    "at least one variable"
  )
})

test_that("projection functions handle one-variable polynomials", {
  p <- mp("x^2 - 2")
  for (f in list(
    project_onto_variety,
    project_onto_variety_lagrange,
    project_onto_variety_newton,
    project_onto_variety_gradient_descent
  )) {
    expect_equal(suppressWarnings(f(3, p)), sqrt(2), tolerance = 1e-5)
  }
})

test_that("lagrange projection uses the starting point at full precision", {
  # the Lagrangian was built by printing x0 into a formula, which rounds it
  # to 15 significant digits (here by about 2e-5)
  x0 <- c(123456789012.3456789, 123456789014.3456789)
  x_proj <- project_onto_variety_lagrange(x0, mp("x - y"))
  expect_lt(max(abs(x_proj - mean(x0))), 1e-6)
})

test_that("batch projections return one row per point", {
  # apply() simplifies one-coordinate results to a vector, which t() turned
  # into a single row holding every point
  p1 <- mp("x^2 - 2")
  pts1 <- matrix(c(2, 3, 4), ncol = 1, dimnames = list(NULL, "x"))
  p2 <- mp("x^2 + y^2 - 1")
  pts2 <- rbind(c(1.2, 0.1), c(-1.3, 0.4), c(0.2, -1.4))

  for (f in list(
    project_onto_variety,
    project_onto_variety_lagrange,
    project_onto_variety_newton,
    project_onto_variety_gradient_descent
  )) {
    mat_out <- suppressWarnings(f(pts1, p1))
    expect_equal(dim(mat_out), c(3L, 1L))
    expect_equal(colnames(mat_out), "x")
    # each point lands on a root (gradient descent can step past x = sqrt(2))
    expect_equal(abs(as.numeric(mat_out)), rep(sqrt(2), 3), tolerance = 1e-5)

    df_out <- suppressWarnings(f(as.data.frame(pts1), p1))
    expect_s3_class(df_out, "data.frame")
    expect_named(df_out, "x")
    expect_equal(df_out$x, as.numeric(mat_out))

    out2 <- suppressWarnings(f(pts2, p2))
    expect_equal(dim(out2), dim(pts2))
    for (i in seq_len(nrow(pts2))) {
      expect_equal(out2[i, ], suppressWarnings(f(pts2[i, ], p2)))
    }
  }
})

test_that("newton projection lands on varieties with mixed-order hessians", {
  # the hessian rows of g^2 for these have variables in different appearance
  # orders (or are partly constant), which broke the hessian evaluation
  for (ps in c("y^2 - x^3 - x", "x y - 1", "x^2 + y")) {
    p <- mp(ps)
    pf <- as.function(p, varorder = c("x", "y"), silent = TRUE)
    x_newt <- suppressWarnings(project_onto_variety_newton(c(0.7, 1.3), p))
    expect_lt(abs(pf(x_newt)), 1e-4)
  }
})

test_that("gradient descent and newton stop on the residual tolerance", {
  # the loops used to stop once g^2 < tol, i.e. |g| ~ sqrt(tol), and then
  # warn that |g| >= tol on their own converged runs
  p <- mp("x^2 + y^2 - 1")
  pf <- as.function(p, varorder = c("x", "y"), silent = TRUE)

  expect_no_warning(x_gd <- project_onto_variety_gradient_descent(c(0.7, 1.3), p))
  expect_no_warning(x_newt <- project_onto_variety_newton(c(0.7, 1.3), p))
  expect_lt(abs(pf(x_gd)), sqrt(.Machine$double.eps))
  expect_lt(abs(pf(x_newt)), sqrt(.Machine$double.eps))
})

test_that("newton projection handles singular hessians of g^2", {
  # for linear g the hessian of g^2 is 2 grad(g) grad(g)', which is singular;
  # the minimum-norm newton step is then the orthogonal projection
  expect_equal(project_onto_variety_newton(c(2, 3), mp("x + y - 1")), c(0, 1))
  expect_equal(
    project_onto_variety_newton(c(2, 3, -1), mp("x + 2 y - z - 1")),
    c(2, 3, -1) - 8 / 6 * c(1, 2, -1)
  )

  # at the origin the hessian is matrix(c(2, 4, 4, 8), 2) and the gradient
  # c(2, -4) leaves its column space: the least-squares step is c(.12, .24),
  # but pdvnorm()'s row-scaled pseudoinverse gave 0, so projection stopped
  # at the start with residual 1
  p <- mp("1 + x - 2 y + 4 x y")
  expect_no_warning(z <- project_onto_variety_newton(c(0, 0), p))
  expect_lt(abs(as.function(p, silent = TRUE)(z)), 1e-8)
})

test_that("project_onto_variety runs no corrector steps when n_correct = 0", {
  # 1:n_correct iterated over c(1, 0) when n_correct = 0
  p <- mp("x^2 + y^2 - 1")
  count_corrections <- function(n_correct) {
    msgs <- testthat::capture_messages(suppressWarnings(project_onto_variety(
      c(1, 1), p, n_correct = n_correct, adaptive = FALSE, dt = 0.5,
      message = TRUE
    )))
    # predictor messages start with a coordinate, corrector messages with
    # two spaces
    sum(startsWith(msgs, "  "))
  }

  expect_equal(count_corrections(0), 0)
  expect_equal(count_corrections(1), 2)
  expect_equal(count_corrections(3), 6)
})
