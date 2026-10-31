test_that("rejection_sampler returns bounded samples for a simple polynomial", {
  set.seed(1)
  p <- mp("x^2 + y^2 - 1")

  samps <- rejection_sampler(
    n = 20, poly = p, sd = 0.05, w = 1.5, output = "simple", message = FALSE
  )

  expect_true(is.matrix(samps) || is.data.frame(samps))
  expect_equal(nrow(samps), 20)
  expect_true(all(c("x", "y") %in% colnames(samps)))
  expect_true(all(abs(samps[, "x"]) <= 1.5 + 1e-8))
  expect_true(all(abs(samps[, "y"]) <= 1.5 + 1e-8))
})

test_that("rejection_sampler supports tibble/unif output", {
  set.seed(2)
  p <- mp("x^2 + y^2 - 1")

  samps <- rejection_sampler(
    n = 10,
    poly = p,
    sd = 0.1,
    w = c(-1.25, 1.25),
    dist = "unif",
    output = "tibble",
    message = FALSE
  )

  expect_true(tibble::is_tibble(samps))
  expect_equal(nrow(samps), 10)
  expect_true(all(c("x", "y") %in% names(samps)))
})

test_that("rejection_sampler mpolyList homo=FALSE works with vector sd", {
  set.seed(3)
  p <- mp(c("x^2 + y^2 - 1", "y"))

  expect_no_error({
    samps <- rejection_sampler(
      n = 10,
      poly = p,
      sd = c(0.05, 0.05),
      w = 1.5,
      homo = FALSE,
      message = FALSE
    )
  })
  expect_equal(nrow(samps), 10)
})

test_that("rejection_sampler validates mpolyList vector sd lengths", {
  p <- mp(c("x^2 + y^2 - 1", "y"))

  expect_error(
    rejection_sampler(
      n = 5, poly = p, sd = c(0.1, 0.2, 0.3), w = 1.5, homo = FALSE, message = FALSE
    ),
    "length\\(`sd`\\) must be 1 or `length\\(poly\\)`"
  )
})

test_that("rejection_sampler validates homo=TRUE mpolyList sd vector length", {
  p <- mp(c("x^2 + y^2 - 1", "y"))

  expect_error(
    rejection_sampler(
      n = 5, poly = p, sd = c(0.1, 0.2, 0.3), w = 1.5, homo = TRUE, message = FALSE
    ),
    "length\\(`sd`\\) must be 1 or `length\\(vars\\)`"
  )
})

test_that("rejection_sampler supports mpolyList matrix sd in homo and non-homo modes", {
  set.seed(4)
  p <- mp(c("x^2 + y^2 - 1", "y"))
  s <- diag(c(0.05, 0.08))

  s1 <- rejection_sampler(n = 5, poly = p, sd = s, w = 1.5, homo = TRUE, message = FALSE)
  s2 <- rejection_sampler(n = 5, poly = p, sd = s, w = 1.5, homo = FALSE, message = FALSE)

  expect_equal(nrow(s1), 5)
  expect_equal(nrow(s2), 5)
})

test_that("rejection_sampler handles single-variable polynomial branch", {
  set.seed(5)
  p <- mp("x")
  samps <- rejection_sampler(n = 6, poly = p, sd = 0.1, w = 1, message = FALSE)

  expect_equal(ncol(samps), 1)
  expect_equal(colnames(samps), "x")
  expect_equal(nrow(samps), 6)
})

test_that("rejection_sampler handles gradient norms missing a variable", {
  # |grad(y - x^2)|^2 = 4 x^2 + 1 has no y, which mpoly's as.function()
  # rejected for varorder c("x", "y")
  set.seed(5)
  samps <- rejection_sampler(n = 20, poly = mp("y - x^2"), sd = 0.05)
  expect_equal(dim(samps), c(20L, 2L))
  expect_true(all(abs(samps[, "y"] - samps[, "x"]^2) < 0.5))
})

test_that("rejection_sampler validates poly class", {
  expect_error(
    rejection_sampler(
      n = 3, poly = 1:3, vars = c("x", "y"), sd = 0.1, w = 1, message = FALSE
    ),
    "`poly` should be a `mpoly` or `mpolyList`"
  )
})

test_that("rejection_sampler validates n, output, windows, and scalar unif sd", {
  p <- mp("x^2 + y^2 - 1")

  expect_error(
    rejection_sampler(n = 0, poly = p, sd = 0.1, w = 1, message = FALSE),
    "`n` must be a positive integer"
  )
  expect_error(
    rejection_sampler(n = 2, poly = p, sd = 0.1, w = 1, output = "bad", message = FALSE),
    "'arg' should be one of"
  )
  expect_error(
    rejection_sampler(n = 2, poly = p, sd = 0.1, w = -1, message = FALSE),
    "Scalar `w` must be positive"
  )
  expect_error(
    rejection_sampler(
      n = 2,
      poly = p,
      sd = 0.1,
      w = list(x = c(-1, 1)),
      message = FALSE
    ),
    "missing bounds"
  )
  expect_error(
    rejection_sampler(
      n = 2,
      poly = p,
      sd = c(0.1, 0.2),
      w = 1,
      dist = "unif",
      message = FALSE
    ),
    "requires scalar `sd`"
  )
})

test_that("rejection_sampler errors promptly when proposal box misses variety", {
  set.seed(20)
  expect_error(
    rejection_sampler(
      n = 2,
      poly = mp("x^2 + y^2 - 1"),
      sd = 0.05,
      w = list(x = c(10, 11), y = c(10, 11)),
      max_batches = 2,
      message = FALSE
    ),
    "accepted 0 of 2 requested draws after 2 proposal batches"
  )
})

test_that("rejection_sampler validates matrix sd dimensions and definiteness", {
  p <- mp(c("x^2 + y^2 - 1", "y"))
  expect_error(
    rejection_sampler(n = 2, poly = p, sd = diag(3), w = 1, homo = TRUE, message = FALSE),
    "`length\\(vars\\)` by `length\\(vars\\)`"
  )
  expect_error(
    rejection_sampler(
      n = 2,
      poly = p,
      sd = matrix(c(1, 2, 2, 4), 2, 2),
      w = 1,
      homo = TRUE,
      message = FALSE
    ),
    "positive definite"
  )
})

test_that("rejection_sampler can print progress messages", {
  set.seed(6)
  p <- mp("x^2 + y^2 - 1")

  expect_output(
    rejection_sampler(n = 5, poly = p, sd = 0.1, w = 1.5, message = TRUE),
    "remaining"
  )
})

test_that("rejection_sampler coefficient normalization switches call normalize_coefficients", {
  p <- mp("x^2 + y^2 - 1")
  calls <- 0L

  testthat::local_mocked_bindings(
    normalize_coefficients = function(obj) {
      calls <<- calls + 1L
      obj
    }
  )

  set.seed(7)
  rejection_sampler(
    n = 4, poly = p, sd = 0.1, w = 1.5,
    correct_p_coefficients = TRUE,
    correct_dp_coefficients = TRUE,
    message = FALSE
  )

  expect_gte(calls, 2L)
})

test_that("rejection_sampler dist='norm' works explicitly", {
  set.seed(10)
  samps <- rejection_sampler(
    n = 10, poly = mp("x^2 + y^2 - 1"), sd = 0.05, w = 1.5,
    dist = "norm", message = FALSE
  )
  expect_equal(nrow(samps), 10)
  expect_true(all(c("x", "y") %in% colnames(samps)))
})

test_that("rejection_sampler w as named list with per-variable bounds", {
  set.seed(11)
  w <- list(x = c(-1, 1), y = c(-2, 2))
  samps <- rejection_sampler(
    n = 10, poly = mp("x^2 + y^2 - 1"), sd = 0.1, w = w, message = FALSE
  )
  expect_equal(nrow(samps), 10)
  expect_true(all(samps[, "x"] >= -1 & samps[, "x"] <= 1))
  expect_true(all(samps[, "y"] >= -2 & samps[, "y"] <= 2))
})

test_that("rejection_sampler single-mpoly with homo=FALSE", {
  set.seed(12)
  samps <- rejection_sampler(
    n = 10, poly = mp("x^2 + y^2 - 1"), sd = 0.1, w = 1.5,
    homo = FALSE, message = FALSE
  )
  expect_equal(nrow(samps), 10)
})

test_that("rejection_sampler mpolyList with scalar sd and homo=TRUE", {
  set.seed(13)
  p <- mp(c("x^2 + y^2 - 1", "y"))
  samps <- rejection_sampler(
    n = 10, poly = p, sd = 0.05, w = 1.5, homo = TRUE, message = FALSE
  )
  expect_equal(nrow(samps), 10)
})

test_that("rejection_sampler handles rectangular mpolyList jacobians", {
  set.seed(15)
  p <- mpolyList(mp("x^2 + y^2 - 1"))
  samps <- rejection_sampler(
    n = 6, poly = p, sd = c(0.05, 0.05), w = 1.5,
    homo = TRUE, message = FALSE
  )

  expect_equal(nrow(samps), 6)
  expect_true(all(c("x", "y") %in% colnames(samps)))
})

test_that("rejection_sampler single-variable mpolyList", {
  set.seed(14)
  p <- mp(c("x"))
  samps <- rejection_sampler(n = 6, poly = p, sd = 0.1, w = 1, message = FALSE)
  expect_equal(nrow(samps), 6)
  expect_equal(ncol(samps), 1)
})

test_that("rejection_sampler rejects constant polynomials", {
  expect_error(rejection_sampler(5, mp("3"), sd = .1), "at least one variable")
})

test_that("rejection_sampler scalar sd works for overdetermined systems", {
  # g'(JJ')^{-1}g needs JJ' invertible, which fails whenever there are more
  # equations than variables; all homoskedastic kernels now use J^+ g
  p <- mp(c("x", "y", "x + y"))

  set.seed(1)
  scalar_sd <- rejection_sampler(20, p, sd = 0.3, w = 1)
  set.seed(1)
  vector_sd <- rejection_sampler(20, p, sd = c(0.3, 0.3), w = 1)
  set.seed(1)
  matrix_sd <- rejection_sampler(20, p, sd = diag(c(0.3, 0.3)^2), w = 1)

  expect_equal(nrow(scalar_sd), 20)
  expect_equal(scalar_sd, vector_sd)
  expect_equal(scalar_sd, matrix_sd)
})

test_that("rejection_sampler follows pdvnorm()'s support for dependent systems", {
  # the kernel took J^+ g without pdvnorm()'s rule that the density is 0
  # where a square or underdetermined system's J loses rank and g leaves its
  # column space; for these equations (J has rank 1 everywhere, g = (s, s^2)
  # leaves its column space off s = 0) every draw had density 0
  p <- mp(c("x + y + z", "(x + y + z)^2"))
  expect_error(rejection_sampler(10, p, sd = 0.1, w = 0.5), "density is 0 off the variety")
  expect_error(
    rvnorm(10, p, sd = 0.1, rejection = TRUE, w = 0.5),
    "density is 0 off the variety"
  )

  # linearly dependent equations keep g in J's column space and a positive
  # density, which the kernel now samples (solve(J) used to fail)
  p_lin <- mp(c("x + y", "2 x + 2 y"))
  set.seed(1)
  samps <- rejection_sampler(20, p_lin, sd = 0.1, w = 1)
  expect_equal(nrow(samps), 20)
  expect_true(all(pdvnorm(samps, p_lin, sd = 0.1) > 0))
})

test_that("rejection_sampler evaluates full-rank jacobians with large condition numbers", {
  # the kernel took (J'J)^-1 J' g, and solve() refused J'J at every point
  # (condition number about 1e18); g = J (x, y), so draws are N(0, 0.1^2 I)
  p <- mpolyList(1e-9 * mp("x"), mp("y"), 1e-9 * mp("x") + mp("y"))
  set.seed(1)
  samps <- rejection_sampler(100, p, sd = 0.1, w = 1)
  expect_equal(nrow(samps), 100)
  expect_lt(max(abs(samps)), 0.6)
})

test_that("rejection_sampler samples badly scaled independent equations", {
  # J = diag(100 x^99, 100 y^99) is about diag(1e-20, 1e-60) at points of
  # [-1, 1]^2, where the old check's numerical rank (1) passed for proof of
  # dependent equations
  p <- mp(c("x^100 - 1", "y^100 - 1"))
  set.seed(1)
  samps <- rejection_sampler(10, p, sd = 0.01, w = c(0.99, 1.01))
  expect_equal(nrow(samps), 10)
  expect_true(all(pdvnorm(samps, p, sd = 0.01) > 0))
})

test_that("rejection_sampler decides the support of dependent equations by its kernel", {
  # J has rank 1, and exactly g leaves its column space (0.1 * 0.9 is not
  # 0.3 * 0.3 in doubles), but only by rounding: pdvnorm()'s tolerance keeps
  # the density positive, and so does the kernel
  p <- mp(c("0.1 x + 0.2 y - 0.3", "0.3 x + 0.6 y - 0.9"))
  expect_true(dependent_equations(p, c("x", "y")))
  set.seed(1)
  samps <- rejection_sampler(20, p, sd = 0.1, w = 2)
  expect_equal(nrow(samps), 20)
  expect_true(all(pdvnorm(samps, p, sd = 0.1) > 0))
})

test_that("dependent_equations decides dependence exactly", {
  dep <- function(p) dependent_equations(p, sort(mpoly::vars(p)))
  # J has rank 1 everywhere
  expect_true(dep(mp(c("x + y + z", "(x + y + z)^2"))))
  expect_true(dep(mp(c("x + y", "2 x + 2 y"))))
  expect_true(dep(mp(c("x + y", "2 x + 2 y", "x + y + 1"))))
  # singular only on proper subsets (the axes here), or overdetermined
  expect_false(dep(mp(c("x^2 - 1", "y^2 - 1"))))
  expect_false(dep(mp(c("x", "y", "x + y"))))
  # badly scaled at the fixed points: diag(1e-20, 1e-60), and a J whose J'J
  # solve() refused there
  expect_false(dep(mp(c("x^100 - 1", "y^100 - 1"))))
  expect_false(dep(mp(c(
    "0.000000001 x^2 - 1", "y^2 - 1", "0.000000001 x^2 + y^2 - 2"
  ))))
  # numerically rank 1 at every point (condition number 2e12); in integers
  # 1e-12 against 1 passes 2^53, so the exact stage is inconclusive
  expect_false(dep(mpolyList(mp("x") + 1e-12 * mp("y"), mp("x") + 2e-12 * mp("y"))))
})

test_that("dependent_equations gives the same verdict for rescaled equations", {
  # J = 11' + 1e-8 I has full rank, but elimination in doubles underflowed
  # products of 1e-20-sized coefficients to 0 and took that for dependence
  vars <- paste0("x", 1:5)
  xs <- lapply(vars, mp)
  s <- Reduce(`+`, xs)
  p <- structure(lapply(xs, function(x) s + 1e-8 * x), class = "mpolyList")
  rescale_first <- function(p, scale) {
    structure(c(list(scale * p[[1]]), unclass(p)[-1]), class = "mpolyList")
  }
  for (scale in c(1, 1e-20, 1e20)) {
    expect_false(dependent_equations(scale * p, vars))
    expect_false(dependent_equations(rescale_first(p, scale), vars))
  }

  # dependent equations stay dependent when rescaling keeps their
  # coefficients' exact ratios, as for equal coefficients or doubling
  for (q in list(mp(c("x + y + z", "(x + y + z)^2")), mp(c("x + y", "2 x + 2 y")))) {
    v <- sort(mpoly::vars(q))
    for (scale in c(3, 1e-20, 2^-60)) {
      expect_true(dependent_equations(scale * q, v))
      expect_true(dependent_equations(rescale_first(q, scale), v))
    }
  }
})

test_that("integer_jacobian_rows clears each equation's binary fractions exactly", {
  rows <- function(p) {
    r <- integer_jacobian_rows(p, sort(mpoly::vars(p)))
    if (is.null(r)) return(NULL)
    lapply(r, function(row) vapply(row, function(q) print(q, silent = TRUE), character(1)))
  }
  # 0.2 is exactly twice 0.1, and 0.6 twice 0.3
  expect_identical(
    rows(mp(c("0.1 x + 0.2 y - 0.3", "0.3 x + 0.6 y - 0.9"))),
    list(c("1", "2"), c("1", "2"))
  )
  # 0.3 is not exactly three times 0.1, and as integers the pair passes 2^53
  expect_null(rows(mp(c("0.1 x + 0.3 y", "x"))))

  x <- c(0.1, 3, -0.5, 2^-1074, 1e300, -6)
  parts <- binary_parts(x)
  expect_identical(parts$m * 2^parts$e, x)
  expect_true(all(parts$m %% 2 == 1))
})

test_that("polynomial_matrix_rank finds generic ranks within its budget", {
  rows <- function(p) integer_jacobian_rows(p, c("x", "y", "z"))
  p <- mp(c("x^2 + y^2 + z^2 - 1", "x y - z", "x + y^2 z - 2"))
  expect_identical(polynomial_matrix_rank(rows(p)), 3L)
  expect_identical(
    polynomial_matrix_rank(rows(mp(c("x + y + z", "(x + y + z)^2")))), 1L
  )
  expect_identical(polynomial_matrix_rank(rows(p), max_work = 10), NA_integer_)

  # doubles stop holding every integer at 2^53, where (2^30 + 1)^2 lies
  big <- 2^30 + 1
  expect_identical(
    polynomial_matrix_rank(list(list(big * mp("x"), mp("1")), list(mp("1"), big * mp("y")))),
    NA_integer_
  )
})

test_that("rejection_sampler handles one-variable systems and 1 x 1 scales", {
  # deriv() returns a bare constant mpoly for x - 1, which as.function()
  # rejects; diag() of a scalar builds an identity matrix of that size
  p <- mp(c("x^2 - 1", "x - 1"))

  set.seed(1)
  expect_equal(nrow(rejection_sampler(10, p, sd = 0.3, w = 2)), 10)
  set.seed(1)
  expect_equal(nrow(rejection_sampler(5, p, sd = matrix(0.01), w = 2)), 5)
})

test_that("rejection_sampler grows proposal batches at low acceptance rates", {
  # batches used to shrink to the number of draws still needed, so a 0.1%
  # acceptance rate needed thousands of batches (about log(n) / rate)
  set.seed(1)
  samps <- rejection_sampler(100, mp("x"), sd = 1e-3, w = 1, max_batches = 50)

  expect_equal(nrow(samps), 100)
  expect_lt(max(abs(samps[, "x"])), 0.01)
})
