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
