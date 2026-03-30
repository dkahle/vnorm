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
    "length\\(`sd`\\) must be 1, `length\\(poly\\)`"
  )
})

test_that("rejection_sampler validates homo=TRUE mpolyList sd vector length", {
  p <- mp(c("x^2 + y^2 - 1", "y"))

  expect_error(
    rejection_sampler(
      n = 5, poly = p, sd = c(0.1, 0.2, 0.3), w = 1.5, homo = TRUE, message = FALSE
    ),
    "length\\(`sd`\\) must be 1, `length\\(vars\\)`"
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

test_that("rejection_sampler validates poly class", {
  expect_error(
    rejection_sampler(
      n = 3, poly = 1:3, vars = c("x", "y"), sd = 0.1, w = 1, message = FALSE
    ),
    "`poly` should be a `mpoly` or `mpolyList`"
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

test_that("rejection_sampler single-variable mpolyList", {
  set.seed(14)
  p <- mp(c("x"))
  samps <- rejection_sampler(n = 6, poly = p, sd = 0.1, w = 1, message = FALSE)
  expect_equal(nrow(samps), 6)
  expect_equal(ncol(samps), 1)
})
