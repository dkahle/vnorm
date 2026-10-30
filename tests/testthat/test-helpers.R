test_that("get_listed_coefficients converts coefficient names to Stan-style names", {
  coefs <- c(1, 2, 3)
  names(coefs) <- c("1", "x", "x^2 y")

  out <- get_listed_coefficients(coefs)

  expect_true(is.list(out))
  expect_equal(names(out), c("b1", "bx", "bx2y"))
  expect_equal(unname(unlist(out)), unname(coefs))
})

test_that("Stan coefficient names canonicalize monomial variable order", {
  p_scrambled <- mp("y x + x^2 - 1")
  p_canonical <- mp("x y + x^2 - 1")

  scrambled_data <- make_coefficients_data(p_scrambled, num_of_vars = 2, deg = 2)
  canonical_data <- make_coefficients_data(p_canonical, num_of_vars = 2, deg = 2)

  expect_equal(scrambled_data, canonical_data)
  expect_equal(scrambled_data[["bxy"]], 1)
  expect_false("byx" %in% names(scrambled_data))
  expect_equal(
    get_coefficients_data(p_scrambled),
    get_coefficients_data(p_canonical)
  )
})

test_that("make_coefficients_data fills missing basis coefficients with zero", {
  p <- mp("x^2 + 2 y")
  out <- make_coefficients_data(p, num_of_vars = 2, deg = 2)

  expect_true(is.list(out))
  expect_true(all(c("bx2", "by", "bxy", "by2", "b1") %in% names(out)))
  expect_equal(out[["bx2"]], 1)
  expect_equal(out[["by"]], 2)
  expect_equal(out[["bxy"]], 0L)
  expect_equal(out[["by2"]], 0L)
})

test_that("get_coefficients_data supports mpoly and mpolyList", {
  p1 <- mp("x^2 + y")
  p2 <- mp(c("x^2 + y", "x - 1"))

  out1 <- get_coefficients_data(p1)
  out2 <- get_coefficients_data(p2)

  expect_true(is.list(out1))
  expect_true(all(c("bx2", "by") %in% names(out1)))

  expect_true(is.list(out2))
  expect_true(any(grepl("_1$", names(out2))))
  expect_true(any(grepl("_2$", names(out2))))
})

test_that("variable remapping helpers round-trip names", {
  p <- mp("a^2 + b - 1")
  remapped <- check_and_replace_vars(p)

  expect_true(inherits(remapped$polynomial, "mpoly"))
  expect_true(setequal(mpoly::vars(remapped$polynomial), c("x", "y")))
  expect_equal(unname(unlist(remapped$mapping[c("x", "y")])), c("a", "b"))

  df <- data.frame(x = 1, y = 2, lp__ = 3, check.names = FALSE)
  renamed <- rename_output_df(df, remapped$mapping)
  expect_true(all(c("a", "b", "lp__") %in% names(renamed)))
})

test_that("Stan string helpers return expected formats", {
  expect_match(mpoly_to_stan(mp("x^2 + y")), "x\\^2")

  ml <- mpolyList(mp("x^2 + y"), mp("x - 1"))
  ml_str <- mpolyList_to_stan(ml)
  expect_type(ml_str, "character")
  expect_match(ml_str, ",")

  d_str <- get_derivative("x", num_of_vars = 2, deg = 2)
  expect_type(d_str, "character")
  expect_true(nchar(d_str) > 0)
})

test_that("get_derivative works for 'y' variable", {
  d_str <- get_derivative("y", num_of_vars = 2, deg = 2)
  expect_type(d_str, "character")
  expect_true(nchar(d_str) > 0)
})

test_that("mpoly_to_stan handles negative coefficients and constants", {
  s1 <- mpoly_to_stan(mp("-3 x^2 + 1"))
  expect_true(grepl("-3", s1))
  s2 <- mpoly_to_stan(mp("5"))
  expect_true(grepl("5", s2))
})

test_that("mpolyList_to_stan handles 3+ polynomials", {
  p1 <- mp("x"); p2 <- mp("y"); p3 <- mp("x + y - 1")
  ml <- structure(list(p1, p2, p3), class = "mpolyList")
  s <- mpolyList_to_stan(ml)
  # should have 3 comma-separated entries
  parts <- strsplit(s, ",")[[1]]
  expect_equal(length(parts), 3)
})

test_that("check_and_replace_vars with 1-variable polynomial", {
  p <- mp("a^2 - 1")
  remapped <- check_and_replace_vars(p)
  expect_equal(sort(mpoly::vars(remapped$polynomial)), "x")
  expect_equal(remapped$mapping[["x"]], "a")
})

test_that("check_and_replace_vars with 3-variable polynomial", {
  p <- mp("a^2 + b^2 + c^2 - 1")
  remapped <- check_and_replace_vars(p)
  expect_true(setequal(mpoly::vars(remapped$polynomial), c("x", "y", "z")))
  expect_equal(length(remapped$mapping), 3)
})

test_that("check_and_replace_vars when variables are already x/y", {
  p <- mp("x^2 + y^2 - 1")
  remapped <- check_and_replace_vars(p)
  expect_equal(length(remapped$mapping), 0)
})

test_that("generate_model_name is deterministic and differentiates structure", {
  p1 <- mp("x^2 + y^2 - 1")
  name1a <- generate_model_name(p1, homo = TRUE, windowed = FALSE)
  name1b <- generate_model_name(p1, homo = TRUE, windowed = FALSE)
  expect_equal(name1a, name1b)

  # homo flag matters
  name1_hetero <- generate_model_name(p1, homo = FALSE, windowed = FALSE)
  expect_false(name1a == name1_hetero)

  # windowed flag matters
  name1_w <- generate_model_name(p1, homo = TRUE, windowed = TRUE)
  expect_false(name1a == name1_w)

  # different structure -> different hash
  p3 <- mp("x^3 + y^2 - 1")
  name3 <- generate_model_name(p3, homo = TRUE, windowed = FALSE)
  expect_false(name1a == name3)
})

test_that("generate_model_name ignores term and in-term variable order", {
  # compile_stan_code() and rvnorm(user_compiled = TRUE) must agree on the
  # cache key for the same polynomial however its terms were written
  expect_identical(
    generate_model_name(mp("x y + x^4 - 1")),
    generate_model_name(mp("y x + x^4 - 1"))
  )
  expect_identical(
    generate_model_name(mp("x^2 y + x y^2 + x^4 - 1")),
    generate_model_name(mp("x y^2 + x^2 y + x^4 - 1"))
  )
  expect_identical(
    get_custom_stan_code(mp("y x + x^4 - 1")),
    get_custom_stan_code(mp("x y + x^4 - 1"))
  )
})

test_that("polynomial-list model keys and code ignore coefficient values", {
  # equations were ordered by their coefficient-bearing terms, so changing a
  # coefficient could reorder them and change the structural cache key
  p <- mp(c("x^2 + y - 1", "x^2 + z - 1"))
  q <- mp(c("2 x^2 + y - 1", "x^2 + z - 1"))

  expect_identical(generate_model_name(p), generate_model_name(q))
  expect_identical(get_custom_stan_code(p), get_custom_stan_code(q))

  # the coefficient data follows the same equation order as the code
  q_data <- get_coefficients_data(q)
  expect_equal(q_data[["bx2_1"]], 2)
  expect_equal(q_data[["by_1"]], 1)
  expect_equal(q_data[["bx2_2"]], 1)
  expect_equal(q_data[["bz_2"]], 1)
  expect_identical(names(q_data), names(get_coefficients_data(p)))
})

test_that("gradient and hessian functions follow varorder", {
  # appearance order of every second-derivative row of this polynomial is
  # (y, x), so building the hessian without varorder swaps coordinates
  p <- mp("y^2 - x^3 - x")
  expect_equal(gradient_function(p, c("x", "y"))(c(0.7, 1.3)), c(-3 * 0.49 - 1, 2.6))
  expect_equal(
    hessian_function(p, c("x", "y"))(c(0.7, 1.3)),
    matrix(c(-4.2, 0, 0, 2), 2, 2)
  )

  # constant and zero entries
  expect_equal(hessian_function(mp("x y"), c("x", "y"))(c(1, 2)), matrix(c(0, 1, 1, 0), 2, 2))
  expect_equal(hessian_function(mp("x^2"), c("x", "y"))(c(3, 4)), matrix(c(2, 0, 0, 0), 2, 2))
})

test_that("gradient and hessian functions handle one-variable polynomials", {
  expect_equal(gradient_function(mp("x - 1"), "x")(5), 1)
  expect_equal(gradient_function(mp("x^3"), "x")(2), 12)
  expect_equal(hessian_function(mp("x^3"), "x")(2), matrix(12, 1, 1))
})

test_that("poly_function evaluates with the exact coefficients", {
  # mpoly's as.function() parses printed coefficients (7 significant digits
  # by default, 15 under the earlier digits = 15 wrapper), and R's parser
  # misreads some 17-digit decimals
  a <- 123456789012345.67
  f <- poly_function(mp("x") - a, "x")
  expect_identical(f(a), 0)
  expect_identical(poly_function(mp("x") - 1 / 3, "x")(1 / 3), 0)
  withr::with_options(list(digits = 3), expect_identical(f(a), 0))
  withr::with_options(list(digits = 17), expect_identical(f(a), 0))

  # mpoly prints coefficients within sqrt(eps) of 1 as 1
  b <- 1 + 1e-10
  expect_identical(poly_function(b * mp("x") - b, "x")(1), 0)
})

test_that("poly_function evaluates points given as vectors or matrix rows", {
  f <- poly_function(mp("x^2 y - 3.5 x + 2"), c("x", "y"))
  expect_identical(f(c(2, 3)), 7)
  expect_identical(f(rbind(c(2, 3), c(1, 1))), c(7, -0.5))

  # varorder orders the coordinates and may name variables poly lacks
  expect_identical(poly_function(mp("x^2 y"), c("y", "x"))(c(3, 2)), 12)
  expect_identical(
    poly_function(mp("x^2 - 1"), c("x", "y"))(cbind(c(0, 2), 5)),
    c(-1, 3)
  )

  # with one variable each vector entry is a point; constants give one
  # value per point
  expect_identical(poly_function(mp("x^2"), "x")(c(1, 2, 3)), c(1, 4, 9))
  expect_identical(
    poly_function(mp("5"), c("x", "y"))(rbind(c(1, 2), c(3, 4))),
    c(5, 5)
  )

  # polynomial lists: one column per polynomial, a vector for a single point
  g <- poly_function(mp(c("x^2 + y", "x y - 1", "5")), c("x", "y"))
  expect_identical(g(c(2, 3)), c(7, 5, 5))
  expect_identical(g(rbind(c(2, 3), c(1, 1))), rbind(c(7, 5, 5), c(2, 0, 5)))
  expect_identical(
    poly_function(mp(c("x^2", "x - 1")), "x")(c(3, 4)),
    rbind(c(9, 2), c(16, 3))
  )

  expect_error(poly_function(mp("x y"), "x"), "missing variable")
  expect_error(f(1:3), "one entry per variable")
  expect_error(f(matrix(1:3, 1)), "one column per variable")
})
