test_that("get_listed_coeficients converts coefficient names to Stan-style names", {
  coefs <- c(1, 2, 3)
  names(coefs) <- c("1", "x", "x^2 y")

  out <- get_listed_coeficients(coefs)

  expect_true(is.list(out))
  expect_equal(names(out), c("b1", "bx", "bx2y"))
  expect_equal(unname(unlist(out)), unname(coefs))
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
  # Should have 3 comma-separated entries
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

  # Different structure -> different hash
  p3 <- mp("x^3 + y^2 - 1")
  name3 <- generate_model_name(p3, homo = TRUE, windowed = FALSE)
  expect_false(name1a == name3)
})
