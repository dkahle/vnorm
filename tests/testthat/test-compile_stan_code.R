skip_if_no_cmdstan <- function() {
  skip_if_not_installed("cmdstanr")
  ver <- tryCatch(cmdstanr::cmdstan_version(), error = function(e) NULL)
  if (is.null(ver)) {
    skip("Skipping compile_stan_code tests: CmdStan path has not been set.")
  }
}

test_that("Compiles model correctly for an mpoly object", {
  skip_if_no_cmdstan()
  p <- mp("x^2 + y^2 - 1")
  result <- compile_stan_code(
    p,
    custom_stan_code = TRUE,
    windowed = FALSE,
    homo = TRUE
  ) |>
    suppressMessages()
  expect_equal(
    result$code(),
    c(
      "data {",
      "  real si;",
      "  real bx2;   real by2;   real b1;",
      "}",
      "parameters {",
      "  real x;",
      "  real y;",
      " }",
      "model {",
      "  real g = bx2*x^2+by2*y^2+b1;",
      "  real dgx = 2*bx2*x;  real dgy = 2*by2*y;",
      "  real ndg = sqrt(dgx^2 + dgy^2);",
      "  target += normal_lpdf(0.00 | g/ndg, si); ",
      "}")
  )
})

test_that("Stops if pre-compiled model already exists", {
  p <- mp("x^2 + y^2")
  expect_error(
    compile_stan_code(p, custom_stan_code = FALSE, windowed = FALSE, homo = TRUE),
    paste0(
      "Pre-compiled model for the general case already exists.",
      "*custom_stan_code = TRUE"
    )
  )
})

test_that("Compiles model correctly for an mpolyList object", {
  skip_if_no_cmdstan()
  p <- mp(c("x^2 + y^2 - 1","y - x"))
  result <- compile_stan_code(
    p,
    custom_stan_code = TRUE,
    windowed = FALSE,
    homo = TRUE
  ) |>
    suppressMessages()
  expect_equal(
    result$code(),
    c(
      "data {",
      "  real si;",
      "  real bx_1;   real by_1;   real bx2_2;   real by2_2;   real b1_2;",
      "}",
      "",
      "parameters {",
      "  real x;",
      "  real y;",
      "}",
      "",
      "transformed parameters {",
      "  vector[2] g = [bx_1*x+by_1*y,bx2_2*x^2+by2_2*y^2+b1_2]';",
      "  matrix[2,2] J = [ ",
      "      [bx_1,by_1],",
      "      [2*bx2_2*x,2*by2_2*y]",
      "    ];",
      "}",
      "",
      "model {",
      "  target += normal_lpdf(0.00 |J \\ g, si);",
      "}"
    )
  )
})

test_that("Creates or updates internal compiled_stan_info cache", {
  skip_if_no_cmdstan()
  poly <- mp("x^2 + y^2 - 1")
  compile_stan_code(poly, custom_stan_code = TRUE, windowed = FALSE, homo = TRUE)
  compiled_info <- vnorm:::get_compiled_stan_info()
  expect_true(is.data.frame(compiled_info))
  expect_gt(nrow(compiled_info), 0)
  expect_true("name" %in% names(compiled_info))
  expect_true("path" %in% names(compiled_info))
})

test_that("Compiles model with box constraints (windowed)", {
  skip_if_no_cmdstan()
  p <- mp("x^2 + y^2 - 1")
  result <- compile_stan_code(p, custom_stan_code = TRUE, windowed = TRUE, homo = TRUE)
  expect_equal(
    result$code(),
    c(
      "data {",
      "  real si;",
      "  real bx2;   real by2;   real b1;  real w;",
      "}",
      "parameters {",
      "  real<lower=-w, upper=w> x;",
      "  real<lower=-w, upper=w> y;",
      " }",
      "model {",
      "  real g = bx2*x^2+by2*y^2+b1;",
      "  real dgx = 2*bx2*x;  real dgy = 2*by2*y;",
      "  real ndg = sqrt(dgx^2 + dgy^2);",
      "  target += normal_lpdf(0.00 | g/ndg, si); ",
      "}"
    )
  )
})

test_that("Compiles model with homoskedastic option", {
  skip_if_no_cmdstan()
  p <- mp("x^2 + y^2 - 1")
  result <- compile_stan_code(
    p,
    custom_stan_code = TRUE,
    windowed = FALSE,
    homo = TRUE
  )
  expect_equal(
    result$code(),
    c(
      "data {",
      "  real si;",
      "  real bx2;   real by2;   real b1;",
      "}",
      "parameters {",
      "  real x;", "  real y;", " }",
      "model {",
      "  real g = bx2*x^2+by2*y^2+b1;",
      "  real dgx = 2*bx2*x;  real dgy = 2*by2*y;",
      "  real ndg = sqrt(dgx^2 + dgy^2);",
      "  target += normal_lpdf(0.00 | g/ndg, si); ",
      "}"
    )
  )
})

test_that("compile_stan_code errors for invalid poly class", {
  expect_error(compile_stan_code(123), "mpoly or mpolyList")
})

test_that("compile_stan_code emits 'cache updated' when adding a new model", {
  skip_if_no_cmdstan()
  on.exit(vnorm:::clear_compiled_stan_info(), add = TRUE)
  vnorm:::clear_compiled_stan_info()

  vnorm:::add_compiled_stan_info("other_model", "/tmp/other.stan")

  expect_message(
    compile_stan_code(mp("x^6 + y^6 - 1"), custom_stan_code = TRUE),
    "cache updated"
  )
})

test_that("compile_stan_code emits 'already contained' on duplicate compile", {
  skip_if_no_cmdstan()
  on.exit(vnorm:::clear_compiled_stan_info(), add = TRUE)
  vnorm:::clear_compiled_stan_info()

  suppressMessages(
    compile_stan_code(mp("x^6 + y^6 - 1"), custom_stan_code = TRUE)
  )

  expect_message(
    compile_stan_code(mp("x^6 + y^6 - 1"), custom_stan_code = TRUE),
    "already contained"
  )
})

test_that("get_custom_stan_code homo=FALSE for mpolyList uses identity Jacobian", {
  code <- vnorm:::get_custom_stan_code(
    mp(c("x^2 + y^2 - 1", "y - x")), windowed = FALSE, homo = FALSE
  )
  expect_match(code, "1, 0", fixed = TRUE)
  expect_match(code, "0, 1", fixed = TRUE)
})
