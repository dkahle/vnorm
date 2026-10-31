skip_if_no_cmdstan <- function() {
  skip_if_not_installed("cmdstanr")
  ver <- tryCatch(cmdstanr::cmdstan_version(), error = function(e) NULL)
  if (is.null(ver) || identical(ver, "none")) skip("CmdStan is not available")
}

read_stan_data <- function(path) {
  jsonlite::read_json(path, simplifyVector = TRUE)
}

fake_cmdstan_fit <- function(draws_df) {
  structure(
    list(
      draws = function(format = "df", inc_warmup = FALSE) {
        draws_df
      }
    ),
    class = "CmdStanMCMC"
  )
}

fake_cmdstan_model <- function(draws_df, capture = NULL) {
  fit <- fake_cmdstan_fit(draws_df)
  exe <- tempfile()
  writeLines("fake", exe)
  list(
    sample = function(...) {
      if (!is.null(capture)) {
        capture$args <- list(...)
      }
      fit
    },
    exe_file = function() exe
  )
}

test_that("rvnorm() respects RNG seeds", {
  skip_if_no_cmdstan()
  set.seed(123)
  samps1 <- rvnorm(
    250, mp("x^2 + y^2 -1"), sd = 0.1, chains = 1, cores = 1
  )
  set.seed(123)
  samps2 <- rvnorm(
    250, mp("x^2 + y^2 -1"), sd = 0.1, chains = 1, cores = 1
  )
  expect_equal(samps1, samps2)
})

test_that("rvnorm() has stated output", {
  skip_if_no_cmdstan()
  samps1 <- rvnorm(
    250, mp("x^2 + y^2 -1"), sd = 0.1, chains = 1, cores = 1
  )
  expect_true(is.data.frame(samps1))
  samps2 <- rvnorm(
    250,
    mp("x^2 + y^2 -1"),
    sd = 0.1,
    output = "tibble",
    chains = 1,
    cores = 1
  )
  expect_true(tibble::is_tibble(samps2))
  samps3 <- rvnorm(
    250,
    mp("x^2 + y^2 -1"),
    sd = 0.1,
    output = "stanfit",
    chains = 1,
    cores = 1
  )
  expect_true(inherits(samps3, "CmdStanMCMC"))
})

test_that("rvnorm(): code_only = TRUE works", {
  skip_on_cran()
  stan_code <- rvnorm(
    1000, mp("x^2 + y^2 -1"), sd = 0.1, w = 5, code_only = TRUE
  )
  expected_result <- "data {
  real<lower=0> si;
}

parameters {
  real<lower=-5,upper=5> x;
  real<lower=-5,upper=5> y;
}

transformed parameters {
  real g = x^2+y^2-1;
  real ndg = sqrt(4*x^2+4*y^2);
}

model {
  target += normal_lpdf(0.00 | g/ndg, si);
}"
  expect_equal(stan_code, expected_result)
})

test_that("rvnorm(): generated code with large or long coefficients compiles", {
  skip_if_no_cmdstan()
  for (p in list(
    mp("x^4 + y^4 - 3000000000"),
    mp("x^4 + y^4") - 123456789012345.67,
    mp("x^4 + y^4") - 1e-20 / 3
  )) {
    stan_code <- rvnorm(10, p, sd = 0.1, code_only = TRUE)
    model <- cmdstanr::cmdstan_model(
      cmdstanr::write_stan_file(stan_code), compile = FALSE
    )
    expect_true(model$check_syntax(quiet = TRUE))
  }
})

test_that("rvnorm() generated code and pdvnorm() agree at full precision", {
  skip_if_no_cmdstan()
  # lp__ of the generated program is normal_lpdf(0 | g/ndg, si), which is
  # exactly pdvnorm(log = TRUE) when both use the same (unrounded) polynomial
  p <- mp("x - 1.23456789012345")
  samps <- rvnorm(
    200, p, sd = 1e-3, output = "tibble", pre_compiled = FALSE,
    chains = 1, cores = 1, seed = 1, sig_figs = 18
  )
  lp <- pdvnorm(samps$x, p, sd = 1e-3, log = TRUE)
  expect_lt(max(abs(samps$lp__ - lp)), 1e-8)
})

test_that("rvnorm(): pre_compiled == TRUE works", {
  skip_if_no_cmdstan()
  poly <- mp("x^2 + y^2 - 1")
  result_precompiled <- rvnorm(
    n = 10,
    poly = poly,
    sd = 0.05,
    pre_compiled = TRUE,
    chains = 1,
    cores = 1
  )
  expect_true(is.data.frame(result_precompiled))
})

test_that(
  paste0(
    "rvnorm() correctly replaces variables in Stan output ",
    "if poly has non-standard vars"
  ),
  {
  skip_if_no_cmdstan()
  poly <- mp("a^2 + b^2 - 1")
  result <- rvnorm(
    n = 10, poly = poly, sd = 0.05, chains = 1, cores = 1
  )
  expect_true(all(mpoly::vars(poly) %in% colnames(result)))
  }
)

test_that("rvnorm() errors for invalid poly input", {
  expect_error(
    rvnorm(n = 10, poly = 123, sd = 0.05),
    "`poly` should be either a character vector, mpoly, or mpolyList."
  )
})

test_that("rvnorm() rejection path returns samples without Stan", {
  set.seed(1)
  out <- rvnorm(
    n = 12,
    poly = mp("x^2 + y^2 - 1"),
    sd = 0.05,
    rejection = TRUE,
    output = "tibble",
    w = 1.25
  )

  expect_true(tibble::is_tibble(out))
  expect_equal(nrow(out), 12)
  expect_true(all(c("x", "y") %in% names(out)))
})

test_that("rvnorm() rejection path returns a data frame for simple output", {
  # matches the Stan path and the documented return value
  set.seed(1)
  out <- rvnorm(
    n = 12,
    poly = mp("x^2 + y^2 - 1"),
    sd = 0.05,
    rejection = TRUE,
    w = 1.25
  )

  expect_s3_class(out, "data.frame", exact = TRUE)
  expect_named(out, c("x", "y"))
  expect_equal(nrow(out), 12)
})

test_that("rvnorm() rejection path defaults to rejection_sampler()'s window", {
  # w = 1 clipped 18% of the unit circle's mass at sd = .05
  capture <- new.env(parent = emptyenv())
  testthat::local_mocked_bindings(
    rejection_sampler = function(...) {
      capture$args <- list(...)
      matrix(0, nrow = 2, ncol = 2, dimnames = list(NULL, c("x", "y")))
    },
    .package = "vnorm"
  )

  rvnorm(2, mp("x^2 + y^2 - 1"), sd = 0.05, rejection = TRUE)
  expect_equal(capture$args$w, 1.25)
})

test_that("rvnorm() rejection path validates sd shape", {
  expect_error(
    rvnorm(
      n = 10,
      poly = mp("x^2 + y^2 - 1"),
      sd = c(0.1, 0.2),
      rejection = TRUE
    ),
    "not supported yet for rejection sampler"
  )
})

test_that("rvnorm() errors when user_compiled cache is missing", {
  clear_compiled_stan_info()
  expect_error(
    rvnorm(
      n = 10,
      poly = mp("x^2 + y^2 - 1"),
      sd = 0.05,
      user_compiled = TRUE
    ),
    "No compiled model cache found"
  )
})

test_that("rvnorm() code_only works for mpolyList templates", {
  code <- rvnorm(
    n = 10,
    poly = mp(c("x^2 + y^2 - 1", "y")),
    sd = 0.05,
    code_only = TRUE,
    homo = FALSE
  )
  expect_type(code, "character")
  expect_true(length(code) == 1)
  expect_match(code, "target \\+=")
})

test_that("rvnorm() pre_compiled path works with mocked package model", {
  draws_df <- data.frame(
    lp__ = c(-1, -2), x = c(0.1, 0.2), y = c(0.3, 0.4), check.names = FALSE
  )
  capture <- new.env(parent = emptyenv())

  testthat::local_mocked_bindings(
    cmdstan_model = function(...) stop("cmdstan_model should not be called"),
    make_coefficients_data = function(...) list(fake_coef = 1),
    check_and_replace_vars = function(p) {
      list(polynomial = mp("x^2 + y^2 - 1"), mapping = list(x = "a", y = "b"))
    },
    rename_output_df = function(df, replacement_list) {
      names(df)[names(df) == "x"] <- replacement_list[["x"]]
      names(df)[names(df) == "y"] <- replacement_list[["y"]]
      df
    },
    .package = "vnorm"
  )
  testthat::local_mocked_bindings(
    stan_package_model = function(...) fake_cmdstan_model(draws_df, capture = capture),
    .package = "instantiate"
  )

  out <- rvnorm(
    n = 2,
    poly = mp("a^2 + b^2 - 1"),
    sd = 0.05,
    output = "tibble",
    pre_compiled = TRUE,
    chains = 1,
    cores = 1
  )

  expect_true(tibble::is_tibble(out))
  expect_equal(nrow(out), 2)
  expect_gte(ncol(out), 2)
  expect_equal(read_stan_data(capture$args$data), list(fake_coef = 1, si = 0.05))
})

test_that("rvnorm() tibble output puts variables before model quantities", {
  # precompiled templates also return g, dgx, dgy, and ndg, which were sorted
  # in among (for x/y, ahead of) the variable columns
  draws_df <- data.frame(
    lp__ = c(-1, -2), x = c(0.1, 0.2), y = c(0.3, 0.4), g = c(0, 0),
    dgx = c(1, 1), dgy = c(1, 1), ndg = c(1, 1),
    .chain = 1L, .iteration = 1:2, .draw = 1:2, check.names = FALSE
  )

  testthat::local_mocked_bindings(
    cmdstan_model = function(...) stop("cmdstan_model should not be called"),
    make_coefficients_data = function(...) list(fake_coef = 1),
    .package = "vnorm"
  )
  testthat::local_mocked_bindings(
    stan_package_model = function(...) fake_cmdstan_model(draws_df),
    .package = "instantiate"
  )

  model_cols <- c("g", "dgx", "dgy", "ndg", ".chain", ".iteration", ".draw", "lp__")
  out_xy <- rvnorm(
    n = 2, poly = mp("x^2 + y^2 - 1"), sd = 0.05, output = "tibble",
    chains = 1, cores = 1
  )
  out_ab <- rvnorm(
    n = 2, poly = mp("b^2 + a^2 - 1"), sd = 0.05, output = "tibble",
    chains = 1, cores = 1
  )

  expect_named(out_xy, c("x", "y", model_cols))
  expect_named(out_ab, c("a", "b", model_cols))
})

test_that("write_stan_data keeps every significant digit", {
  # cmdstanr's writer keeps 15 significant digits
  data <- list(
    b1 = -123456789012345.67, bx = 1 / 3,
    si = matrix(c(1, 0.1, 0.1, 2), 2), w = 1.25
  )
  expect_identical(read_stan_data(write_stan_data(data)), data)
})

test_that("rvnorm() passes precompiled coefficients at full precision", {
  draws_df <- data.frame(lp__ = c(-1, -2), x = c(0.1, 0.2), check.names = FALSE)
  capture <- new.env(parent = emptyenv())
  testthat::local_mocked_bindings(
    cmdstan_model = function(...) stop("cmdstan_model should not be called"),
    .package = "vnorm"
  )
  testthat::local_mocked_bindings(
    stan_package_model = function(...) fake_cmdstan_model(draws_df, capture = capture),
    .package = "instantiate"
  )

  a <- 123456789012345.67
  rvnorm(2, mp("x") - a, sd = 0.1, pre_compiled = TRUE, chains = 1, cores = 1)
  stan_data <- read_stan_data(capture$args$data)
  expect_identical(stan_data$b1, -a)
  expect_equal(stan_data$bx, 1)
})

test_that("rvnorm() keeps scalar w on precompiled path", {
  draws_df <- data.frame(
    lp__ = c(-1, -2), x = c(0.1, 0.2), y = c(0.3, 0.4), check.names = FALSE
  )
  capture <- new.env(parent = emptyenv())
  instantiate_called <- new.env(parent = emptyenv())
  instantiate_called$flag <- FALSE

  testthat::local_mocked_bindings(
    cmdstan_model = function(...) stop("cmdstan_model should not be called"),
    make_coefficients_data = function(...) list(fake_coef = 1),
    .package = "vnorm"
  )
  testthat::local_mocked_bindings(
    stan_package_model = function(...) {
      instantiate_called$flag <- TRUE
      fake_cmdstan_model(draws_df, capture = capture)
    },
    .package = "instantiate"
  )

  out <- rvnorm(
    n = 2,
    poly = mp("x^2 + y^2 - 1"),
    sd = 0.05,
    w = 1.5,
    pre_compiled = TRUE,
    chains = 1,
    cores = 1
  )

  expect_true(instantiate_called$flag)
  expect_true(is.data.frame(out))
  expect_equal(read_stan_data(capture$args$data)$w, 1.5)
})

test_that("rvnorm() uses generated Stan for non-scalar w", {
  draws_df <- data.frame(
    lp__ = c(-1, -2), x = c(0.1, 0.2), y = c(0.3, 0.4), check.names = FALSE
  )
  capture <- new.env(parent = emptyenv())
  instantiate_called <- new.env(parent = emptyenv())
  instantiate_called$flag <- FALSE

  testthat::local_mocked_bindings(
    create_stan_code = function(poly, sd, n_eqs, w, homo, vars) {
      capture$w <- w
      "fake_stan_code"
    },
    write_stan_file = function(...) "fake.stan",
    cmdstan_model = function(...) fake_cmdstan_model(draws_df, capture = capture),
    .package = "vnorm"
  )
  testthat::local_mocked_bindings(
    stan_package_model = function(...) {
      instantiate_called$flag <- TRUE
      stop("precompiled path should not be used")
    },
    .package = "instantiate"
  )

  expect_message(
    out <- rvnorm(
      n = 2,
      poly = mp("x^2 + y^2 - 1"),
      sd = 0.05,
      w = list(x = c(-1.5, 1.5)),
      pre_compiled = TRUE,
      verbose = TRUE,
      chains = 1,
      cores = 1
    ),
    "non-scalar `w` requires generated Stan code"
  )

  expect_false(instantiate_called$flag)
  expect_true(is.data.frame(out))
  expect_equal(capture$w, list(x = c(-1.5, 1.5)))
  expect_false("w" %in% names(read_stan_data(capture$args$data)))
})

test_that("rvnorm() falls back when pre-compiled model is unavailable", {
  draws_df <- data.frame(
    lp__ = c(-1, -2), x = c(0.1, 0.2), y = c(0.3, 0.4), check.names = FALSE
  )
  capture <- new.env(parent = emptyenv())

  testthat::local_mocked_bindings(
    create_stan_code = function(...) "fake_stan_code",
    write_stan_file = function(...) "fake.stan",
    cmdstan_model = function(...) fake_cmdstan_model(draws_df, capture = capture),
    .package = "vnorm"
  )
  testthat::local_mocked_bindings(
    stan_package_model = function(...) stop("Stan model file \"\" not found."),
    .package = "instantiate"
  )

  expect_message(
    out <- rvnorm(
      n = 2,
      poly = mp("x^2 + y^2 - 1"),
      sd = 0.05,
      pre_compiled = TRUE,
      chains = 1,
      cores = 1
    ),
    "falling back to regular rvnorm sampling"
  )

  expect_true(is.data.frame(out))
  expect_equal(nrow(out), 2)
  expect_equal(capture$args$iter_sampling, 2)
})

test_that("rvnorm() pre_compiled auto-disables for mpolyList and uses regular path", {
  draws_df <- data.frame(
    lp__ = c(-1, -2), x = c(0.1, 0.2), y = c(0.3, 0.4), check.names = FALSE
  )
  instantiate_called <- new.env(parent = emptyenv())
  instantiate_called$flag <- FALSE

  testthat::local_mocked_bindings(
    create_stan_code = function(...) "fake_stan_code",
    write_stan_file = function(...) "fake.stan",
    cmdstan_model = function(...) fake_cmdstan_model(draws_df),
    .package = "vnorm"
  )
  testthat::local_mocked_bindings(
    stan_package_model = function(...) {
      instantiate_called$flag <- TRUE
      stop("should not be called")
    },
    .package = "instantiate"
  )

  out <- rvnorm(
    n = 2,
    poly = mp(c("x^2 + y^2 - 1", "y")),
    sd = 0.05,
    pre_compiled = TRUE,
    chains = 1,
    cores = 1
  )

  expect_false(instantiate_called$flag)
  expect_true(is.data.frame(out))
})

test_that("rvnorm() user_compiled path works with cached model metadata", {
  draws_df <- data.frame(
    lp__ = c(-1, -2), x = c(0.1, 0.2), y = c(0.3, 0.4), check.names = FALSE
  )
  model_name <- generate_model_name(mp("x^2 + y^2 - 1"), homo = TRUE, windowed = FALSE)
  set_compiled_stan_info(data.frame(
    name = model_name, path = "cached-model.stan", stringsAsFactors = FALSE
  ))
  on.exit(clear_compiled_stan_info(), add = TRUE)

  testthat::local_mocked_bindings(
    cmdstan_model = function(...) fake_cmdstan_model(draws_df),
    get_coefficients_data = function(...) list(fake_coef = 1),
    .package = "vnorm"
  )

  out <- rvnorm(
    n = 2,
    poly = mp("x^2 + y^2 - 1"),
    sd = 0.05,
    user_compiled = TRUE,
    chains = 1,
    cores = 1
  )

  expect_true(is.data.frame(out))
  expect_equal(nrow(out), 2)
})

test_that("rvnorm() errors when requested user_compiled model is missing", {
  set_compiled_stan_info(data.frame(
    name = "other_model", path = "cached-model.stan", stringsAsFactors = FALSE
  ))
  on.exit(clear_compiled_stan_info(), add = TRUE)

  expect_error(
    rvnorm(
      n = 2,
      poly = mp("x^2 + y^2 - 1"),
      sd = 0.05,
      user_compiled = TRUE
    ),
    "Requested compiled model not found in cache"
  )
})

test_that("rvnorm() validates Sigma dimensions before sampling", {
  expect_error(
    rvnorm(
      n = 10,
      poly = mp("x^2 + y^2 - 1"),
      sd = 0.05,
      Sigma = c(1, 2, 3),
      code_only = TRUE
    ),
    "single positive variance"
  )
})

test_that("rvnorm() computes refresh from n when verbose=TRUE and refresh is default", {
  draws_df <- data.frame(
    lp__ = seq_len(20),
    x = seq(0.1, 2.0, length.out = 20),
    y = seq(0.2, 2.1, length.out = 20),
    check.names = FALSE
  )
  capture <- new.env(parent = emptyenv())

  testthat::local_mocked_bindings(
    create_stan_code = function(...) "fake_stan_code",
    write_stan_file = function(...) "fake.stan",
    cmdstan_model = function(...) fake_cmdstan_model(draws_df, capture = capture),
    .package = "vnorm"
  )

  rvnorm(
    n = 20,
    poly = mp(c("x^2 + y^2 - 1", "y")),
    sd = 0.05,
    pre_compiled = TRUE,
    verbose = TRUE,
    chains = 2,
    cores = 1
  )

  expect_true(is.numeric(capture$args$refresh))
  expect_gte(capture$args$refresh, 1)
})

test_that("rvnorm() preserves user-supplied refresh value", {
  draws_df <- data.frame(
    lp__ = seq_len(20),
    x = seq(0.1, 2.0, length.out = 20),
    y = seq(0.2, 2.1, length.out = 20),
    check.names = FALSE
  )
  capture <- new.env(parent = emptyenv())

  testthat::local_mocked_bindings(
    create_stan_code = function(...) "fake_stan_code",
    write_stan_file = function(...) "fake.stan",
    cmdstan_model = function(...) fake_cmdstan_model(draws_df, capture = capture),
    .package = "vnorm"
  )

  rvnorm(
    n = 20,
    poly = mp(c("x^2 + y^2 - 1", "y")),
    sd = 0.05,
    pre_compiled = TRUE,
    refresh = 500,
    verbose = TRUE,
    chains = 2,
    cores = 1
  )

  expect_equal(capture$args$refresh, 500)
})

test_that("rvnorm() accepts character poly and narrows single-polynomial Sigma", {
  expect_type(
    rvnorm(
      n = 10,
      poly = "x^2 + y^2 - 1",
      sd = 0.05,
      code_only = TRUE
    ),
    "character"
  )

  expect_type(
    rvnorm(
      n = 10,
      poly = mp("x^2 + y^2 - 1"),
      Sigma = 0.25,
      code_only = TRUE
    ),
    "character"
  )

  expect_error(
    rvnorm(
      n = 10,
      poly = mp("x^2 + y^2 - 1"),
      sd = 0.05,
      Sigma = c(1, 2),
      code_only = TRUE
    ),
    "single positive variance"
  )

  expect_error(
    rvnorm(
      n = 10,
      poly = mp("x^2 + y^2 - 1"),
      sd = 0.05,
      Sigma = diag(c(1, 2)),
      code_only = TRUE
    ),
    "single positive variance"
  )
})

test_that("rvnorm() accepts mpolyList Sigma shapes in code_only mode", {
  expect_type(
    rvnorm(
      n = 10,
      poly = mp(c("x^2 + y^2 - 1", "y")),
      sd = 0.05,
      Sigma = c(1, 2),
      code_only = TRUE
    ),
    "character"
  )

  expect_type(
    rvnorm(
      n = 10,
      poly = mp(c("x^2 + y^2 - 1", "y")),
      sd = 0.05,
      Sigma = diag(c(1, 2)),
      code_only = TRUE
    ),
    "character"
  )
})

test_that("rvnorm() emits stanfit rename warning in precompiled path", {
  draws_df <- data.frame(
    lp__ = c(-1, -2), x = c(0.1, 0.2), y = c(0.3, 0.4), check.names = FALSE
  )

  testthat::local_mocked_bindings(
    make_coefficients_data = function(...) list(fake_coef = 1),
    check_and_replace_vars = function(p) {
      list(polynomial = mp("x^2 + y^2 - 1"), mapping = list(x = "a", y = "b"))
    },
    .package = "vnorm"
  )
  testthat::local_mocked_bindings(
    stan_package_model = function(...) fake_cmdstan_model(draws_df),
    .package = "instantiate"
  )

  expect_message(
    out <- rvnorm(
      n = 2,
      poly = mp("a^2 + b^2 - 1"),
      sd = 0.05,
      output = "stanfit",
      pre_compiled = TRUE,
      chains = 1,
      cores = 1
    ),
    "variable names as x,y,z"
  )
  expect_true(inherits(out, "CmdStanMCMC"))
})

test_that("rvnorm() returns all draws when inc_warmup = TRUE", {
  warmup_df <- data.frame(
    lp__ = seq_len(6),
    x = seq(0.1, 0.6, length.out = 6),
    y = seq(0.2, 0.7, length.out = 6),
    check.names = FALSE
  )
  post_df <- data.frame(
    lp__ = seq_len(4),
    x = seq(0.7, 1.0, length.out = 4),
    y = seq(0.8, 1.1, length.out = 4),
    check.names = FALSE
  )
  combined_df <- rbind(warmup_df, post_df)

  fake_fit <- structure(
    list(
      draws = function(format = "df", inc_warmup = FALSE) {
        if (inc_warmup) combined_df else post_df
      }
    ),
    class = "CmdStanMCMC"
  )
  exe <- tempfile()
  writeLines("fake", exe)
  fake_model <- list(
    sample = function(...) fake_fit,
    exe_file = function() exe
  )

  testthat::local_mocked_bindings(
    create_stan_code = function(...) "fake_stan_code",
    write_stan_file = function(...) "fake.stan",
    cmdstan_model = function(...) fake_model,
    .package = "vnorm"
  )

  out_with <- rvnorm(
    n = 4,
    poly = mp(c("x^2 + y^2 - 1", "y")),
    sd = 0.05,
    output = "tibble",
    inc_warmup = TRUE,
    pre_compiled = TRUE,
    chains = 1,
    cores = 1
  )
  expect_equal(nrow(out_with), nrow(combined_df))

  out_without <- rvnorm(
    n = 4,
    poly = mp(c("x^2 + y^2 - 1", "y")),
    sd = 0.05,
    output = "tibble",
    inc_warmup = FALSE,
    pre_compiled = TRUE,
    chains = 1,
    cores = 1
  )
  expect_equal(nrow(out_without), 4)
})

test_that("rvnorm() errors for non-positive n", {
  expect_error(rvnorm(n = 0, poly = mp("x"), sd = 0.1), "`n` must be a positive integer")
  expect_error(rvnorm(n = -5, poly = mp("x"), sd = 0.1), "`n` must be a positive integer")
})

test_that("rvnorm() errors for non-integer n", {
  expect_error(rvnorm(n = 2.5, poly = mp("x"), sd = 0.1), "`n` must be a positive integer")
  expect_error(rvnorm(n = "a", poly = mp("x"), sd = 0.1), "`n` must be a positive integer")
})

test_that("rvnorm() Sigma dimension mismatch with matrix", {
  expect_error(
    rvnorm(n = 10, poly = mp("x^2 + y^2 - 1"), sd = 0.05, Sigma = matrix(1:9, 3, 3)),
    "single positive variance"
  )
})

test_that("rvnorm() validates output, missing scale, chains, warmup, and thin", {
  expect_error(
    rvnorm(n = 2, poly = mp("x"), sd = 0.1, output = "bad", code_only = TRUE),
    "'arg' should be one of"
  )
  expect_error(
    rvnorm(n = 2, poly = mp("x"), code_only = TRUE),
    "`sd` must be supplied"
  )
  expect_error(
    rvnorm(n = 2, poly = mp("x"), sd = 0.1, chains = 0, code_only = TRUE),
    "`chains` must be a positive integer"
  )
  expect_error(
    rvnorm(n = 2, poly = mp("x"), sd = 0.1, warmup = -1, code_only = TRUE),
    "`warmup` must be a non-negative integer"
  )
  expect_error(
    rvnorm(n = 2, poly = mp("x"), sd = 0.1, thin = 0, code_only = TRUE),
    "`thin` must be a positive integer"
  )
})

test_that("rvnorm() seed parameter is forwarded to Stan sampler", {
  draws_df <- data.frame(
    lp__ = seq_len(4), x = seq(0.1, 0.4, length.out = 4),
    y = seq(0.2, 0.5, length.out = 4), check.names = FALSE
  )
  capture <- new.env(parent = emptyenv())

  testthat::local_mocked_bindings(
    create_stan_code = function(...) "fake_stan_code",
    write_stan_file = function(...) "fake.stan",
    cmdstan_model = function(...) fake_cmdstan_model(draws_df, capture = capture),
    .package = "vnorm"
  )

  rvnorm(
    n = 4, poly = mp(c("x^2 + y^2 - 1", "y")), sd = 0.05,
    pre_compiled = TRUE, chains = 1, cores = 1, seed = 42
  )
  expect_equal(capture$args$seed, 42)
})

test_that("rvnorm() thin parameter is forwarded to Stan sampler", {
  draws_df <- data.frame(
    lp__ = seq_len(4), x = seq(0.1, 0.4, length.out = 4),
    y = seq(0.2, 0.5, length.out = 4), check.names = FALSE
  )
  capture <- new.env(parent = emptyenv())

  testthat::local_mocked_bindings(
    create_stan_code = function(...) "fake_stan_code",
    write_stan_file = function(...) "fake.stan",
    cmdstan_model = function(...) fake_cmdstan_model(draws_df, capture = capture),
    .package = "vnorm"
  )

  rvnorm(
    n = 4, poly = mp(c("x^2 + y^2 - 1", "y")), sd = 0.05,
    pre_compiled = TRUE, chains = 1, cores = 1, thin = 2
  )
  expect_equal(capture$args$thin, 2)
  # n counts returned draws, so thinning by 2 doubles the iterations
  expect_equal(capture$args$iter_sampling, 8)
})

test_that("rvnorm() extra ... args are forwarded to Stan sampler", {
  draws_df <- data.frame(
    lp__ = seq_len(4), x = seq(0.1, 0.4, length.out = 4),
    y = seq(0.2, 0.5, length.out = 4), check.names = FALSE
  )
  capture <- new.env(parent = emptyenv())

  testthat::local_mocked_bindings(
    create_stan_code = function(...) "fake_stan_code",
    write_stan_file = function(...) "fake.stan",
    cmdstan_model = function(...) fake_cmdstan_model(draws_df, capture = capture),
    .package = "vnorm"
  )

  rvnorm(
    n = 4, poly = mp(c("x^2 + y^2 - 1", "y")), sd = 0.05,
    pre_compiled = TRUE, chains = 1, cores = 1, init = 0.5
  )
  expect_equal(capture$args$init, 0.5)
})

test_that("rvnorm() forwards sampler control defaults and overrides", {
  draws_df <- data.frame(
    lp__ = seq_len(4), x = seq(0.1, 0.4, length.out = 4),
    y = seq(0.2, 0.5, length.out = 4), check.names = FALSE
  )
  capture <- new.env(parent = emptyenv())

  testthat::local_mocked_bindings(
    create_stan_code = function(...) "fake_stan_code",
    write_stan_file = function(...) "fake.stan",
    cmdstan_model = function(...) fake_cmdstan_model(draws_df, capture = capture),
    .package = "vnorm"
  )

  rvnorm(
    n = 4, poly = mp(c("x^2 + y^2 - 1", "y")), sd = 0.05,
    pre_compiled = TRUE, chains = 1, cores = 1
  )
  expect_equal(capture$args$adapt_delta, .999)
  expect_equal(capture$args$max_treedepth, 20L)

  rvnorm(
    n = 4, poly = mp(c("x^2 + y^2 - 1", "y")), sd = 0.05,
    pre_compiled = TRUE, chains = 1, cores = 1,
    adapt_delta = .9, max_treedepth = 12
  )
  expect_equal(capture$args$adapt_delta, .9)
  expect_equal(capture$args$max_treedepth, 12L)
})

test_that("rvnorm() default warmup scales with per-chain draws", {
  draws_df <- data.frame(
    lp__ = seq_len(4), x = seq(0.1, 0.4, length.out = 4),
    y = seq(0.2, 0.5, length.out = 4), check.names = FALSE
  )
  capture <- new.env(parent = emptyenv())

  testthat::local_mocked_bindings(
    create_stan_code = function(...) "fake_stan_code",
    write_stan_file = function(...) "fake.stan",
    cmdstan_model = function(...) fake_cmdstan_model(draws_df, capture = capture),
    .package = "vnorm"
  )

  rvnorm(
    n = 10000, poly = mp(c("x^2 + y^2 - 1", "y")), sd = 0.05,
    pre_compiled = TRUE, chains = 4, cores = 1
  )
  expect_equal(capture$args$iter_warmup, 1250)
})

test_that("rvnorm() show_messages parameter is forwarded", {
  draws_df <- data.frame(
    lp__ = seq_len(4), x = seq(0.1, 0.4, length.out = 4),
    y = seq(0.2, 0.5, length.out = 4), check.names = FALSE
  )
  capture <- new.env(parent = emptyenv())

  testthat::local_mocked_bindings(
    create_stan_code = function(...) "fake_stan_code",
    write_stan_file = function(...) "fake.stan",
    cmdstan_model = function(...) fake_cmdstan_model(draws_df, capture = capture),
    .package = "vnorm"
  )

  rvnorm(
    n = 4, poly = mp(c("x^2 + y^2 - 1", "y")), sd = 0.05,
    pre_compiled = TRUE, chains = 1, cores = 1, show_messages = TRUE
  )
  expect_true(capture$args$show_messages)
})

test_that("rvnorm() handles character poly input through sampling path", {
  draws_df <- data.frame(
    lp__ = c(-1, -2), x = c(0.7, 0.8), y = c(0.7, 0.6), check.names = FALSE
  )

  testthat::local_mocked_bindings(
    create_stan_code = function(...) "fake_stan_code",
    write_stan_file = function(...) "fake.stan",
    cmdstan_model = function(...) fake_cmdstan_model(draws_df),
    .package = "vnorm"
  )

  out <- rvnorm(
    n = 2, poly = "x^2 + y^2 - 1", sd = 0.05,
    pre_compiled = TRUE, chains = 1, cores = 1
  )
  expect_true(is.data.frame(out))
  expect_equal(nrow(out), 2)
  expect_true(all(c("x", "y") %in% names(out)))
})

test_that("rvnorm() homo=FALSE through sampling path", {
  draws_df <- data.frame(
    lp__ = c(-1, -2), x = c(0.7, 0.8), y = c(0.7, 0.6), check.names = FALSE
  )

  testthat::local_mocked_bindings(
    create_stan_code = function(...) "fake_stan_code",
    write_stan_file = function(...) "fake.stan",
    cmdstan_model = function(...) fake_cmdstan_model(draws_df),
    .package = "vnorm"
  )

  out <- rvnorm(
    n = 2, poly = mp("x^2 + y^2 - 1"), sd = 0.05,
    homo = FALSE, pre_compiled = TRUE, chains = 1, cores = 1
  )
  expect_true(is.data.frame(out))
  expect_equal(nrow(out), 2)
})

test_that("rvnorm() rejects dependent equations before building Stan models", {
  # the homoskedastic programs solve with J, J*J', or J'*J, which are
  # singular everywhere for dependent equations: no chain could start
  draws_df <- data.frame(
    lp__ = c(-1, -2), x = c(0.1, 0.2), y = c(0.3, 0.4), z = c(0, 0),
    check.names = FALSE
  )
  testthat::local_mocked_bindings(
    create_stan_code = function(...) "fake_stan_code",
    write_stan_file = function(...) "fake.stan",
    cmdstan_model = function(...) fake_cmdstan_model(draws_df),
    .package = "vnorm"
  )

  for (p in list(
    mp(c("x + y + z", "(x + y + z)^2")),
    mp(c("x + y", "2 x + 2 y")),
    mp(c("x + y", "2 x + 2 y", "x + y + 1"))
  )) {
    expect_error(
      rvnorm(2, p, sd = 0.1, chains = 1, cores = 1),
      "loses rank everywhere"
    )
  }

  # homo = FALSE never uses J
  out <- rvnorm(
    2, mp(c("x + y + z", "(x + y + z)^2")), sd = 0.1, homo = FALSE,
    chains = 1, cores = 1
  )
  expect_equal(nrow(out), 2)

  # badly scaled but independent: J = diag(100 x^99, 100 y^99) is tiny at the
  # old check's fixed points, and solve() refused the second one's J'J there
  for (p in list(
    mp(c("x^100 - 1", "y^100 - 1")),
    mp(c("0.000000001 x^2 - 1", "y^2 - 1", "0.000000001 x^2 + y^2 - 2"))
  )) {
    out <- rvnorm(2, p, sd = 0.1, chains = 1, cores = 1)
    expect_equal(nrow(out), 2)
  }

  # rescaled by 1e-20: elimination in doubles underflowed to 0 and refused
  # these independent equations
  vars <- paste0("x", 1:5)
  xs <- lapply(vars, mp)
  s <- Reduce(`+`, xs)
  q <- 1e-20 * structure(lapply(xs, function(x) s + 1e-8 * x), class = "mpolyList")
  draws5 <- data.frame(lp__ = c(-1, -2), matrix(0.1, 2, 5, dimnames = list(NULL, vars)))
  testthat::local_mocked_bindings(
    cmdstan_model = function(...) fake_cmdstan_model(draws5),
    .package = "vnorm"
  )
  out <- rvnorm(2, q, sd = 0.1, chains = 1, cores = 1)
  expect_equal(nrow(out), 2)
})

test_that("rvnorm() full non-diagonal Sigma passed to Stan data", {
  draws_df <- data.frame(
    lp__ = c(-1, -2), x = c(0.7, 0.8), y = c(0.7, 0.6), check.names = FALSE
  )
  capture <- new.env(parent = emptyenv())

  testthat::local_mocked_bindings(
    create_stan_code = function(...) "fake_stan_code",
    write_stan_file = function(...) "fake.stan",
    cmdstan_model = function(...) fake_cmdstan_model(draws_df, capture = capture),
    .package = "vnorm"
  )

  S <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
  rvnorm(
    n = 2, poly = mp(c("x^2 + y^2 - 1", "y")), Sigma = S, sd = 0.05,
    pre_compiled = TRUE, chains = 1, cores = 1
  )
  si <- read_stan_data(capture$args$data)$si
  expect_true(is.matrix(si))
  expect_equal(si, S)
})

test_that("rvnorm() n < chains edge case sets iter_sampling to 1", {
  draws_df <- data.frame(
    lp__ = seq_len(4), x = seq(0.1, 0.4, length.out = 4),
    y = seq(0.2, 0.5, length.out = 4), check.names = FALSE
  )
  capture <- new.env(parent = emptyenv())

  testthat::local_mocked_bindings(
    create_stan_code = function(...) "fake_stan_code",
    write_stan_file = function(...) "fake.stan",
    cmdstan_model = function(...) fake_cmdstan_model(draws_df, capture = capture),
    .package = "vnorm"
  )

  rvnorm(
    n = 1, poly = mp(c("x^2 + y^2 - 1", "y")), sd = 0.05,
    pre_compiled = TRUE, chains = 4, cores = 1
  )
  expect_equal(capture$args$iter_sampling, 1)
})

test_that("rvnorm() Sigma as diagonal vector is converted to matrix", {
  draws_df <- data.frame(
    lp__ = c(-1, -2), x = c(0.7, 0.8), y = c(0.7, 0.6), check.names = FALSE
  )
  capture <- new.env(parent = emptyenv())

  testthat::local_mocked_bindings(
    create_stan_code = function(...) "fake_stan_code",
    write_stan_file = function(...) "fake.stan",
    cmdstan_model = function(...) fake_cmdstan_model(draws_df, capture = capture),
    .package = "vnorm"
  )

  rvnorm(
    n = 2, poly = mp(c("x^2 + y^2 - 1", "y")), sd = 0.05, Sigma = c(0.5, 0.5),
    pre_compiled = TRUE, chains = 1, cores = 1
  )
  si <- read_stan_data(capture$args$data)$si
  expect_true(is.matrix(si))
  expect_equal(si, diag(c(0.5, 0.5)))
})

test_that("rvnorm() heteroskedastic mpolyList Sigma uses equation dimension", {
  code <- rvnorm(
    n = 2,
    poly = mp(c("x", "y", "x + y")),
    sd = c(0.1, 0.2, 0.3),
    homo = FALSE,
    code_only = TRUE
  )
  expect_match(code, "cov_matrix\\[3\\] si")
  expect_match(code, "vector\\[3\\] g")
  expect_false(grepl("matrix[3,2] J", code, fixed = TRUE))

  expect_error(
    rvnorm(
      n = 2,
      poly = mp(c("x", "y", "x + y")),
      Sigma = diag(2),
      homo = FALSE,
      code_only = TRUE
    ),
    "number of polynomials"
  )
})

test_that("rvnorm() rejection path rejects Sigma explicitly", {
  expect_error(
    rvnorm(
      n = 2,
      poly = mp("x"),
      sd = 0.1,
      Sigma = 1,
      rejection = TRUE
    ),
    "supports `sd`, not `Sigma`"
  )
})

test_that("rvnorm() rejects constant polynomials", {
  expect_error(rvnorm(10, mp("3"), sd = .1), "at least one variable")
})

test_that("tibble output moves lp__ to the last column by name", {
  draws_df <- data.frame(
    lp__ = c(-1, -2), x = c(0.1, 0.2), y = c(0.3, 0.4), check.names = FALSE
  )

  testthat::local_mocked_bindings(
    create_stan_code = function(...) "// stan code",
    write_stan_file = function(...) "fake.stan",
    cmdstan_model = function(...) fake_cmdstan_model(draws_df),
    .package = "vnorm"
  )

  out <- rvnorm(
    n = 2,
    poly = mp("x^2 + y^2 - 1"),
    sd = 0.05,
    output = "tibble",
    pre_compiled = FALSE,
    chains = 1,
    cores = 1
  )

  expect_equal(names(out), c("x", "y", "lp__"))
})

test_that("Stan-path output columns follow sorted variable order", {
  draws_df <- data.frame(
    lp__ = c(-1, -2), y = c(0.3, 0.4), x = c(0.1, 0.2), check.names = FALSE
  )

  testthat::local_mocked_bindings(
    create_stan_code = function(...) "// stan code",
    write_stan_file = function(...) "fake.stan",
    cmdstan_model = function(...) fake_cmdstan_model(draws_df),
    .package = "vnorm"
  )

  out_simple <- rvnorm(
    n = 2, poly = mp("y^2 - x"), sd = .1,
    pre_compiled = FALSE, chains = 1, cores = 1
  )
  expect_equal(names(out_simple), c("x", "y"))

  out_tbl <- rvnorm(
    n = 2, poly = mp("y^2 - x"), sd = .1, output = "tibble",
    pre_compiled = FALSE, chains = 1, cores = 1
  )
  expect_equal(names(out_tbl), c("x", "y", "lp__"))
})
