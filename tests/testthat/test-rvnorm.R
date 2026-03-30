skip_if_no_cmdstan <- function() {
  skip_if_not_installed("cmdstanr")
  ver <- tryCatch(cmdstanr::cmdstan_version(), error = function(e) NULL)
  if (is.null(ver) || identical(ver, "none")) skip("CmdStan is not available")
}

rvnorm_or_skip <- function(...) {
  tryCatch(
    rvnorm(...),
    error = function(e) {
      skip(paste("Skipping rvnorm sampling test:", conditionMessage(e)))
    }
  )
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
  samps1 <- rvnorm_or_skip(
    250, mp("x^2 + y^2 -1"), sd = 0.1, chains = 1, cores = 1
  )
  set.seed(123)
  samps2 <- rvnorm_or_skip(
    250, mp("x^2 + y^2 -1"), sd = 0.1, chains = 1, cores = 1
  )
  expect_equal(samps1, samps2)
})

test_that("rvnorm() has stated output", {
  skip_if_no_cmdstan()
  samps1 <- rvnorm_or_skip(
    250, mp("x^2 + y^2 -1"), sd = 0.1, chains = 1, cores = 1
  )
  expect_true(is.data.frame(samps1))
  samps2 <- rvnorm_or_skip(
    250,
    mp("x^2 + y^2 -1"),
    sd = 0.1,
    output = "tibble",
    chains = 1,
    cores = 1
  )
  expect_true(tibble::is_tibble(samps2))
  samps3 <- rvnorm_or_skip(
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

test_that("rvnorm(): pre_compiled == TRUE works", {
  skip_if_no_cmdstan()
  poly <- mp("x^2 + y^2 - 1")
  result_precompiled <- rvnorm_or_skip(
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
  result <- rvnorm_or_skip(
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
  expect_true(is.list(capture$args$data))
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
    "`Sigma` should be a scalar"
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

test_that("rvnorm() accepts character poly and valid Sigma shapes in code_only mode", {
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
      sd = 0.05,
      Sigma = c(1, 2),
      code_only = TRUE
    ),
    "character"
  )

  expect_type(
    rvnorm(
      n = 10,
      poly = mp("x^2 + y^2 - 1"),
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
    "`Sigma` should be a scalar"
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
  expect_true(is.matrix(capture$args$data$si))
  expect_equal(dim(capture$args$data$si), c(2, 2))
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
  expect_true(is.matrix(capture$args$data$si))
  expect_equal(capture$args$data$si, diag(c(0.5, 0.5)))
})
