skip_if_no_cmdstan <- function() {
  skip_if_not_installed("cmdstanr")
  ver <- tryCatch(cmdstanr::cmdstan_version(), error = function(e) NULL)
  if (is.null(ver) || identical(ver, "none")) skip("CmdStan is not available")
}

near_circle <- function(n, chains = 2L, sd = 0.01, seed = 1L, z = NULL) {
  # pilot-like draws near the unit circle, n / chains per chain
  withr::with_seed(seed, {
    th <- stats::runif(n, 0, 2 * pi)
    r <- 1 + stats::rnorm(n, sd = sd)
    out <- data.frame(x = r * cos(th), y = r * sin(th))
    if (!is.null(z)) out$z <- stats::runif(n, z[1], z[2])
    out$.chain <- rep(seq_len(chains), each = n / chains)
    out$.iteration <- rep(seq_len(n / chains), chains)
    out
  })
}

mock_pilot <- function(draws) {
  # replace rvnorm_domain_pilot(), recording its arguments; draws is a data
  # frame, or a function of the call number returning one
  calls <- new.env(parent = emptyenv())
  calls$args <- list()
  testthat::local_mocked_bindings(
    rvnorm_domain_pilot = function(n, poly, pilot_sd, rejection, seed, pilot_args) {
      calls$args[[length(calls$args) + 1L]] <- list(
        n = n, pilot_sd = pilot_sd, rejection = rejection, seed = seed,
        pilot_args = pilot_args
      )
      if (is.function(draws)) draws(length(calls$args)) else draws
    },
    .package = "vnorm",
    .env = parent.frame()
  )
  calls
}

circle <- mp("x^2 + y^2 - 1")

test_that("rvnorm_domain_covariance() builds the noise covariance over the variables", {
  vars <- c("x", "y", "z")
  expect_equal(rvnorm_domain_covariance(0.1, NULL, vars), diag(0.01, 3))
  expect_equal(
    rvnorm_domain_covariance(c(0.1, 0.2, 0.3), NULL, vars),
    diag(c(0.01, 0.04, 0.09))
  )
  expect_equal(
    rvnorm_domain_covariance(c(z = 0.3, x = 0.1, y = 0.2), NULL, vars),
    diag(c(0.01, 0.04, 0.09))
  )
  expect_equal(rvnorm_domain_covariance(NULL, 0.04, vars), diag(0.04, 3))
  expect_equal(rvnorm_domain_covariance(NULL, c(y = 2, x = 1, z = 3), vars), diag(c(1, 2, 3)))
  S <- matrix(
    c(2, 0.5, 0, 0.5, 1, 0, 0, 0, 3), 3,
    dimnames = list(c("y", "x", "z"), c("y", "x", "z"))
  )
  expect_equal(
    rvnorm_domain_covariance(NULL, S, vars),
    matrix(c(1, 0.5, 0, 0.5, 2, 0, 0, 0, 3), 3)
  )

  expect_error(rvnorm_domain_covariance(c(0.1, 0.2), NULL, vars), "Length of `sd`")
  expect_error(
    rvnorm_domain_covariance(c(a = 0.1, b = 0.2, c = 0.3), NULL, vars),
    "Names of `sd` must match"
  )
  expect_error(rvnorm_domain_covariance(diag(3), NULL, vars), "Use `Sigma`")
  expect_error(rvnorm_domain_covariance(NULL, NULL, vars), "`sd` must be supplied")
  expect_error(rvnorm_domain_covariance(NULL, diag(2), vars), "dimensions")
  # chol() alone would accept this, reading only the upper triangle
  expect_error(
    rvnorm_domain_covariance(NULL, matrix(c(1, 0.9, 0, 1), 2), c("x", "y")),
    "symmetric"
  )
  expect_error(
    rvnorm_domain_covariance(NULL, matrix(c(1, 2, 2, 1), 2), c("x", "y")),
    "positive definite"
  )
})

test_that("project_rows_closest() finds the closest points on the variety", {
  p <- mp("0.25 x^2 + y^2 - 1")
  x <- rbind(c(1.9, 0.4), c(-0.3, -1.08), c(1.2, 0.75), c(-2.05, 0.02))
  proj <- project_rows_closest(x, p, c("x", "y"))
  u <- proj$points
  expect_true(all(proj$ok))
  expect_lt(max(abs(poly_function(p, c("x", "y"))(u))), 1e-12)
  # x - u is normal to the ellipse, parallel to the gradient (x / 2, 2 y);
  # Gauss-Newton alone leaves about 1e-3 here
  r <- x - u
  expect_lt(max(abs(r[, 1] * 2 * u[, 2] - r[, 2] * u[, 1] / 2)), 1e-8)
  expect_equal(u, suppressWarnings(project_onto_variety(x, p)), tolerance = 1e-8)

  # a system: the unit circle in the plane z = 0
  x3 <- rbind(c(1.1, 0.2, 0.3), c(-0.5, 0.9, -0.2))
  proj3 <- project_rows_closest(
    x3, mp(c("x^2 + y^2 + z^2 - 1", "z")), c("x", "y", "z")
  )
  expect_true(all(proj3$ok))
  expect_equal(
    proj3$points, cbind(x3[, 1:2] / sqrt(rowSums(x3[, 1:2]^2)), 0),
    tolerance = 1e-10
  )

  # variables the polynomial does not involve stay put
  proj_z <- project_rows_closest(cbind(x, 0.7), p, c("x", "y", "z"))
  expect_equal(proj_z$points[, 3], rep(0.7, 4))
  expect_equal(proj_z$points[, 1:2], u, tolerance = 1e-10)
})

test_that("project_rows_closest() accepts only points where every equation vanishes", {
  # the variety of this overdetermined system is {1}; near -1 the
  # least-squares Gauss-Newton step vanishes at a point where g does not
  over <- mp(c("x^2 - 1", "0.1 x - 0.1"))
  proj <- project_rows_closest(matrix(c(-0.9, 0.9), ncol = 1), over, "x")
  expect_identical(proj$ok, c(FALSE, TRUE))
  expect_equal(proj$points[2, 1], 1)

  # an empty variety has no points to project to
  empty <- mp(c("x", "x - 1"))
  expect_false(any(project_rows_closest(matrix(c(0.3, 0.6), ncol = 1), empty, "x")$ok))
})

test_that("rvnorm(noise = \"domain\") places no points off the variety", {
  over <- mp(c("x^2 - 1", "0.1 x - 0.1"))
  expect_warning(
    z <- rvnorm(
      100, over, sd = 1, pilot_sd = 0.1, noise = "domain", rejection = TRUE,
      output = "tibble", w = 1.5, seed = 43
    ),
    "pilot draws were dropped"
  )
  expect_true(all(abs(z$x_proj - 1) < 1e-8))

  expect_error(
    rvnorm(
      8, mp(c("x", "x - 1")), sd = 0.2, noise = "domain", rejection = TRUE,
      w = 1, seed = 5
    ),
    "Only 0 of 8 pilot draws projected onto the variety"
  )
})

test_that("translating a system leaves the projection's tolerances alone", {
  # tolerances grew with |u|, so x - 1e9 and x - 1e9 - 1, an empty
  # variety, passed at 1e9 + 0.5
  expect_error(
    rvnorm(
      8, mp(c("x - 1000000000", "x - 1000000001")), sd = 0.2,
      noise = "domain", rejection = TRUE, w = 1e9 + c(-1, 2), seed = 5
    ),
    "Only 0 of 8 pilot draws projected onto the variety"
  )
  # as did a free coordinate far away
  expect_error(
    rvnorm(
      8, mp(c("x", "x - 1")), sd = 0.2, noise = "domain", rejection = TRUE,
      vars = c("x", "z"), w = list(x = c(-1, 2), z = 1e9 + c(-1, 1)), seed = 5
    ),
    "Only 0 of 8 pilot draws projected onto the variety"
  )

  # a translated variety with points still gets them
  z <- rvnorm(
    20, mp("x - 1000000000"), sd = 0.2, noise = "domain", rejection = TRUE,
    output = "tibble", w = 1e9 + c(-1, 1), seed = 6
  )
  expect_true(all(abs(z$x_proj - 1e9) < 1e-6))
})

test_that("free variables leave a system's geometry alone", {
  # x and x^2 cut out the line x = 0, and z is free. Counting z, the system
  # looked square, so its Jacobian's rank of 1 read as dependence
  p <- mp(c("x", "x^2"))
  proj <- project_rows_closest(rbind(c(0.05, 0.3), c(-0.02, -0.7)), p, c("x", "z"))
  expect_true(all(proj$ok))
  expect_equal(proj$points, rbind(c(0, 0.3), c(0, -0.7)))

  z <- rvnorm(
    50, p, sd = 0.1, vars = c("x", "z"), w = 1, noise = "domain",
    rejection = TRUE, output = "tibble", seed = 1
  )
  expect_true(all(abs(z$x_proj) < 1e-8))
  expect_true(all(abs(z$z_proj) <= 1))

  # the codomain rejection sampler too
  withr::local_seed(2)
  expect_equal(nrow(rvnorm(20, p, sd = 0.1, vars = c("x", "z"), w = 1, rejection = TRUE)), 20)
})

test_that("free variables leave a system's least squares alone", {
  # pdvnorm_pseudoinverse() row-scales a system that is not overdetermined.
  # Counting the free z, (x, 2 x - 2) looked square, and its least-squares
  # residual x - 0.8 became x - 0.5
  p <- mp(c("x", "2 x - 2"))
  withr::local_seed(1)
  b <- rvnorm(1000, p, sd = 0.1, w = 1.5, rejection = TRUE, vars = c("x", "z"))
  expect_lt(abs(mean(b$x) - 0.8), 0.015)

  # domain noise on a system whose two components weigh the same. Its
  # least-squares pilot also puts draws at x = 0.552, a least-squares point
  # off the variety, whose projections are dropped
  expect_warning(
    z <- rvnorm(
      2000, mp(c("x^2 - x", "x^3 - x^2")), sd = 0.1, noise = "domain",
      rejection = TRUE, output = "tibble", w = 1.5, seed = 53, vars = c("x", "z")
    ),
    "dropped"
  )
  expect_equal(mean(z$x_proj > 0.5), 0.5, tolerance = 0.06)

  J <- cbind(c(1, 2), 0)
  expect_equal(
    pdvnorm_pseudoinverse_cols(J, 1),
    structure(rbind(as.numeric(pdvnorm_pseudoinverse(J[, 1, drop = FALSE])), 0), rank = 1L)
  )
})

test_that("free variables leave a system's geometry alone with CmdStan", {
  skip_if_no_cmdstan()
  z <- rvnorm(
    200, mp(c("x", "x^2")), sd = 0.1, vars = c("x", "z"), w = 1,
    noise = "domain", output = "tibble", chains = 2, cores = 1, seed = 3
  )
  expect_true(all(abs(z$x_proj) < 1e-8))
  expect_true(all(abs(z$z_proj) <= 1))
})

test_that("rvnorm(noise = \"domain\") corrects for repeated factors", {
  # x^2 - x and x (x - 1)^2 both cut out {0, 1}, but the pilot's tube around
  # the double root is twice as wide, so 1 got two thirds of the points at
  # any pilot_sd
  share <- function(expr) {
    z <- rvnorm(
      2000, mp(expr), sd = 0.1, noise = "domain", rejection = TRUE,
      output = "tibble", w = 1.5, seed = 53
    )
    mean(z$x_proj > 0.5)
  }
  expect_equal(share("x^2 - x"), 0.5, tolerance = 0.06)
  expect_equal(share("x^3 - 2 x^2 + x"), 0.5, tolerance = 0.06)

  # a squared circle is still uniform in angle. Its gradient vanishes on
  # the circle, which rejected the projections that landed on it exactly
  expect_true(project_rows_closest(matrix(c(1, 0), 1), mp("x^2 + y^2 - 1")^2, c("x", "y"))$ok)
  expect_no_warning(z <- rvnorm(
    1000, mp("x^2 + y^2 - 1")^2, sd = 0.1, noise = "domain", rejection = TRUE,
    output = "tibble", w = 1.5, seed = 54
  ))
  angle <- atan2(z$y_proj, z$x_proj) %% (2 * pi) / (2 * pi)
  expect_gt(stats::ks.test(angle, "punif")$p.value, 0.01)

  # on a cube, Gauss-Newton steps of rounding error over a vanishing
  # gradient threw points across the circle
  expect_no_warning(z <- rvnorm(
    1000, mp("x^2 + y^2 - 1")^3, sd = 0.1, noise = "domain", rejection = TRUE,
    output = "tibble", w = 1.5, seed = 55
  ))
  expect_lt(max(abs(sqrt(z$x_proj^2 + z$y_proj^2) - 1)), 1e-4)
  angle <- atan2(z$y_proj, z$x_proj) %% (2 * pi) / (2 * pi)
  expect_gt(stats::ks.test(angle, "punif")$p.value, 0.01)

  # a system's correction cannot be estimated from one draw: it stops
  expect_error(
    rvnorm(
      200, mp(c("x^2 + y^2 - 1", "y^2")), sd = 0.1, noise = "domain",
      rejection = TRUE, w = 1.5, seed = 1
    ),
    "repeated factor"
  )
})

test_that("project_rows_closest() flags rows it cannot project", {
  # the gradient vanishes at the center
  expect_false(project_rows_closest(matrix(0, 1, 2), circle, c("x", "y"))$ok)
  far <- matrix(c(3, 0), 1)
  expect_false(project_rows_closest(far, circle, c("x", "y"), max_dist = 1)$ok)
  expect_true(project_rows_closest(far, circle, c("x", "y"), max_dist = 3)$ok)
})

test_that("rvnorm(noise = \"domain\") projects pilot draws and adds noise", {
  calls <- mock_pilot(near_circle(200))
  S <- matrix(c(0.01, 0.012, 0.012, 0.04), 2)

  out <- rvnorm(
    100, circle, Sigma = S, noise = "domain", output = "tibble",
    chains = 2, seed = 3, sig_figs = 12
  )
  expect_named(out, c("x", "y", "x_proj", "y_proj", ".chain", ".iteration", ".draw"))
  expect_equal(nrow(out), 100)
  u <- as.matrix(out[, c("x_proj", "y_proj")])
  expect_equal(sqrt(rowSums(u^2)), rep(1, 100), tolerance = 1e-10)
  # an evenly spaced subset of the pilot's draws keeps each chain's share
  expect_equal(as.vector(table(out$.chain)), c(50L, 50L))
  expect_equal(out$.draw, 1:100)

  # the pilot: a tenth of the noise scale, the seed, sampler settings
  args <- calls$args[[1]]
  expect_equal(args$pilot_sd, sqrt(0.025) / 10)
  expect_equal(args$seed, 3)
  expect_false(args$rejection)
  expect_equal(args$pilot_args$chains, 2L)
  expect_equal(args$pilot_args$sig_figs, 12)

  # "simple" output holds the same noisy points
  simple <- rvnorm(100, circle, Sigma = S, noise = "domain", chains = 2, seed = 3)
  expect_equal(simple, as.data.frame(out[, c("x", "y")]))
})

test_that("rvnorm(noise = \"domain\") noise has covariance Sigma", {
  mock_pilot(near_circle(6000))
  S <- matrix(c(0.01, 0.012, 0.012, 0.04), 2)
  out <- rvnorm(5000, circle, Sigma = S, noise = "domain", output = "tibble", seed = 1)
  e <- as.matrix(out[, c("x", "y")] - out[, c("x_proj", "y_proj")])
  # within 4 standard errors, entry by entry
  se <- sqrt((outer(diag(S), diag(S)) + S^2) / nrow(e))
  expect_true(all(abs(stats::cov(e) - S) < 4 * se))
  expect_true(all(abs(colMeans(e)) < 4 * sqrt(diag(S) / nrow(e))))
  # Mahalanobis distances are chi-squared with 2 degrees of freedom
  d2 <- rowSums((e %*% solve(S)) * e)
  expect_gt(stats::ks.test(d2, "pchisq", 2)$p.value, 0.001)
})

test_that("rvnorm(noise = \"domain\") widens the pilot's window and keeps points inside w", {
  calls <- mock_pilot(near_circle(400))

  out <- rvnorm(
    100, circle, sd = 0.1, noise = "domain", output = "tibble", chains = 2,
    w = list(x = c(-0.5, 2), y = c(-2, 2))
  )
  # pilot_sd = 0.01, so the pilot's window is wider by 0.06
  expect_equal(calls$args[[1]]$pilot_args$w, list(x = c(-0.56, 2.06), y = c(-2.06, 2.06)))
  expect_true(all(out$x_proj >= -0.5))
  expect_equal(nrow(out), 100)

  # a scalar window stays scalar, so the templates still apply
  rvnorm(100, circle, sd = 0.1, noise = "domain", chains = 2, w = 1.5)
  expect_equal(calls$args[[2]]$pilot_args$w, 1.56)

  # variables the polynomial does not involve keep their bounds
  calls_z <- mock_pilot(near_circle(400, z = c(0, 1)))
  out_z <- rvnorm(
    100, circle, sd = 0.1, noise = "domain", output = "tibble", chains = 2,
    vars = c("x", "y", "z"), w = list(x = c(-2, 2), y = c(-2, 2), z = c(0, 1))
  )
  args_z <- calls_z$args[[1]]$pilot_args
  expect_equal(args_z$w$z, c(0, 1))
  expect_equal(args_z$w$x, c(-2.06, 2.06))
  expect_equal(args_z$vars, c("x", "y", "z"))
  expect_named(out_z, c("x", "y", "z", "x_proj", "y_proj", "z_proj", ".chain", ".iteration", ".draw"))

  # a scalar window too: it becomes per-variable bounds
  rvnorm(
    100, circle, sd = 0.1, noise = "domain", chains = 2,
    vars = c("x", "y", "z"), w = 1.5
  )
  expect_equal(
    calls_z$args[[2]]$pilot_args$w,
    list(x = c(-1.56, 1.56), y = c(-1.56, 1.56), z = c(-1.5, 1.5))
  )
})

test_that("rvnorm(noise = \"domain\") keeps a narrow window on free variables", {
  # widening y and z with x made the pilot draw them from (-0.61, 0.61),
  # and all but about 1 in 3721 latent points left the window
  z <- rvnorm(
    20, "x", sd = 1, noise = "domain", rejection = TRUE,
    vars = c("x", "y", "z"), w = 0.01, seed = 1, output = "tibble"
  )
  expect_equal(nrow(z), 20)
  expect_true(all(abs(z$x_proj) < 1e-8))
  expect_true(all(abs(as.matrix(z[, c("y_proj", "z_proj")])) <= 0.01))
})

test_that("rvnorm(noise = \"domain\") replaces failed projections from new pilot runs", {
  center <- data.frame(x = 0, y = 0, .chain = 1L, .iteration = 0L)
  calls <- mock_pilot(function(call) {
    if (call == 1L) {
      rbind(near_circle(40, seed = 1), center[rep(1, 3), ])
    } else {
      near_circle(100, seed = 2)
    }
  })

  expect_warning(
    out <- rvnorm(80, circle, sd = 0.1, noise = "domain", output = "tibble", seed = 7),
    "3 of 143 pilot draws were dropped: 3 projections did not converge"
  )
  expect_equal(vapply(calls$args, `[[`, numeric(1), "seed"), c(7, 8))
  expect_equal(nrow(out), 80)
  expect_equal(anyDuplicated(as.data.frame(out[, c("x_proj", "y_proj")])), 0L)
  # the second run's chains follow the first run's two, so a chain and
  # iteration identify one pilot draw
  expect_equal(anyDuplicated(as.data.frame(out[, c(".chain", ".iteration")])), 0L)
  expect_setequal(unique(out$.chain), 1:4)

  mock_pilot(center[rep(1, 20), ])
  expect_error(
    rvnorm(10, circle, sd = 0.1, noise = "domain"),
    "Only 0 of 10 pilot draws projected onto the variety"
  )

  # many projections that do not converge stop, as parts of the variety
  # would be missing, instead of a warning over a full sample
  mock_pilot(rbind(near_circle(70, seed = 3), center[rep(1, 30), ]))
  expect_error(
    rvnorm(50, circle, sd = 0.1, noise = "domain"),
    "30 of 100 projections onto the variety did not converge"
  )
})

test_that("rvnorm(noise = \"domain\") corrects multiplicities up to the degree", {
  # a cap of 4 took the draws at x^6 - x^5's root 0, of multiplicity 5, for
  # jumps to another branch, and every point went to 1
  z <- rvnorm(
    2000, mp("x^6 - x^5"), sd = 0.1, pilot_sd = 0.01, noise = "domain",
    rejection = TRUE, output = "tibble", w = 1.5, seed = 53
  )
  expect_equal(mean(z$x_proj > 0.5), 0.5, tolerance = 0.06)
})

test_that("rvnorm(noise = \"domain\") seeds its noise without changing R's stream", {
  mock_pilot(near_circle(200))
  withr::local_seed(10)
  before <- get(".Random.seed", envir = globalenv())
  a <- rvnorm(100, circle, sd = 0.1, noise = "domain", seed = 5)
  expect_identical(get(".Random.seed", envir = globalenv()), before)
  expect_identical(rvnorm(100, circle, sd = 0.1, noise = "domain", seed = 5), a)
  expect_false(identical(rvnorm(100, circle, sd = 0.1, noise = "domain", seed = 6), a))

  # without a seed, set.seed() reproduces the draws
  set.seed(4)
  b <- rvnorm(100, circle, sd = 0.1, noise = "domain")
  set.seed(4)
  expect_identical(rvnorm(100, circle, sd = 0.1, noise = "domain"), b)
})

test_that("rvnorm(noise = \"domain\") rejects options it does not support", {
  expect_error(
    rvnorm(10, circle, sd = 0.1, noise = "domain", homo = FALSE),
    "`homo` describes codomain noise"
  )
  expect_error(
    rvnorm(10, circle, sd = 0.1, noise = "domain", inc_warmup = TRUE),
    "`inc_warmup = TRUE` is not available"
  )
  expect_error(
    rvnorm(10, circle, sd = 0.1, noise = "domain", output = "stanfit"),
    "`output = \"stanfit\"` is not available"
  )
  expect_error(
    rvnorm(10, circle, sd = 0.1, noise = "domain", pilot_sd = -1),
    "`pilot_sd` must be NULL or a positive"
  )
  expect_error(rvnorm(10, circle, sd = 0.1, noise = "range"), "should be one of")
  expect_error(
    rvnorm(
      10, circle, sd = 0.1, noise = "domain", output = "tibble",
      vars = c("x", "x_proj", "y"), w = 1
    ),
    "clash"
  )
})

test_that("rvnorm(noise = \"domain\", code_only = TRUE) returns the pilot's program", {
  code <- rvnorm(10, circle, sd = 0.1, noise = "domain", w = 1.5, code_only = TRUE)
  expect_match(code, "real<lower=-1.56,upper=1.56> x;", fixed = TRUE)
  expect_match(code, "g/ndg", fixed = TRUE)
})

test_that("rvnorm(noise = \"domain\") samples its pilot with the window templates", {
  draws <- near_circle(200)
  capture <- new.env(parent = emptyenv())
  exe <- withr::local_tempfile()
  writeLines("fake", exe)
  fit <- structure(
    list(draws = function(format = "df", inc_warmup = FALSE) draws),
    class = "CmdStanMCMC"
  )
  testthat::local_mocked_bindings(
    cmdstan_model = function(...) stop("cmdstan_model should not be called"),
    .package = "vnorm"
  )
  testthat::local_mocked_bindings(
    stan_package_model = function(name, ...) {
      capture$name <- name
      list(
        sample = function(...) {
          capture$args <- list(...)
          fit
        },
        exe_file = function() exe
      )
    },
    .package = "instantiate"
  )

  out <- rvnorm(100, circle, sd = 0.1, noise = "domain", w = 1.5, chains = 2)
  data <- jsonlite::read_json(capture$args$data, simplifyVector = TRUE)
  expect_identical(capture$name, "2_2_vn_w")
  expect_equal(data$si, 0.01)
  expect_equal(data$w, 1.56)
  expect_equal(nrow(out), 100)
})

test_that("rvnorm(noise = \"domain\", rejection = TRUE) is uniform in arclength", {
  # the rejection pilot's draws are independent, so a KS test applies
  ell <- mp("0.25 x^2 + y^2 - 1")
  out <- withr::with_seed(1, rvnorm(
    1500, ell, sd = 0.05, noise = "domain", rejection = TRUE,
    w = list(x = c(-2.3, 2.3), y = c(-1.3, 1.3)), output = "tibble"
  ))
  # normalized arclength from the parameter angle t of (2 cos t, sin t)
  t_grid <- seq(0, 2 * pi, length.out = 20001)
  speed <- sqrt(4 * sin(t_grid)^2 + cos(t_grid)^2)
  s_grid <- c(0, cumsum((speed[-1] + speed[-length(speed)]) / 2 * diff(t_grid)))
  angle <- atan2(out$y_proj, out$x_proj / 2) %% (2 * pi)
  arclength <- stats::approx(t_grid, s_grid / max(s_grid), xout = angle)$y

  expect_gt(stats::ks.test(arclength, "punif")$p.value, 0.01)
  # points uniform in the angle instead would be detected
  expect_lt(stats::ks.test(angle / (2 * pi), "punif")$p.value, 0.01)

  e <- as.matrix(out[, c("x", "y")] - out[, c("x_proj", "y_proj")])
  expect_gt(stats::ks.test(rowSums(e^2) / 0.05^2, "pchisq", 2)$p.value, 0.01)
})

test_that("rvnorm(noise = \"domain\") samples with CmdStan", {
  skip_if_no_cmdstan()
  S <- matrix(c(0.01, 0.012, 0.012, 0.04), 2)
  out <- rvnorm(
    2000, circle, Sigma = S, noise = "domain", output = "tibble",
    chains = 2, cores = 1, seed = 2
  )
  u <- as.matrix(out[, c("x_proj", "y_proj")])
  expect_equal(sqrt(rowSums(u^2)), rep(1, 2000), tolerance = 1e-10)
  expect_setequal(unique(out$.chain), 1:2)
  e <- as.matrix(out[, c("x", "y")]) - u
  se <- sqrt((outer(diag(S), diag(S)) + S^2) / nrow(e))
  expect_true(all(abs(stats::cov(e) - S) < 4 * se))
})
