test_that("poly_to_df evaluates on a regular grid with shift", {
  p <- mp("x + y")
  out <- poly_to_df(p, xlim = c(0, 1), ylim = c(0, 1), nx = 3, ny = 2, shift = 1)

  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), 6)
  expect_true(all(c("x", "y", "z") %in% names(out)))
  expect_equal(out$z[1], 1)
})

test_that("geom_variety sign helpers ignore near-zero noise but detect real sign changes", {
  z_noise <- c(1e-20, 0, -1e-20, 2e-20)
  z_change <- c(-0.5, 0, 0.25)

  expect_false(has_strict_sign_change(z_noise))
  expect_true(has_strict_sign_change(z_change))
})

test_that("split_large_projected_jumps splits path groups at large jumps", {
  df <- data.frame(
    x = c(0, 0.1, 0.2, 5, 5.1),
    y = c(0, 0.1, 0.2, 5, 5.1),
    group = factor(rep("1-001", 5)),
    piece = rep(1L, 5)
  )

  out <- split_large_projected_jumps(df)

  expect_s3_class(out, "data.frame")
  expect_gt(length(unique(out$group)), 1)
  expect_equal(nrow(out), nrow(df))
})

test_that("effective_variety_grid boosts hard shifted no-sign-change cases", {
  ordinary <- effective_variety_grid(201, 201, shift = 0, no_sign_change0 = FALSE)
  shifted <- effective_variety_grid(201, 201, shift = -0.001, no_sign_change0 = FALSE)
  hard <- effective_variety_grid(201, 201, shift = -0.001, no_sign_change0 = TRUE)

  expect_equal(ordinary$nx, 201L)
  expect_equal(ordinary$ny, 201L)
  expect_equal(shifted$nx, 301L)
  expect_equal(shifted$ny, 301L)
  expect_equal(hard$nx, 401L)
  expect_equal(hard$ny, 401L)
})

test_that("check_sign_warning emits message for all-positive with shift=0", {
  expect_message(
    vnorm:::check_sign_warning(c(1, 2, 3, 4), shift = 0),
    "All values are positive"
  )
})

test_that("check_sign_warning emits message for all-negative with shift=0", {
  expect_message(
    vnorm:::check_sign_warning(c(-1, -2, -3), shift = 0),
    "All values are negative"
  )
})

test_that("check_sign_warning emits message for positive after nonzero shift", {
  expect_message(
    vnorm:::check_sign_warning(c(0.5, 1, 2), shift = -0.1),
    "positive after applying shift"
  )
})

test_that("snap_shifted_contours_to_variety reduces polynomial residuals", {
  p <- mp("x^2 + y^2 - 1")
  theta <- seq(0, 2 * pi, length.out = 50)
  df <- data.frame(
    x = 1.05 * cos(theta), y = 1.05 * sin(theta),
    group = factor(1), piece = 1L
  )
  out <- vnorm:::snap_shifted_contours_to_variety(df, p)
  pf <- as.function(p, varorder = c("x", "y"), silent = TRUE)
  residuals <- apply(as.matrix(out[, c("x", "y")]), 1, pf)
  expect_true(max(abs(residuals)) < max(abs(apply(as.matrix(df[, c("x", "y")]), 1, pf))))
})

test_that("collapse_near_duplicate_contours merges near-coincident paths", {
  theta <- seq(0, 2 * pi, length.out = 30)
  g1 <- data.frame(x = cos(theta), y = sin(theta), group = factor("a"), piece = 1L)
  g2 <- data.frame(x = cos(theta) + 1e-5, y = sin(theta) + 1e-5, group = factor("b"), piece = 2L)
  df <- rbind(g1, g2)
  out <- vnorm:::collapse_near_duplicate_contours(df, tol = 0.01)
  expect_equal(length(unique(out$group)), 1)
})

test_that("is_fragmented_paths detects fragmented vs non-fragmented", {
  df_ok <- data.frame(x = 1:100, y = 1:100, group = factor(rep(1:2, each = 50)))
  expect_false(vnorm:::is_fragmented_paths(df_ok))

  df_frag <- data.frame(x = 1:100, y = 1:100, group = factor(rep(1:20, each = 5)))
  expect_true(vnorm:::is_fragmented_paths(df_frag))
})
