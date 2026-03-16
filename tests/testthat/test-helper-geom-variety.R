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
