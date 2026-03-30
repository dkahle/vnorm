test_that("make_stan_files writes expected univariate Stan template", {
  path <- file.path("src", "stan", "2_2_vn_w.stan")
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  old <- if (file.exists(path)) readLines(path, warn = FALSE) else NULL
  withr::defer({
    if (is.null(old)) {
      if (file.exists(path)) unlink(path)
    } else {
      writeLines(old, con = path)
    }
  })

  expect_no_error(make_stan_files(num_of_vars = 2, totaldeg = 2, homo = TRUE, windowed = TRUE))
  expect_true(file.exists(path))
  txt <- paste(readLines(path, warn = FALSE), collapse = "\n")
  expect_match(txt, "data \\{")
  expect_match(txt, "real si;")
  expect_match(txt, "real<lower=-w, upper=w> x;")
  expect_match(txt, "normal_lpdf\\(0\\.00 \\| g/ndg, si\\)")
})

test_that("make_stan_files supports heteroskedastic and no-window template", {
  path <- file.path("src", "stan", "1_3_hvn.stan")
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  old <- if (file.exists(path)) readLines(path, warn = FALSE) else NULL
  withr::defer({
    if (is.null(old)) {
      if (file.exists(path)) unlink(path)
    } else {
      writeLines(old, con = path)
    }
  })

  make_stan_files(num_of_vars = 1, totaldeg = 3, homo = FALSE, windowed = FALSE)
  txt <- paste(readLines(path, warn = FALSE), collapse = "\n")
  expect_false(grepl("real w;", txt, fixed = TRUE))
  expect_match(txt, "real ndg = 1;")
})
