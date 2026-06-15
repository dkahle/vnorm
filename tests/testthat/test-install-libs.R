test_that("R CMD INSTALL succeeds with an unwritable user cache", {
  skip_if_not(
    identical(Sys.getenv("VNORM_RUN_INSTALL_TEST"), "true"),
    "Set VNORM_RUN_INSTALL_TEST=true to run install integration test."
  )

  pkg_root <- normalizePath(file.path(test_path(), "..", ".."), mustWork = TRUE)
  lib <- tempfile("vnorm-install-lib-")
  cache_parent <- tempfile("vnorm-cache-parent-")
  dir.create(lib)
  dir.create(cache_parent)
  unwritable_cache <- file.path(cache_parent, "cache")
  dir.create(unwritable_cache)
  Sys.chmod(unwritable_cache, mode = "0555")
  withr::defer(Sys.chmod(unwritable_cache, mode = "0755"))
  withr::defer(unlink(c(lib, cache_parent), recursive = TRUE, force = TRUE))

  env <- c(
    paste0("XDG_CACHE_HOME=", unwritable_cache),
    "VNORM_PRECOMPILE_STAN=false"
  )
  result <- system2(
    file.path(R.home("bin"), "R"),
    c("CMD", "INSTALL", paste0("--library=", lib), pkg_root),
    stdout = TRUE,
    stderr = TRUE,
    env = env
  )

  status <- attr(result, "status")
  if (is.null(status)) status <- 0L
  expect_equal(status, 0L)
})

test_that("Stan precompilation defaults to auto during installation", {
  install_script_path <- test_path("../../src/install.libs.R")
  skip_if_not(
    file.exists(install_script_path),
    "src/install.libs.R is only available in source-tree tests."
  )

  install_script <- readLines(install_script_path, warn = FALSE)
  expect_true(any(grepl('Sys.getenv\\("VNORM_PRECOMPILE_STAN", "auto"\\)', install_script)))
})
