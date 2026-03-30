test_that("compiled Stan cache deduplicates entries", {
  on.exit(vnorm:::clear_compiled_stan_info(), add = TRUE)
  vnorm:::clear_compiled_stan_info()

  vnorm:::add_compiled_stan_info("m1", "/tmp/model1.stan")
  vnorm:::add_compiled_stan_info("m1", "/tmp/model1.stan")
  info <- vnorm:::get_compiled_stan_info()

  expect_true(is.data.frame(info))
  expect_equal(nrow(info), 1)
  expect_equal(info$name[[1]], "m1")
  expect_equal(info$path[[1]], "/tmp/model1.stan")
})

test_that("remove_stan_files deletes tracked files and clears cache", {
  on.exit(vnorm:::clear_compiled_stan_info(), add = TRUE)
  vnorm:::clear_compiled_stan_info()

  td <- tempfile("vnorm-remove-test-")
  dir.create(td)
  stan_path <- file.path(td, "tmp_model.stan")
  exe_path <- sub("\\.stan$", "", stan_path)
  file.create(stan_path)
  file.create(exe_path)

  vnorm:::set_compiled_stan_info(
    data.frame(name = "tmp_model", path = stan_path, stringsAsFactors = FALSE)
  )

  expect_true(file.exists(stan_path))
  expect_true(file.exists(exe_path))

  expect_message(remove_stan_files(), "Files deleted")
  expect_false(file.exists(stan_path))
  expect_false(file.exists(exe_path))
  expect_equal(nrow(vnorm:::get_compiled_stan_info()), 0)
})

test_that("rvnorm(user_compiled=TRUE) errors when cache is empty", {
  on.exit(vnorm:::clear_compiled_stan_info(), add = TRUE)
  vnorm:::clear_compiled_stan_info()

  expect_error(
    rvnorm(
      n = 10,
      poly = mp("x^2 + y^2 - 1"),
      sd = 0.1,
      user_compiled = TRUE
    ),
    "No compiled model cache found"
  )
})

test_that("rvnorm(user_compiled=TRUE) errors when model is missing in cache", {
  on.exit(vnorm:::clear_compiled_stan_info(), add = TRUE)
  vnorm:::clear_compiled_stan_info()

  vnorm:::set_compiled_stan_info(
    data.frame(
      name = "other_model",
      path = "/tmp/other.stan",
      stringsAsFactors = FALSE
    )
  )

  expect_error(
    rvnorm(
      n = 10,
      poly = mp("x^2 + y^2 - 1"),
      sd = 0.1,
      user_compiled = TRUE
    ),
    "Requested compiled model not found in cache"
  )
})

test_that("multiple distinct models stored in cache", {
  on.exit(vnorm:::clear_compiled_stan_info(), add = TRUE)
  vnorm:::clear_compiled_stan_info()
  vnorm:::add_compiled_stan_info("model_a", "/tmp/a.stan")
  vnorm:::add_compiled_stan_info("model_b", "/tmp/b.stan")
  info <- vnorm:::get_compiled_stan_info()
  expect_equal(nrow(info), 2)
  expect_true(all(c("model_a", "model_b") %in% info$name))
})

test_that("clear_compiled_stan_info on already-empty cache does not error", {
  vnorm:::clear_compiled_stan_info()
  expect_no_error(vnorm:::clear_compiled_stan_info())
})

test_that("get_compiled_stan_info returns zero-row data.frame on fresh state", {
  vnorm:::clear_compiled_stan_info()
  info <- vnorm:::get_compiled_stan_info()
  expect_true(is.data.frame(info))
  expect_equal(nrow(info), 0)
  expect_true(all(c("name", "path") %in% names(info)))
})
