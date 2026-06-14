test_that(".onLoad sets cmdstanr version-check option", {
  old <- getOption("cmdstanr_no_ver_check")
  withr::defer(options(cmdstanr_no_ver_check = old))

  options(cmdstanr_no_ver_check = NULL)
  get(".onLoad", asNamespace("vnorm"), inherits = FALSE)("vnorm", "vnorm")

  expect_true(isTRUE(getOption("cmdstanr_no_ver_check")))
})
