test_that("geom_variety works with mpoly objects", {
  poly <- mp("y - x")
  p <- ggplot() + geom_variety(poly = poly, xlim = c(-2, 2), ylim = c(-2, 2))
  expect_true(inherits(p, "ggplot"))
})

test_that("geom_variety works with mpolyList objects", {
  poly1 <- mp("x^2 + y^2 - 1")
  poly2 <- mp("y - x")
  poly_list <- mpolyList(poly1, poly2)
  p <- ggplot() +
    geom_variety(poly = poly_list, xlim = c(-2, 2), ylim = c(-2, 2))
  expect_true(inherits(p, "ggplot"))
})

test_that("geom_variety maps linetype by polynomial for mpolyList", {
  poly_list <- mpolyList(mp("x^2 + y^2 - 1"), mp("y - x"))
  p <- ggplot() +
    geom_variety(poly = poly_list, xlim = c(-2, 2), ylim = c(-2, 2))

  built <- ggplot2::ggplot_build(p)
  expect_gt(length(unique(built$data[[1]]$linetype)), 1)
})

test_that("geom_variety keeps a single default colour for mpolyList", {
  poly_list <- mpolyList(mp("x^2 + y^2 - 1"), mp("y - x"))
  p <- ggplot() +
    geom_variety(poly = poly_list, xlim = c(-2, 2), ylim = c(-2, 2))

  built <- ggplot2::ggplot_build(p)
  expect_equal(length(unique(built$data[[1]]$colour)), 1)
})

test_that("geom_variety supports scale_colour_manual with vary_colour", {
  poly_list <- mpolyList(mp("x^2 + y^2 - 1"), mp("y - x"))
  p <- NULL
  expect_no_warning(
    p <- ggplot() +
      geom_variety(
        poly = poly_list,
        xlim = c(-2, 2),
        ylim = c(-2, 2),
        vary_colour = TRUE
      ) +
      ggplot2::scale_colour_manual(values = c("steelblue", "tomato"))
  )

  built <- ggplot2::ggplot_build(p)
  expect_gt(length(unique(built$data[[1]]$colour)), 1)
})

test_that("geom_variety handles xlim and ylim parameters", {
  poly <- mp("x^2 + y^2 - 1")
  p <- ggplot() + geom_variety(poly = poly, xlim = c(-1, 1), ylim = c(-1, 1))
  expect_equal(p$layers[[1]]$stat_params$xlim, c(-1, 1))
  expect_equal(p$layers[[1]]$stat_params$ylim, c(-1, 1))
})

test_that("geom_variety computes correct paths for simple varieties", {
  poly <- mp("y - x")
  p <- ggplot() + geom_variety(poly = poly, xlim = c(-2, 2), ylim = c(-2, 2))
  ggplot2::ggplot_build(p)
  df <- poly_to_df(poly, xlim = c(-2, 2), ylim = c(-2, 2), nx = 100, ny = 100)
  expect_true("x" %in% names(df) && "y" %in% names(df) && "z" %in% names(df))
})

test_that("geom_variety generates valid non-empty contour data", {
  poly <- mp("x^2 + y^2 - 1")
  plot_valid <- ggplot() +
    geom_variety(poly = poly, xlim = c(-2, 2), ylim = c(-2, 2)) +
    ggplot2::coord_equal() +
    labs(title = "Circle Variety Plot")

  built <- ggplot2::ggplot_build(plot_valid)
  dat <- built$data[[1]]
  expect_gt(nrow(dat), 0)
  expect_true(all(c("x", "y", "group") %in% names(dat)))
})

test_that("geom_variety projection off returns raw shifted contours", {
  poly <- mp("x^2 + y^2 - 1")^2
  dat <- ggplot2::ggplot_build(
    ggplot() +
      geom_variety(
        poly = poly,
        xlim = c(-2, 2),
        ylim = c(-2, 2),
        shift = -0.000576,
        projection = "off"
      )
  )$data[[1]]

  expect_gt(nrow(dat), 0)
  expect_gte(length(unique(dat$group)), 2)
  radii <- sqrt(dat$x^2 + dat$y^2)
  expect_true(min(radii, na.rm = TRUE) < 1)
  expect_true(max(radii, na.rm = TRUE) > 1)
})

test_that("geom_variety default resolution keeps heart contour reasonably smooth", {
  poly <- mp("(x^2 + y^2 - 1)^3 - x^2 y^3")
  dat <- ggplot2::ggplot_build(
    ggplot() + geom_variety(poly = poly, xlim = c(-2, 2), ylim = c(-2, 2))
  )$data[[1]]

  dx <- diff(dat$x)
  dy <- diff(dat$y)
  seg <- sqrt(dx * dx + dy * dy)
  seg <- seg[is.finite(seg) & seg > 0]

  expect_gt(nrow(dat), 400)
  expect_lt(max(seg), 0.035)
})

test_that("geom_variety projects ordinary contours back onto the zero set", {
  poly <- mp("(x^2 + y^2 - 1)^3 - x^2 y^3")
  pf <- as.function(poly, varorder = c("x", "y"), silent = TRUE)
  dat <- ggplot2::ggplot_build(
    ggplot() + geom_variety(poly = poly, xlim = c(-2, 2), ylim = c(-2, 2))
  )$data[[1]]

  resid <- abs(pf(as.matrix(dat[, c("x", "y")])))
  expect_lt(max(resid, na.rm = TRUE), 1e-5)
})

test_that("geom_variety projection off leaves shifted repeated crossing away from origin", {
  dat <- ggplot2::ggplot_build(
    ggplot() +
      geom_variety(
        poly = mp("y^2 - x^2")^2,
        xlim = c(-2, 2),
        ylim = c(-2, 2),
        shift = -0.004,
        projection = "off"
      )
  )$data[[1]]

  r <- sqrt(dat$x^2 + dat$y^2)
  expect_gt(min(r, na.rm = TRUE), 0.05)
})

test_that("projection improves shifted repeated crossing compared with raw contour", {
  p <- mp("y^2 - x^2")^2
  off_dat <- ggplot2::ggplot_build(
    ggplot() +
      geom_variety(
        poly = p,
        xlim = c(-2, 2),
        ylim = c(-2, 2),
        shift = -0.004,
        projection = "off"
      )
  )$data[[1]]
  auto_dat <- ggplot2::ggplot_build(
    ggplot() +
      geom_variety(
        poly = p,
        xlim = c(-2, 2),
        ylim = c(-2, 2),
        shift = -0.004,
        projection = "auto"
      )
  )$data[[1]]

  off_r <- sqrt(off_dat$x^2 + off_dat$y^2)
  auto_r <- sqrt(auto_dat$x^2 + auto_dat$y^2)
  expect_lt(min(auto_r, na.rm = TRUE), min(off_r, na.rm = TRUE))
})

test_that("projection collapses shifted repeated-line duplicates to one visible path", {
  dat <- ggplot2::ggplot_build(
    ggplot() +
      geom_variety(
        poly = mp("y - x")^2,
        xlim = c(-2, 2),
        ylim = c(-2, 2),
        shift = -0.001936,
        projection = "auto"
      )
  )$data[[1]]

  expect_lte(length(unique(dat$group)), 1)
  expect_lt(diff(range(dat$y - dat$x, na.rm = TRUE)), 0.02)
})

test_that("shifted repeated crossings close back to the singular point", {
  dat <- ggplot2::ggplot_build(
    ggplot() +
      geom_variety(
        poly = mp("y^2 - x^2")^2,
        xlim = c(-2, 2),
        ylim = c(-2, 2),
        shift = -0.004,
        projection = "auto"
      )
  )$data[[1]]

  r <- sqrt(dat$x^2 + dat$y^2)
  expect_lt(min(r, na.rm = TRUE), 1e-6)
})

test_that(
  "geom_variety shift generates continuous contours for squared polynomials",
  {
  poly <- mp("x^2 + y^2 - 1")^2
  p <- ggplot() +
    geom_variety(
      poly = poly,
      xlim = c(-2, 2),
      ylim = c(-2, 2),
      shift = -0.000576
    )

  built <- ggplot2::ggplot_build(p)
  dat <- built$data[[1]]

  expect_gt(nrow(dat), 0)

  points_per_group <- table(dat$group)
  expect_lte(length(points_per_group), 2)
  expect_true(all(points_per_group >= 20))
  }
)

test_that("fragmented shifted contours trigger refinement", {
  poly <- mp("x^2 + y^2 - 1")^2
  df <- variety_paths_with_refinement(
    poly = poly,
    rangex = c(-2, 2),
    rangey = c(-2, 2),
    nx = 101,
    ny = 101,
    shift = -0.000576,
    group = 1
  )

  expect_gt(nrow(df), 0)
  expect_false(is_fragmented_paths(df))
})

test_that("shift suggestion is emitted as message when shift is zero", {
  poly <- mp("(x^2 + y^2 - 1)^2")
  p <- ggplot() + geom_variety(poly = poly, xlim = c(-2, 2), ylim = c(-2, 2))
  expect_message(ggplot2::ggplot_build(p), "try shift =")
})

test_that("shifted repeated-factor recovery prints a caution message", {
  poly <- mp("((x^2 + y^2 - 1)^3 - x^2 y^3)^2")
  p <- ggplot() +
    geom_variety(
      poly = poly,
      xlim = c(-2, 2),
      ylim = c(-2, 2),
      shift = -1e-4,
      projection = "auto"
    )

  expect_message(
    ggplot2::ggplot_build(p),
    "projected result may still miss branches"
  )
})

test_that("no-shift repeated crossing polynomial does not emit bogus point contour", {
  poly <- mp("y^2 - x^2")^2
  p <- ggplot() + geom_variety(poly = poly, xlim = c(-2, 2), ylim = c(-2, 2))
  built <- NULL
  expect_message(built <- ggplot2::ggplot_build(p), "Zero contours were generated")
  dat <- built$data[[1]]
  expect_equal(nrow(dat), 0)
})

test_that("shift suggestion for repeated crossing polynomial is nonzero", {
  poly <- mp("y^2 - x^2")^2
  p <- ggplot() + geom_variety(poly = poly, xlim = c(-2, 2), ylim = c(-2, 2))
  msgs <- testthat::capture_messages(ggplot2::ggplot_build(p))
  suggestion <- msgs[grepl("try shift =", msgs, fixed = TRUE)]
  expect_true(length(suggestion) >= 1)
  expect_false(any(grepl("try shift = 0\\.?$", suggestion)))
})
