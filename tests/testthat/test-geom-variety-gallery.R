test_that("geom_variety gallery cases build for projection off and auto", {
  cases <- list(
    list(name = "circle", poly = mp("x^2 + y^2 - 1"), xlim = c(-2, 2), ylim = c(-2, 2), shift = 0),
    list(name = "line", poly = mp("y - x"), xlim = c(-2, 2), ylim = c(-2, 2), shift = 0),
    list(name = "heart", poly = mp("(x^2 + y^2 - 1)^3 - x^2 y^3"), xlim = c(-2, 2), ylim = c(-2, 2), shift = 0),
    list(name = "folium", poly = mp("x^3 + y^3 - 3*x*y"), xlim = c(-2, 3), ylim = c(-2, 3), shift = 0),
    list(name = "squared_circle", poly = mp("x^2 + y^2 - 1")^2, xlim = c(-2, 2), ylim = c(-2, 2), shift = -0.000576),
    list(name = "squared_line", poly = mp("y - x")^2, xlim = c(-2, 2), ylim = c(-2, 2), shift = -0.001936),
    list(name = "squared_crossing", poly = mp("y^2 - x^2")^2, xlim = c(-2, 2), ylim = c(-2, 2), shift = -0.004),
    list(name = "squared_heart", poly = mp("((x^2 + y^2 - 1)^3 - x^2 y^3)")^2, xlim = c(-2, 2), ylim = c(-2, 2), shift = -0.0001),
    list(name = "squared_folium", poly = mp("(x^3 + y^3 - 3*x*y)")^2, xlim = c(-2, 3), ylim = c(-2, 3), shift = -0.001)
  )

  for (case in cases) {
    off_plot <- ggplot() +
      geom_variety(
        poly = case$poly,
        xlim = case$xlim,
        ylim = case$ylim,
        shift = case$shift,
        projection = "off"
      )
    auto_plot <- ggplot() +
      geom_variety(
        poly = case$poly,
        xlim = case$xlim,
        ylim = case$ylim,
        shift = case$shift,
        projection = "auto"
      )

    off_dat <- ggplot2::ggplot_build(off_plot)$data[[1]]
    auto_dat <- ggplot2::ggplot_build(auto_plot)$data[[1]]

    expect_true(nrow(off_dat) > 0, info = paste(case$name, "projection off"))
    expect_true(nrow(auto_dat) > 0, info = paste(case$name, "projection auto"))
  }
})
