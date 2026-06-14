make_vnorm_logo_magick <- function(
    preview_path,
    logo_path,
    bg,
    hex_border,
    text_color,
    label = "vnorm"
  ) {
  magick <- Sys.which("magick")
  if (!nzchar(magick)) {
    return(FALSE)
  }

  tmp <- tempfile(pattern = "vnorm-logo-", fileext = rep(".png", 5))
  on.exit(unlink(tmp), add = TRUE)

  curve_path <- tmp[[1]]
  curve_small_path <- tmp[[2]]
  hex_mask_path <- tmp[[3]]
  curve_layer_path <- tmp[[4]]
  curve_clipped_path <- tmp[[5]]

  polygon <- "polygon 518,30 1007,314 1007,886 518,1170 30,886 30,314"
  font <- Sys.getenv("VNORM_LOGO_FONT", unset = "/System/Library/Fonts/HelveticaNeue.ttc")
  font_args <- if (nzchar(font) && file.exists(font)) c("-font", font) else character()

  run_magick <- function(args) {
    status <- system2(magick, args = shQuote(args))
    if (!identical(status, 0L)) {
      stop("ImageMagick logo assembly failed.", call. = FALSE)
    }
  }

  run_magick(c(preview_path, "-fuzz", "2%", "-trim", "+repage", curve_path))
  run_magick(c(curve_path, "-resize", "850x", curve_small_path))
  run_magick(c(
    "-size", "1037x1200", "xc:none",
    "-fill", "white", "-stroke", "none", "-draw", polygon,
    hex_mask_path
  ))
  run_magick(c(
    "-size", "1037x1200", "xc:none",
    curve_small_path, "-geometry", "+94+570", "-composite",
    curve_layer_path
  ))
  run_magick(c(
    curve_layer_path, hex_mask_path,
    "-compose", "DstIn", "-composite",
    curve_clipped_path
  ))
  run_magick(c(
    "-size", "1037x1200", "xc:none",
    "-fill", bg, "-stroke", "none", "-draw", polygon,
    curve_clipped_path, "-composite",
    "-fill", "none", "-stroke", hex_border, "-strokewidth", "26",
    "-draw", polygon,
    "-fill", text_color, "-stroke", "none",
    font_args,
    "-pointsize", "100", "-gravity", "north", "-annotate", "+0+245", label,
    "-depth", "8",
    logo_path
  ))

  TRUE
}

make_vnorm_torus_layer <- function(
    filename,
    seed = 20260531,
    width = 2200,
    height = 900,
    dpi = 300
  ) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required.", call. = FALSE)
  }

  old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  } else {
    NULL
  }
  on.exit({
    if (is.null(old_seed)) {
      if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    } else {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(seed)

  project_torus <- function(theta, phi, r_noise = 0, xyz_noise = 0) {
    major_r <- 1
    minor_r <- 0.32 + r_noise
    x <- (major_r + minor_r * cos(phi)) * cos(theta)
    y <- (major_r + minor_r * cos(phi)) * sin(theta)
    z <- minor_r * sin(phi)
    if (xyz_noise > 0) {
      x <- x + stats::rnorm(length(x), 0, xyz_noise)
      y <- y + stats::rnorm(length(y), 0, xyz_noise)
      z <- z + stats::rnorm(length(z), 0, xyz_noise)
    }
    data.frame(
      x = 1.08 * x,
      y = 0.42 * y + 1.14 * z,
      depth = 0.72 * y - 0.24 * z
    )
  }

  rescale01 <- function(x) {
    rng <- range(x, finite = TRUE)
    if (diff(rng) == 0) return(rep(0.5, length(x)))
    (x - rng[[1]]) / diff(rng)
  }

  mesh_theta <- seq(0, 2 * pi, length.out = 300)
  mesh_phi <- seq(0, 2 * pi, length.out = 20)[-20]
  mesh_a <- do.call(rbind, lapply(seq_along(mesh_phi), function(i) {
    out <- project_torus(mesh_theta, rep(mesh_phi[[i]], length(mesh_theta)))
    out$group <- paste0("phi-", i)
    out
  }))

  tube_phi <- seq(0, 2 * pi, length.out = 120)
  tube_theta <- seq(0, 2 * pi, length.out = 38)[-38]
  mesh_b <- do.call(rbind, lapply(seq_along(tube_theta), function(i) {
    out <- project_torus(rep(tube_theta[[i]], length(tube_phi)), tube_phi)
    out$group <- paste0("theta-", i)
    out
  }))

  mesh <- rbind(mesh_a, mesh_b)
  mesh$tint <- rescale01(mesh$x)

  n <- 7800
  pts <- project_torus(
    theta = stats::runif(n, 0, 2 * pi),
    phi = stats::runif(n, 0, 2 * pi),
    r_noise = stats::rnorm(n, 0, 0.055),
    xyz_noise = 0.022
  )
  pts <- pts[order(pts$depth), ]
  pts$tint <- rescale01(pts$x)
  pts$alpha <- 0.14 + 0.24 * rescale01(pts$depth)

  theta <- seq(0, 2 * pi, length.out = 600)
  top_rim <- project_torus(theta, rep(pi / 2, length(theta)))
  bottom_rim <- project_torus(theta, rep(3 * pi / 2, length(theta)))
  outer_rim <- project_torus(theta, rep(0, length(theta)))

  p <- ggplot2::ggplot() +
    ggplot2::geom_point(
      data = pts,
      ggplot2::aes(x, y, colour = tint, alpha = alpha),
      size = 1.8,
      stroke = 0,
      show.legend = FALSE
    ) +
    ggplot2::geom_path(
      data = mesh,
      ggplot2::aes(x, y, group = group, colour = tint),
      linewidth = 0.22,
      alpha = 0.34,
      lineend = "round",
      show.legend = FALSE
    ) +
    ggplot2::geom_point(
      data = pts,
      ggplot2::aes(x, y, colour = tint, alpha = alpha),
      size = 0.45,
      stroke = 0,
      show.legend = FALSE
    ) +
    ggplot2::geom_path(
      data = outer_rim,
      ggplot2::aes(x, y),
      colour = "#FFE7A3",
      linewidth = 1.15,
      alpha = 0.90,
      lineend = "round"
    ) +
    ggplot2::geom_path(
      data = top_rim,
      ggplot2::aes(x, y),
      colour = "#FF8A1F",
      linewidth = 0.90,
      alpha = 0.78,
      lineend = "round"
    ) +
    ggplot2::geom_path(
      data = bottom_rim,
      ggplot2::aes(x, y),
      colour = "#F037A5",
      linewidth = 0.90,
      alpha = 0.78,
      lineend = "round"
    ) +
    ggplot2::scale_colour_gradientn(
      colours = c("#FF9A1F", "#FF2D7A", "#8B5CF6"),
      limits = c(0, 1)
    ) +
    ggplot2::scale_alpha_identity() +
    ggplot2::coord_equal(xlim = c(-1.7, 1.7), ylim = c(-1.02, 1.00), expand = FALSE) +
    ggplot2::theme_void() +
    ggplot2::theme(
      plot.background = ggplot2::element_rect(fill = "transparent", colour = NA),
      panel.background = ggplot2::element_rect(fill = "transparent", colour = NA),
      plot.margin = ggplot2::margin(0, 0, 0, 0)
    )

  ggplot2::ggsave(
    filename = filename,
    plot = p,
    width = width / dpi,
    height = height / dpi,
    dpi = dpi,
    bg = "transparent",
    device = if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_png else "png"
  )

  invisible(filename)
}

make_vnorm_formula_layer <- function(filename) {
  grDevices::png(filename, width = 1000, height = 180, bg = "transparent", res = 200)
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::par(mar = c(0, 0, 0, 0), xaxs = "i", yaxs = "i")
  graphics::plot.new()
  graphics::plot.window(xlim = c(0, 1), ylim = c(0, 1))
  graphics::text(
    0.5, 0.50,
    labels = expression(italic(x) ~ plain("~") ~ italic(N)[p](italic(g), Sigma)),
    col = "#FFF4C7",
    cex = 3.05,
    family = "serif"
  )
  invisible(filename)
}

make_vnorm_logo_neon <- function(
    out_dir = "man/figures",
    seed = 20260531,
    show_plot = TRUE
  ) {
  magick <- Sys.which("magick")
  if (!nzchar(magick)) {
    stop("ImageMagick is required for the neon logo assembly.", call. = FALSE)
  }
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  tmp <- tempfile(pattern = "vnorm-neon-", fileext = rep(".png", 14))
  on.exit(unlink(tmp), add = TRUE)
  names(tmp) <- c(
    "torus", "torus_small", "torus_glow", "formula", "formula_small", "formula_glow",
    "hex_mask", "gradient", "hex_fill", "glow", "title_mask",
    "title_gradient", "title_fill", "title_glow"
  )

  run_magick <- function(args) {
    status <- system2(magick, args = shQuote(args))
    if (!identical(status, 0L)) {
      stop("ImageMagick logo assembly failed.", call. = FALSE)
    }
  }

  polygon <- "polygon 518,30 1007,314 1007,886 518,1170 30,886 30,314"
  font <- Sys.getenv("VNORM_LOGO_FONT", unset = "/System/Library/Fonts/HelveticaNeue.ttc")
  font_args <- if (nzchar(font) && file.exists(font)) c("-font", font) else character()

  make_vnorm_torus_layer(tmp[["torus"]], seed = seed)
  make_vnorm_formula_layer(tmp[["formula"]])

  run_magick(c(tmp[["torus"]], "-resize", "1030x", "-depth", "8", tmp[["torus_small"]]))
  run_magick(c(tmp[["torus_small"]], "-blur", "0x8", tmp[["torus_glow"]]))
  run_magick(c(tmp[["formula"]], "-resize", "600x", "-depth", "8", tmp[["formula_small"]]))
  run_magick(c(tmp[["formula_small"]], "-blur", "0x6", tmp[["formula_glow"]]))
  run_magick(c(
    "-size", "1037x1200", "xc:none",
    "-fill", "white", "-stroke", "none", "-draw", polygon,
    tmp[["hex_mask"]]
  ))
  run_magick(c(
    "-size", "1037x1200", "gradient:#FF9A1F-#B026D9",
    "-rotate", "90", "-resize", "1037x1200!",
    tmp[["gradient"]]
  ))
  run_magick(c(
    tmp[["gradient"]], tmp[["hex_mask"]],
    "-compose", "DstIn", "-composite",
    tmp[["hex_fill"]]
  ))
  run_magick(c(
    "-size", "1037x1200", "xc:none",
    "-fill", "none", "-stroke", "#F037A5", "-strokewidth", "56",
    "-draw", polygon, "-blur", "0x24",
    tmp[["glow"]]
  ))
  run_magick(c(
    "-size", "1037x1200", "xc:none",
    "-fill", "white", font_args,
    "-pointsize", "142", "-gravity", "north", "-annotate", "+0+198", "vnorm",
    tmp[["title_mask"]]
  ))
  run_magick(c(
    "-size", "1037x1200", "gradient:#FF9A1F-#FF2D7A",
    "-rotate", "90", "-resize", "1037x1200!",
    tmp[["title_gradient"]]
  ))
  run_magick(c(
    tmp[["title_gradient"]], tmp[["title_mask"]],
    "-compose", "CopyOpacity", "-composite",
    tmp[["title_fill"]]
  ))
  run_magick(c(tmp[["title_fill"]], "-blur", "0x9", tmp[["title_glow"]]))

  logo_path <- file.path(out_dir, "logo.png")
  run_magick(c(
    "-size", "1037x1200", "xc:none",
    tmp[["glow"]], "-composite",
    tmp[["hex_fill"]], "-composite",
    "-fill", "#00000088", "-stroke", "none", "-draw", polygon,
    tmp[["torus_glow"]], "-geometry", "+4+525", "-composite",
    tmp[["torus_small"]], "-geometry", "+4+525", "-composite",
    tmp[["formula_glow"]], "-geometry", "+218+360", "-composite",
    tmp[["formula_small"]], "-geometry", "+218+360", "-composite",
    tmp[["title_glow"]], "-composite",
    tmp[["title_fill"]], "-composite",
    "-fill", "none", "-stroke", "#050505", "-strokewidth", "20",
    "-draw", polygon,
    "-stroke", "#FF3D9A", "-strokewidth", "8",
    "-draw", polygon,
    "-depth", "8",
    logo_path
  ))

  if (isTRUE(show_plot) && requireNamespace("png", quietly = TRUE)) {
    grid::grid.newpage()
    grid::grid.raster(grDevices::as.raster(png::readPNG(logo_path)))
  }

  invisible(list(logo = logo_path))
}

make_vnorm_logo <- function(
    out_dir = "man/figures",
    bg = "#F3F4F6",
    curve_color = "#EAF2FF",
    curve_understroke_color = "#FFFFFF",
    curve_double_stroke = FALSE,
    points_color = "#F59E0B",
    hex_border = "#4FD1C5",
    text_color = "#EAF2FF",
    point_alpha = 0.18,
    point_size = 0.22,
    curve_linewidth = 0.78,
    curve_understroke_linewidth = 1.10,
    curve_n = 701,
    n_points = 2600,
    sample_sd = 0.045,
    sample_w = 1.6,
    dpi = 600,
    seed = 20260530,
    show_plot = TRUE,
    text_on_top = TRUE,
    plot_xlim = c(-1.43, 1.43),
    plot_ylim = c(-0.78, 0.78),
    subplot_x = 1,
    subplot_y_top = 0.735,
    subplot_y_bottom = 0.95,
    subplot_width_top = 1.56,
    subplot_width_bottom = 0.96,
    subplot_height_top = 1.06,
    subplot_height_bottom = 0.90,
    text_size = 18
  ) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required.", call. = FALSE)
  }
  if (!requireNamespace("mpoly", quietly = TRUE)) {
    stop("Package 'mpoly' is required.", call. = FALSE)
  }

  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  theme_logo <- ggplot2::theme_void() +
    ggplot2::theme(
      # Keep subplot transparent so no rectangular panel bleeds outside the hex.
      plot.background = ggplot2::element_rect(fill = "transparent", colour = NA),
      panel.background = ggplot2::element_rect(fill = "transparent", colour = NA),
      legend.position = "none",
      plot.margin = ggplot2::margin(0, 0, 0, 0)
    )

  p <- mpoly::mp("(x^2 + y^2)^2 - 2*(x^2 - y^2)")

  if (!is.null(seed)) {
    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    } else {
      NULL
    }
    on.exit({
      if (is.null(old_seed)) {
        if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
          rm(".Random.seed", envir = .GlobalEnv)
        }
      } else {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      }
    }, add = TRUE)
    set.seed(seed)
  }

  pts <- tryCatch(
    vnorm::rvnorm(
      n = n_points,
      poly = p,
      sd = sample_sd,
      output = "tibble",
      rejection = TRUE,
      w = sample_w,
      verbose = FALSE
    ),
    error = function(e) {
      message("rvnorm point generation failed for logo: ", conditionMessage(e))
      NULL
    }
  )

  layers <- list()
  if (!is.null(pts)) {
    layers <- c(layers, list(
      ggplot2::geom_point(
        data = pts,
        ggplot2::aes(x, y),
        colour = points_color,
        alpha = point_alpha,
        size = point_size,
        shape = 16
      )
    ))
  }

  if (isTRUE(curve_double_stroke)) {
    layers <- c(layers, list(
      vnorm::geom_variety(
        poly = p,
        xlim = plot_xlim,
        ylim = plot_ylim,
        n = curve_n,
        colour = curve_understroke_color,
        linewidth = curve_understroke_linewidth,
        lineend = "round",
        linejoin = "round",
        show.legend = FALSE
      )
    ))
  }

  layers <- c(layers, list(
    vnorm::geom_variety(
      poly = p,
      xlim = plot_xlim,
      ylim = plot_ylim,
      n = curve_n,
      colour = curve_color,
      linewidth = curve_linewidth,
      lineend = "round",
      linejoin = "round",
      show.legend = FALSE
    ),
    ggplot2::coord_equal(),
    theme_logo
  ))

  g <- Reduce(`+`, c(list(ggplot2::ggplot()), layers))

  preview_path <- file.path(out_dir, "logo-preview.png")
  ggplot2::ggsave(
    filename = preview_path,
    plot = g,
    width = 4,
    height = 4.6,
    dpi = dpi,
    bg = bg,
    device = if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_png else "png"
  )

  logo_path <- file.path(out_dir, "logo.png")
  if (requireNamespace("hexSticker", quietly = TRUE)) {
    p_y <- if (isTRUE(text_on_top)) 1.53 else 0.18
    s_y <- if (isTRUE(text_on_top)) subplot_y_top else subplot_y_bottom
    s_h <- if (isTRUE(text_on_top)) subplot_height_top else subplot_height_bottom
    s_w <- if (isTRUE(text_on_top)) subplot_width_top else subplot_width_bottom
    hexSticker::sticker(
      subplot = g,
      package = "vnorm",
      filename = logo_path,
      s_x = subplot_x,
      s_y = s_y,
      s_width = s_w,
      s_height = s_h,
      p_x = 1,
      p_y = p_y,
      p_color = text_color,
      p_family = "sans",
      p_size = if (isTRUE(text_on_top)) text_size else 22,
      h_fill = bg,
      h_color = hex_border,
      dpi = dpi
    )
    message("Saved preview to ", preview_path, " and hex logo to ", logo_path)
  } else if (make_vnorm_logo_magick(
    preview_path = preview_path,
    logo_path = logo_path,
    bg = bg,
    hex_border = hex_border,
    text_color = text_color
  )) {
    message("Saved preview to ", preview_path, " and hex logo to ", logo_path)
  } else {
    message(
      "Saved preview to ", preview_path, ". ",
      "Install 'hexSticker' or ImageMagick to generate the final hex logo PNG."
    )
  }

  if (isTRUE(show_plot)) {
    print(g)
  }

  invisible(list(plot = g, preview = preview_path, logo = logo_path))
}

make_vnorm_logo_palette <- function(
    palette = c("neon-torus", "clean-dark", "teal-gold", "crimson-cyan", "mint-copper"),
    ...
  ) {
  palette <- match.arg(palette)

  if (identical(palette, "neon-torus")) {
    return(make_vnorm_logo_neon(...))
  }

  args <- switch(
    palette,
    "clean-dark" = list(
      bg = "#111827",
      curve_color = "#F8FAFC",
      curve_understroke_color = "#111827",
      curve_double_stroke = TRUE,
      points_color = "#F97316",
      hex_border = "#14B8A6",
      text_color = "#F8FAFC",
      point_alpha = 0.42,
      point_size = 0.32,
      curve_linewidth = 0.86,
      curve_understroke_linewidth = 1.25,
      sample_sd = 0.042,
      subplot_y_top = 0.74,
      subplot_width_top = 1.60,
      subplot_height_top = 1.06,
      text_size = 19
    ),
    "teal-gold" = list(
      bg = "#2B2F36",
      curve_color = "#EAF2FF",
      points_color = "#F59E0B",
      hex_border = "#4FD1C5",
      text_color = "#EAF2FF"
    ),
    "crimson-cyan" = list(
      bg = "#111827",
      curve_color = "#F9FAFB",
      points_color = "#F43F5E",
      hex_border = "#22D3EE",
      text_color = "#F9FAFB"
    ),
    "mint-copper" = list(
      bg = "#0F172A",
      curve_color = "#E2FDF7",
      points_color = "#C08457",
      hex_border = "#7DD3C7",
      text_color = "#E2FDF7"
    )
  )

  do.call(make_vnorm_logo, c(args, list(...)))
}
