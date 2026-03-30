# unexported functions from ggplot2

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

tibble0 <- function(...) {
  tibble::tibble(..., .name_repair = "minimal")
}

unique0 <- function(x, ...) {
  if (is.null(x)) x else vctrs::vec_unique(x, ...)
}

isoband_z_matrix <- function(data) {
  # convert long xyz table into a regular z-matrix for contouring
  x_pos <- as.integer(factor(data$x, levels = sort(unique0(data$x))))
  y_pos <- as.integer(factor(data$y, levels = sort(unique0(data$y))))
  nrow <- max(y_pos)
  ncol <- max(x_pos)
  raster <- matrix(NA_real_, nrow = nrow, ncol = ncol)
  raster[cbind(y_pos, x_pos)] <- data$z
  raster
}

xyz_to_isolines <- function(data, breaks) {
  # compute contour lines at the requested break levels
  x <- sort(unique0(data$x))
  y <- sort(unique0(data$y))
  z <- isoband_z_matrix(data)
  z <- break_contour_ties(z, x = x, y = y, breaks = breaks)

  # contourLines expects rows to index x and columns to index y
  grDevices::contourLines(
    x = x,
    y = y,
    z = t(z),
    levels = breaks
  )
}

break_contour_ties <- function(z, x, y, breaks) {
  # contourLines is unstable when grid nodes land exactly on the contour level
  # add a tiny deterministic perturbation only at exact/near-exact ties
  if (length(breaks) != 1 || !is.finite(breaks)) return(z)
  lev <- breaks[[1]]
  if (!any(is.finite(z))) return(z)

  scale_z <- max(abs(z - lev), na.rm = TRUE)
  tol <- sqrt(.Machine$double.eps) * max(1, scale_z)
  ties <- is.finite(z) & abs(z - lev) <= tol
  if (!any(ties)) return(z)

  # z has rows indexed by y and columns indexed by x
  # irrational offset (sqrt(2)) ensures no two grid points get the same perturbation sign
  phase <- outer(y, x, function(yy, xx) xx + sqrt(2) * yy)
  sgn <- sign(phase)
  sgn[sgn == 0] <- 1

  z[ties] <- z[ties] + 4 * tol * sgn[ties]
  z
}

iso_to_path <- function(iso, group = 1) {
  # convert contourLines output into a path tibble for GeomPath
  if (length(iso) == 0) {
    message("Zero contours were generated")
    return(tibble0())
  }

  lengths <- vapply(iso, function(x) length(x$x), integer(1))
  levels <- vapply(iso, function(x) as.character(x$level), character(1))
  xs <- unlist(lapply(iso, "[[", "x"), use.names = FALSE)
  ys <- unlist(lapply(iso, "[[", "y"), use.names = FALSE)
  item_id <- rep(seq_along(iso), lengths)

  groups <- paste(group, sprintf("%03d", item_id), sep = "-")
  groups <- factor(groups)

  tibble0(
    level = rep(levels, lengths),
    x = xs,
    y = ys,
    piece = as.integer(groups),
    group = groups,
    .size = length(xs)
  )
}

empty <- function(df) {
  is.null(df) || nrow(df) == 0 || ncol(df) == 0 || inherits(df, "waiver")
}

ensure_nonempty_data <- function(data) {
  if (empty(data)) {
    tibble0(group = 1, .size = 1)
  } else {
    data
  }
}
