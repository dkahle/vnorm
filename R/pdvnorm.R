#' Pseudo-Density for the Variety Normal Distribution
#'
#' Evaluate the variety normal pseudo-density in either the homoskedastic or
#' heteroskedastic setting.
#'
#' @details A single `mpoly` is normalized as a one-dimensional normal kernel in
#'   the polynomial value, after optional homoskedastic gradient normalization.
#'   A length-one `mpolyList` is treated identically to the same single `mpoly`.
#'
#'   Coordinates are interpreted in sorted variable order,
#'   `sort(mpoly::vars(poly))`, regardless of the order in which variables
#'   appear in the polynomial's terms. If `x` has names (vector) or column
#'   names (matrix/data frame) covering all of the variables, entries are
#'   matched by name instead.
#'
#' @param x A numeric vector of length equal to the number of variables in
#'   `poly`, or a numeric matrix/data frame with that many columns (one row per
#'   evaluation point). Entries/columns are interpreted in sorted variable
#'   order (`sort(mpoly::vars(poly))`); if `x` carries names covering all of
#'   the variables, they are matched by name instead.
#' @param poly An `mpoly` object (single polynomial) or an `mpolyList` object
#'   (multiple polynomials).
#' @param sd Scale parameter for the normal kernel. For polynomial systems,
#'   a scalar is recycled across the relevant dimension, while a vector supplies
#'   component-wise standard deviations.
#' @param Sigma Full covariance matrix, scalar covariance, or a diagonal vector
#'   of covariance terms. If supplied, `Sigma` replaces `sd`. When
#'   `homo = TRUE`, its rows/columns correspond to the variables in sorted
#'   order (`sort(mpoly::vars(poly))`).
#' @param homo Logical; default is `TRUE`. If `TRUE`, compute the homoskedastic
#'   variety normal pseudo-density. If `FALSE`, compute the heteroskedastic
#'   pseudo-density.
#' @param log Logical. If `TRUE`, returns the log of the density.
#' @param ... Deprecated. A named `sigma` argument is accepted for backward
#'   compatibility.
#' @return A numeric scalar or vector containing the pseudo-density evaluated at
#'   `x`.
#'
#' @examples
#'
#' # m = 1 polynomial in n = 1 variable, 0d variety at x = 0
#' p <- mp("x")
#' pdvnorm(c(-1, 0, 1), poly = p, sd = 1)
#' pdvnorm(c(-1, 0, 1), poly = p, sd = 1, log = TRUE)
#'
#'
#'
#' # m = 2 polynomials in n = 2 variables (square system), 0d variety at (0, 0)
#' p <- mp(c("x", "y"))
#' X <- rbind(c(-.25, 0), c(1, 2), c(-1, 3))
#' pdvnorm(X, poly = p, sd = 1)
#'
#'
#'
#' # m = 1 polynomial in n = 2 variables, 1d variety at unit circle
#' p <- mp("x^2 + y^2 - 1")
#' X <- rbind(c(-.25, 0), c(1, 2), c(-1, 3), sqrt(2)/2 * c(1,1))
#' pdvnorm(X, poly = p, sd = 1)
#'
#'
#'
#' # different ways to specify dispersion
#' p <- mp(c("x", "y"))
#' X <- rbind(c(0, 0), c(1, 2), c(-1, 3))
#' pdvnorm(X, p, sd = 1)
#' pdvnorm(X, p, sd = c(1, 2))
#' pdvnorm(X, p, Sigma = diag(c(1, 4)))
#'
#'
#'
#' # multivariate (underdetermined): one polynomial in two variables
#' p <- mp("x^2 + y^2 - 1")
#' X <- rbind(c(1, 1), c(2, -1), c(0, 3))
#' pdvnorm(X, p, sd = 1)
#' pdvnorm(as.data.frame(X), p, sd = 1)
#' pdvnorm(c(1, 1), p, Sigma = diag(c(1, 4)))
#'
#'
#'
#' # multivariate (overdetermined): three polynomials in two variables
#' p <- mp(c("x", "y", "x + y"))
#' X <- rbind(c(1, 2), c(0, -1), c(2, 2))
#' pdvnorm(X, p, Sigma = diag(2), homo = TRUE)
#' pdvnorm(X, p, sd = c(1, 2, 3), homo = FALSE)
#'
#'
#'
#' \dontrun{
#'
#' library("ggplot2")
#'
#' # ggplot2-only helpers for visualizing pdvnorm in 1d and 2d
#' pdf_curve <- function(poly, xlim, n = 501, ...) {
#'   df <- data.frame(x = seq(xlim[1], xlim[2], length.out = n))
#'   df$density <- pdvnorm(df$x, poly, ...)
#'   geom_line(aes(x, density), data = df)
#' }
#'
#' pdf_raster <- function(poly, xlim, ylim, n = 201, ...) {
#'   grid <- expand.grid(
#'     x = seq(xlim[1], xlim[2], length.out = n),
#'     y = seq(ylim[1], ylim[2], length.out = n)
#'   )
#'   grid$density <- pdvnorm(as.matrix(grid), poly, ...)
#'   list(
#'     geom_raster(aes(x, y, fill = density), data = grid),
#'     scale_fill_viridis_c(option = "magma", guide = "none")
#'   )
#' }
#'
#'
#' # m = 1 polynomial in n = 1 variable, 0d variety at x = 0
#' p <- mp("x")
#'
#' ggplot() +
#'   pdf_curve(p, xlim = c(-5, 5), sd = 1)
#'
#'
#'
#' # m = 2 polynomials in n = 2 variables (square system), 0d variety at (0, 0)
#' p <- mp(c("x", "y"))
#'
#' ggplot() +
#'   pdf_raster(p, xlim = c(-5, 5), ylim = c(-5, 5), sd = 1) +
#'   coord_equal()
#'
#'
#'
#' # m = 1 polynomial in n = 2 variables, 1d variety at unit circle
#' p <- mp("x^2 + y^2 - 1")
#'
#' ggplot() +
#'   pdf_raster(p, xlim = c(-1.5, 1.5), ylim = c(-1.5, 1.5), sd = .1) +
#'   coord_equal()
#'
#'
#'
#' # m = 1 polynomial in n = 2 variables, 1d variety
#' p <- mp("y^2 - x^2 (x + 1)")
#'
#' ggplot() +
#'   pdf_raster(p, xlim = c(-1.5, 1.5), ylim = c(-1.5, 1.5), sd = .1) +
#'   coord_equal()
#'
#' ggplot() +
#'   pdf_raster(p, xlim = c(-1.5, 1.5), ylim = c(-1.5, 1.5), sd = .1, homo = FALSE) +
#'   coord_equal()
#'
#'
#'
#' # m = 1 polynomial in n = 2 variables, 1d varieties
#' si <- .025
#' w <- 1.5
#'
#' p1 <- ggplot() +
#'   pdf_raster(mp("x^2 + (4 y)^2 - 1"), c(-w, w), c(-w, w), sd = si) +
#'   coord_equal(xlim = c(-w, w), ylim = c(-w, w))
#'
#' p2 <- ggplot() +
#'   pdf_raster(mp("-1 (x - y) (x + y)"), c(-w, w), c(-w, w), sd = si) +
#'   coord_equal(xlim = c(-w, w), ylim = c(-w, w))
#'
#' p3 <- ggplot() +
#'   pdf_raster(mp("(x^2 + y^2)^3 - 4 x^2 y^2"), c(-w, w), c(-w, w), sd = si) +
#'   coord_equal(xlim = c(-w, w), ylim = c(-w, w))
#'
#' p4 <- ggplot() +
#'   pdf_raster(mp("(x^2 + y^2 - 1)^3 - x^2 y^3"), c(-w, w), c(-w, w), sd = si) +
#'   coord_equal(xlim = c(-w, w), ylim = c(-w, w))
#'
#' library("patchwork")
#' p1 + p2 + p3 + p4 +
#'   plot_layout(nrow = 1, widths = 1, guides = "collect")
#'
#'
#'
#' # m = 2 polynomials in n = 2 variables, 0d varieties
#' # different representations and "covariances"
#'
#' Sigma <- .2^2
#' w <- 1.5
#'
#' p1 <- ggplot() +
#'   pdf_raster(mp(c("x", "y")), c(-w, w), c(-w, w), Sigma = Sigma) +
#'   coord_equal(xlim = c(-w, w), ylim = c(-w, w)) +
#'   labs(title = latex2exp::TeX(r"(V(x, y))"))
#'
#' p2 <- ggplot() +
#'   pdf_raster(mp(c("-1 (x - y)", "(x + y)")), c(-w, w), c(-w, w), Sigma = Sigma) +
#'   coord_equal(xlim = c(-w, w), ylim = c(-w, w)) +
#'   labs(title = latex2exp::TeX(r"(V(y - x, y + x))"))
#'
#' p3 <- ggplot() +
#'   pdf_raster(mp(c("(-1 (x - y))^2", "(x + y)^2")), c(-w, w), c(-w, w), Sigma = Sigma) +
#'   coord_equal(xlim = c(-w, w), ylim = c(-w, w)) +
#'   labs(title = latex2exp::TeX(r"(V((y - x)^2, (y + x)^2))"))
#'
#' p4 <- ggplot() +
#'   pdf_raster(mp(c("x y^3 - x^3 y", "x^2 + y^2 - 1")), c(-w, w), c(-w, w), Sigma = Sigma) +
#'   coord_equal(xlim = c(-w, w), ylim = c(-w, w)) +
#'   labs(title = latex2exp::TeX(r"(V(x y^3 - x^3 y, x^2 + y^2 - 1))"))
#'
#' p1 + p2 + p3 + p4 +
#'   plot_layout(nrow = 1, widths = 1, guides = "collect")
#'
#'
#'
#' S <- (.20 * diag(c(1, 1)))
#' R <- matrix(c(1, .90, .90, 1), nrow = 2)
#' Sigma <- S %*% R %*% S
#' w <- 1.5
#'
#' p1 <- ggplot() +
#'   pdf_raster(mp(c("x", "y")), c(-w, w), c(-w, w), Sigma = Sigma) +
#'   coord_equal(xlim = c(-w, w), ylim = c(-w, w)) +
#'   labs(title = latex2exp::TeX(r"(V(x, y))"))
#'
#' p2 <- ggplot() +
#'   pdf_raster(mp(c("-1 (x - y)", "(x + y)")), c(-w, w), c(-w, w), Sigma = Sigma) +
#'   coord_equal(xlim = c(-w, w), ylim = c(-w, w)) +
#'   labs(title = latex2exp::TeX(r"(V(y - x, y + x))"))
#'
#' p3 <- ggplot() +
#'   pdf_raster(mp(c("(-1 (x - y))^2", "(x + y)^2")), c(-w, w), c(-w, w), Sigma = Sigma) +
#'   coord_equal(xlim = c(-w, w), ylim = c(-w, w)) +
#'   labs(title = latex2exp::TeX(r"(V((y - x)^2, (y + x)^2))"))
#'
#' p4 <- ggplot() +
#'   pdf_raster(mp(c("x y^3 - x^3 y", "x^2 + y^2 - 1")), c(-w, w), c(-w, w), Sigma = Sigma) +
#'   coord_equal(xlim = c(-w, w), ylim = c(-w, w)) +
#'   labs(title = latex2exp::TeX(r"(V(x y^3 - x^3 y, x^2 + y^2 - 1))"))
#'
#' p1 + p2 + p3 + p4 +
#'   plot_layout(nrow = 1, widths = 1, guides = "collect")
#'
#'
#'
#' # geom_hdr_fun accepts functions f(x, y) with vectors x and y,
#' # but pdvnorm accepts the packed f(matrix), so a wrapper is needed
#' library("ggdensity")
#'
#' p <- mp("x^2 + y^2 - 1")
#' f <- function(x, y, ...) pdvnorm(cbind(x, y), ...)
#'
#' f(1, 2, poly = p, sd = .05)
#' f(1:2, 2:3, poly = p, sd = .05)
#'
#' ggplot() +
#'   geom_hdr_fun(
#'     fun = f, xlim = c(-1.25, 1.25), ylim = c(-1.25, 1.25),
#'     args = list(poly = p, sd = .10)
#'   ) +
#'   coord_equal()
#'
#' S <- (.10 * diag(c(1, 1)))
#' R <- matrix(c(1, .95, .95, 1), nrow = 2)
#'
#' ggplot() +
#'   geom_hdr_fun(
#'     fun = f, xlim = c(-1.25, 1.25), ylim = c(-1.25, 1.25),
#'     args = list(poly = p, Sigma = S %*% R %*% S), n = 250
#'   ) +
#'   coord_equal()
#' }
#'
#' @export
pdvnorm <- function(x, poly, sd, homo = TRUE, log = FALSE, Sigma = NULL, ...) {
  # dispatch between single-polynomial and polynomial-list density evaluation
  is_uni <- inherits(poly, "mpoly")
  is_multi <- inherits(poly, "mpolyList")
  if (!(is_uni || is_multi)) {
    stop(
      paste0(
        "'poly' must be either an 'mpoly' (univariate) or an ",
        "'mpolyList' (multivariate)."
      )
    )
  }

  if (length(mpoly::vars(poly)) == 0L) {
    stop("`poly` must contain at least one variable.", call. = FALSE)
  }

  has_sd <- !missing(sd)
  legacy_sigma <- NULL
  legacy_sigma_supplied <- FALSE
  dots <- list(...)
  if (length(dots) > 0L) {
    dot_names <- names(dots)
    if (
      is.null(dot_names) ||
        any(!nzchar(dot_names)) ||
        any(dot_names != "sigma") ||
        length(dots) > 1L
    ) {
      stop(
        "Unused argument(s): ",
        paste(ifelse(nzchar(dot_names), dot_names, "<unnamed>"), collapse = ", "),
        call. = FALSE
      )
    }
    if (has_sd || !is.null(Sigma)) {
      stop("Specify only one of `sd`, `Sigma`, or `sigma`.", call. = FALSE)
    }
    legacy_sigma <- dots$sigma
    legacy_sigma_supplied <- TRUE
  }

  if (is_multi && length(poly) == 1L) {
    args <- list(
      x = x,
      poly = poly[[1]],
      homo = homo,
      log = log,
      Sigma = Sigma
    )
    if (has_sd) args$sd <- sd
    if (legacy_sigma_supplied) args$sigma <- legacy_sigma
    return(do.call(pdvnorm, args))
  }

  if (is_uni) {
    # univariate/single-polynomial path
    # coordinates follow sorted variable order unless x carries variable names
    vs <- sort(mpoly::vars(poly))
    n <- length(vs)
    if (is.data.frame(x)) {
      x <- as.matrix(x)
    }
    if (is.vector(x)) {
      if (n == 1) {
        x <- matrix(x, ncol = 1)
      } else {
        x <- pdvnorm_align_x(x, vs)
        if (length(x) != n) {
          stop("x must be length n or a matrix with n columns.")
        }
        x <- matrix(x, nrow = 1)
      }
    } else if (is.matrix(x)) {
      x <- pdvnorm_align_x(x, vs)
      if (ncol(x) != n) {
        stop("x must be length n or a matrix with n columns.")
      }
    } else {
      stop("x must be length n or a matrix with n columns.")
    }
    if (!is.numeric(x) || any(!is.finite(x))) {
      stop("'x' must be finite numeric.")
    }

    g_func <- poly_function(poly, varorder = vs, silent = TRUE)

    if (homo) {
      grad_g_obj <- deriv(poly, var = vs)
      if (all(is.constant(grad_g_obj))) {
        const <- unname(unlist(grad_g_obj))

        # force() ensures const is evaluated in this iteration, not lazily
        force(const)
        grad_g <- function(...) const
      } else {
        grad_g <- poly_function(grad_g_obj, varorder = vs, silent = TRUE)
      }
    }

    if (
      !legacy_sigma_supplied &&
        !is.null(Sigma) &&
        !pdvnorm_is_scalar(Sigma)
    ) {
      if (!homo) {
        stop(
          "`Sigma` matrices for single-polynomial densities require `homo = TRUE`.",
          call. = FALSE
        )
      }
      covariance <- pdvnorm_covariance_from_Sigma(
        Sigma, homo = TRUE, n = n, m = 1L, label = "`Sigma`"
      )
      tryCatch(
        chol(covariance$Sigma),
        error = function(e) {
          stop("`Sigma` must be positive definite.", call. = FALSE)
        }
      )

      log_density <- apply(x, 1, function(row_vec) {
        g_val <- g_func(row_vec)
        grad_g_val <- as.numeric(grad_g(row_vec))
        grad_g_norm <- sqrt(sum(grad_g_val^2))
        if (grad_g_norm == 0) {
          # at a critical point the normalized residual g/||grad g|| diverges
          # unless the point is on the variety, where the density blows up
          return(if (g_val == 0) Inf else -Inf)
        }
        normal_direction <- grad_g_val / grad_g_norm
        normal_var <- as.numeric(
          t(normal_direction) %*% covariance$Sigma %*% normal_direction
        )
        normal_sd <- sqrt(normal_var)
        z <- (g_val / grad_g_norm) / normal_sd
        -(0.5 * z^2 + base::log(normal_sd) + 0.5 * base::log(2 * base::pi))
      })
      return(if (log) log_density else base::exp(log_density))
    }

    sd_value <- if (legacy_sigma_supplied) {
      legacy_sigma
    } else if (!is.null(Sigma)) {
      pdvnorm_single_sd_from_Sigma(Sigma)
    } else {
      if (!has_sd) {
        stop("`sd` must be supplied unless `Sigma` is supplied.", call. = FALSE)
      }
      sd
    }
    if (
      !is.numeric(sd_value) ||
        length(sd_value) != 1L ||
        !is.finite(sd_value) ||
        sd_value <= 0
    ) {
      stop(
        "For the single-polynomial case, `sd` must be a single positive numeric.",
        call. = FALSE
      )
    }
    sd <- as.numeric(sd_value)

    log_density <- apply(x, 1, function(row_vec) {
      g_val <- g_func(row_vec)
      if (homo) {
        grad_g_val <- sqrt(sum(grad_g(row_vec)^2))
        if (grad_g_val == 0) {
          # at a critical point the normalized residual g/||grad g|| diverges
          # unless the point is on the variety, where the density blows up
          return(if (g_val == 0) Inf else -Inf)
        }
        g_val <- g_val / grad_g_val
      }
      z <- g_val / sd
      -(0.5 * z^2 + base::log(sd) + 0.5 * base::log(2 * base::pi))
    })
    return(if (log) log_density else base::exp(log_density))
  }

  # multivariate/polynomial-list path
  # coordinates follow sorted variable order unless x carries variable names
  vars <- sort(mpoly::vars(poly))
  n <- length(vars)
  m <- length(poly)

  X <- if (is.null(dim(x))) {
    matrix(as.numeric(pdvnorm_align_x(x, vars)), nrow = 1)
  } else {
    pdvnorm_align_x(as.matrix(x), vars)
  }
  if (!is.numeric(X) || any(!is.finite(X))) {
    stop("'x' must be finite numeric.")
  }
  if (ncol(X) != n) {
    stop("'x' must have length n (vector) or n columns (matrix).")
  }

  if (legacy_sigma_supplied) {
    covariance <- pdvnorm_covariance_from_Sigma(
      legacy_sigma, homo = homo, n = n, m = m, label = "`sigma`"
    )
  } else if (!is.null(Sigma)) {
    covariance <- pdvnorm_covariance_from_Sigma(
      Sigma, homo = homo, n = n, m = m, label = "`Sigma`"
    )
  } else {
    if (!has_sd) {
      stop("`sd` must be supplied unless `Sigma` is supplied.", call. = FALSE)
    }
    covariance <- pdvnorm_covariance_from_sd(sd, homo = homo, n = n, m = m)
  }
  Sigma <- covariance$Sigma
  dispersion_label <- covariance$label

  g_fns <- poly_function(poly, varorder = vars, silent = TRUE)
  g_vals_mat <- g_fns(X)
  if (m == 1) {
    g_vals_mat <- matrix(g_vals_mat, ncol = 1)
  }

  # cholesky factor for efficient Mahalanobis distance
  L <- tryCatch(
    chol(Sigma),
    error = function(e) {
      stop(dispersion_label, " must be positive definite.", call. = FALSE)
    }
  )
  log_det_sigma <- 2 * sum(base::log(diag(L)))

  if (!homo) {
    solved <- backsolve(L, t(g_vals_mat), transpose = TRUE)
    quad <- colSums(solved^2)
    out_log <- -0.5 * (m * base::log(2 * base::pi) + log_det_sigma + quad)
    return(if (log) out_log else base::exp(out_log))
  }

  grad_fun <- vector("list", m)
  constant_jacobian <- TRUE
  constant_J <- matrix(NA_real_, nrow = m, ncol = n)
  for (i in seq_len(m)) {
    grad_obj <- deriv(poly[[i]], var = vars)
    if (all(is.constant(grad_obj))) {
      const <- unname(unlist(grad_obj))
      constant_J[i, ] <- const
      grad_fun[[i]] <- local({
        const_i <- const
        force(const_i)
        function(...) const_i
      })
    } else {
      constant_jacobian <- FALSE
      grad_fun[[i]] <- poly_function(grad_obj, varorder = vars, silent = TRUE)
    }
  }

  if (constant_jacobian) {
    Jp <- pdvnorm_pseudoinverse(constant_J)
    v_mat <- t(Jp %*% t(g_vals_mat))
    solved <- backsolve(L, t(v_mat), transpose = TRUE)
    quad <- colSums(solved^2)
    out_log <- -0.5 * (n * base::log(2 * base::pi) + log_det_sigma + quad)
    return(if (log) out_log else base::exp(out_log))
  }

  out_log <- numeric(nrow(X))
  for (i in seq_len(nrow(X))) {
    xi <- as.numeric(X[i, ])
    g_vals <- as.numeric(g_vals_mat[i, ])

    # jacobian of g(x): rows are equations, columns are variables
    J <- matrix(NA_real_, nrow = m, ncol = n)
    for (j in seq_len(m)) {
      J[j, ] <- grad_fun[[j]](xi)
    }

    Jp <- pdvnorm_pseudoinverse(J)
    v <- Jp %*% g_vals
    quad <- sum(backsolve(L, v, transpose = TRUE)^2)
    out_log[i] <- -0.5 * (n * base::log(2 * base::pi) + log_det_sigma + quad)
  }

  if (log) out_log else base::exp(out_log)
}

pdvnorm_pseudoinverse <- function(J) {
  m <- nrow(J)
  n <- ncol(J)
  sv <- svd(J)
  tol <- max(dim(J)) *
    .Machine$double.eps *
    ifelse(length(sv$d) > 0, sv$d[1], 0)
  r <- sum(sv$d > tol)

  if (n > m && r == m) {
    # full row rank (underdetermined): right pseudoinverse
    t(J) %*% solve(J %*% t(J))
  } else if (m > n && r == n) {
    # full column rank (overdetermined): left pseudoinverse
    solve(t(J) %*% J) %*% t(J)
  } else if (m == n && r == n) {
    # square and full rank: exact inverse
    solve(J)
  } else {
    # rank-deficient fallback: SVD pseudoinverse with tolerance cutoff
    dplus <- ifelse(sv$d > tol, 1 / sv$d, 0)
    sv$v %*% (dplus * t(sv$u))
  }
}

pdvnorm_align_x <- function(x, vs) {
  # match x entries/columns to variables by name when variable names are
  # present; names disjoint from the variables are treated as positional
  nms <- if (is.null(dim(x))) names(x) else colnames(x)
  if (is.null(nms) || length(intersect(nms, vs)) == 0L) {
    return(x)
  }
  if (anyDuplicated(nms) > 0L) {
    stop("`x` has duplicated column names.", call. = FALSE)
  }
  missing_vars <- setdiff(vs, nms)
  if (length(missing_vars) > 0L) {
    stop(
      "`x` names match some but not all variables; missing: ",
      paste0("`", missing_vars, "`", collapse = ", "), ".",
      call. = FALSE
    )
  }
  extra <- setdiff(nms, vs)
  if (length(extra) > 0L) {
    stop(
      "`x` has unexpected column(s): ",
      paste0("`", extra, "`", collapse = ", "), ".",
      call. = FALSE
    )
  }
  if (is.null(dim(x))) x[vs] else x[, vs, drop = FALSE]
}

pdvnorm_is_scalar <- function(x) {
  length(x) == 1L
}

pdvnorm_single_sd_from_Sigma <- function(Sigma) {
  if (!is.numeric(Sigma) || any(!is.finite(Sigma))) {
    stop("`Sigma` must be finite numeric.", call. = FALSE)
  }
  if (length(Sigma) != 1L) {
    stop(
      "For the single-polynomial case, `Sigma` must be a single positive variance.",
      call. = FALSE
    )
  }
  if (Sigma <= 0) {
    stop("`Sigma` must be positive.", call. = FALSE)
  }
  sqrt(as.numeric(Sigma))
}

pdvnorm_covariance_from_sd <- function(sd, homo, n, m) {
  q <- if (homo) n else m
  dim_label <- if (homo) "n" else "m"

  if (!is.numeric(sd) || any(!is.finite(sd))) {
    stop("`sd` must be finite numeric.", call. = FALSE)
  }
  if (is.matrix(sd)) {
    stop("Use `Sigma` to supply a covariance matrix.", call. = FALSE)
  }
  if (length(sd) == 1L) {
    if (sd <= 0) {
      stop("`sd` must be positive.", call. = FALSE)
    }
    return(list(Sigma = diag(as.numeric(sd)^2, q), label = "`sd`"))
  }
  if (!is.vector(sd)) {
    stop("`sd` must be a scalar or vector.", call. = FALSE)
  }
  if (any(sd <= 0)) {
    stop("`sd` vector entries must be positive.", call. = FALSE)
  }
  if (length(sd) != q) {
    stop(
      "When homo=", homo, ", length(`sd`) must be ", dim_label, ".",
      call. = FALSE
    )
  }
  list(Sigma = diag(as.numeric(sd)^2), label = "`sd`")
}

pdvnorm_covariance_from_Sigma <- function(Sigma, homo, n, m, label) {
  q <- if (homo) n else m
  dim_label <- if (homo) "n" else "m"

  if (!is.numeric(Sigma) || any(!is.finite(Sigma))) {
    stop(label, " must be finite numeric.", call. = FALSE)
  }
  if (is.null(dim(Sigma)) && length(Sigma) == 1L) {
    if (Sigma <= 0) {
      stop(label, " must be positive.", call. = FALSE)
    }
    return(list(Sigma = diag(as.numeric(Sigma), q), label = label))
  }
  if (is.null(dim(Sigma)) && is.vector(Sigma)) {
    if (any(Sigma <= 0)) {
      stop(label, " vector entries must be positive.", call. = FALSE)
    }
    if (length(Sigma) != q) {
      stop(
        "When homo=", homo, ", length(", label, ") must be ", dim_label, ".",
        call. = FALSE
      )
    }
    return(list(Sigma = diag(as.numeric(Sigma)), label = label))
  }

  Sigma <- as.matrix(Sigma)
  if (!all(dim(Sigma) == c(q, q))) {
    stop(
      label, " matrix must be ", dim_label, " x ", dim_label,
      " when homo=", homo, ".",
      call. = FALSE
    )
  }
  list(Sigma = Sigma, label = label)
}
