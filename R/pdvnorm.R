#' Pseudo-Density for the Variety Normal Distribution
#'
#' Evaluate the variety normal pseudo-density in either the homoskedastic or
#' heteroskedastic setting.
#'
#' @param x A numeric vector of length equal to the number of variables in
#'   `poly`, or a numeric matrix/data frame with that many columns (one row per
#'   evaluation point).
#' @param poly An `mpoly` object (single polynomial) or an `mpolyList` object
#'   (multiple polynomials).
#' @param sd Scale parameter for the normal kernel. For polynomial systems,
#'   a scalar is recycled across the relevant dimension, while a vector supplies
#'   component-wise standard deviations.
#' @param Sigma Full covariance matrix, scalar covariance, or a diagonal vector
#'   of covariance terms. If supplied, `Sigma` replaces `sd`.
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
#' # m = 1 polynomial in n = 1 variable, 0d variety
#' p <- mp("x")
#'
#' pdvnorm(c(-1, 0, 1), poly = p, sd = 1)
#' pdvnorm(c(-1, 0, 1), poly = p, sd = 1, log = TRUE)
#'
#'
#'
#' # m = 2 polynomials in n = 2 variables (square system), 0d variety
#' p <- mp(c("x", "y"))
#'
#' X <- rbind(c(-.25, 0), c(1, 2), c(-1, 3))
#' pdvnorm(X, poly = p, sd = 1)
#'
#'
#'
#' # m = 1 polynomial in n = 2 variables, 1d variety
#' p <- mp("x^2 + y^2 - 1")
#'
#' X <- rbind(c(-.25, 0), c(1, 2), c(-1, 3))
#' pdvnorm(X, poly = p, sd = 1)
#'
#'
#'
#' # different dispersion forms
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
#' library("ggfunction")
#'
#'
#'
#' # m = 1 polynomial in n = 1 variable, 0d variety
#' p <- mp("x")
#'
#' ggplot() +
#'   geom_pdf(fun = pdvnorm, args = list(poly = p, sd = 1), xlim = c(-5, 5))
#'
#'
#'
#' # m = 2 polynomials in n = 2 variables (square system), 0d variety
#' p <- mp(c("x", "y"))
#'
#' ggplot() +
#'   geom_pdf_2d(
#'     fun = pdvnorm, args = list(poly = p, sd = 1),
#'     xlim = c(-5, 5), ylim = c(-5, 5)
#'   ) +
#'   coord_equal()
#'
#'
#'
#' # m = 1 polynomial in n = 2 variables, 1d variety
#' p <- mp("x^2 + y^2 - 1")
#'
#' ggplot() +
#'   geom_pdf_2d(
#'     fun = pdvnorm, args = list(poly = p, sd = .1),
#'     xlim = c(-1.5, 1.5), ylim = c(-1.5, 1.5)
#'   ) +
#'   coord_equal()
#'
#'
#'
#' # m = 1 polynomial in n = 2 variables, 1d variety
#' p <- mp("-1 x^2 (x + 1) + y^2", varorder = c("x", "y"))
#'
#' ggplot() +
#'   geom_pdf_2d(
#'     fun = pdvnorm, args = list(poly = p, sd = .1),
#'     xlim = c(-1.5, 1.5), ylim = c(-1.5, 1.5)
#'   ) +
#'   coord_equal()
#'
#' ggplot() +
#'   geom_pdf_2d(
#'     fun = pdvnorm, args = list(poly = p, sd = .1, homo = FALSE),
#'     xlim = c(-1.5, 1.5), ylim = c(-1.5, 1.5)
#'   ) +
#'   coord_equal()
#'
#'
#'
#' # m = 1 polynomial in n = 2 variables, 1d varieties
#' si <- .025
#' w <- 1.5
#'
#' p1 <- ggplot() +
#'   geom_pdf_2d(
#'     fun = pdvnorm, args = list(poly = mp("x^2 + (4 y)^2 - 1"), sd = si),
#'     xlim = c(-w, w), ylim = c(-w, w)
#'   ) +
#'   coord_equal(xlim = c(-w, w), ylim = c(-w, w))
#'
#' p2 <- ggplot() +
#'   geom_pdf_2d(
#'     fun = pdvnorm, args = list(poly = mp("-1 (x - y) (x + y)"), sd = si),
#'     xlim = c(-w, w), ylim = c(-w, w)
#'   ) +
#'   coord_equal(xlim = c(-w, w), ylim = c(-w, w))
#'
#' p3 <- ggplot() +
#'   geom_pdf_2d(
#'     fun = pdvnorm, args = list(poly = mp("(x^2 + y^2)^3 - 4 x^2 y^2"), sd = si),
#'     xlim = c(-w, w), ylim = c(-w, w)
#'   ) +
#'   coord_equal(xlim = c(-w, w), ylim = c(-w, w))
#'
#' p4 <- ggplot() +
#'   geom_pdf_2d(
#'     fun = pdvnorm, args = list(poly = mp("(x^2 + y^2 - 1)^3 - x^2 y^3"), sd = si),
#'     xlim = c(-w, w), ylim = c(-w, w)
#'   ) +
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
#'   geom_pdf_2d(
#'     fun = pdvnorm, args = list(poly = mp(c("x", "y")), Sigma = Sigma),
#'     xlim = c(-w, w), ylim = c(-w, w)
#'   ) +
#'   coord_equal(xlim = c(-w, w), ylim = c(-w, w)) +
#'   labs(title = latex2exp::TeX(r"(V(x, y))"))
#'
#' p2 <- ggplot() +
#'   geom_pdf_2d(
#'     fun = pdvnorm, args = list(poly = mp(c("-1 (x - y)", "(x + y)")), Sigma = Sigma),
#'     xlim = c(-w, w), ylim = c(-w, w)
#'   ) +
#'   coord_equal(xlim = c(-w, w), ylim = c(-w, w)) +
#'   labs(title = latex2exp::TeX(r"(V(y - x, y + x))"))
#'
#' p3 <- ggplot() +
#'   geom_pdf_2d(
#'     fun = pdvnorm, args = list(poly = mp(c("(-1 (x - y))^2", "(x + y)^2")), Sigma = Sigma),
#'     xlim = c(-w, w), ylim = c(-w, w)
#'   ) +
#'   coord_equal(xlim = c(-w, w), ylim = c(-w, w)) +
#'   labs(title = latex2exp::TeX(r"(V((y - x)^2, (y + x)^2))"))
#'
#' p4 <- ggplot() +
#'   geom_pdf_2d(
#'     fun = pdvnorm, args = list(poly = mp(c("x y^3 - x^3 y", "x^2 + y^2 - 1")), Sigma = Sigma),
#'     xlim = c(-w, w), ylim = c(-w, w)
#'   ) +
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
#'   geom_pdf_2d(
#'     fun = pdvnorm, args = list(poly = mp(c("x", "y")), Sigma = Sigma),
#'     xlim = c(-w, w), ylim = c(-w, w)
#'   ) +
#'   coord_equal(xlim = c(-w, w), ylim = c(-w, w)) +
#'   labs(title = latex2exp::TeX(r"(V(x, y))"))
#'
#' p2 <- ggplot() +
#'   geom_pdf_2d(
#'     fun = pdvnorm, args = list(poly = mp(c("-1 (x - y)", "(x + y)")), Sigma = Sigma),
#'     xlim = c(-w, w), ylim = c(-w, w)
#'   ) +
#'   coord_equal(xlim = c(-w, w), ylim = c(-w, w)) +
#'   labs(title = latex2exp::TeX(r"(V(y - x, y + x))"))
#'
#' p3 <- ggplot() +
#'   geom_pdf_2d(
#'     fun = pdvnorm, args = list(poly = mp(c("(-1 (x - y))^2", "(x + y)^2")), Sigma = Sigma),
#'     xlim = c(-w, w), ylim = c(-w, w)
#'   ) +
#'   coord_equal(xlim = c(-w, w), ylim = c(-w, w)) +
#'   labs(title = latex2exp::TeX(r"(V((y - x)^2, (y + x)^2))"))
#'
#' p4 <- ggplot() +
#'   geom_pdf_2d(
#'     fun = pdvnorm, args = list(poly = mp(c("x y^3 - x^3 y", "x^2 + y^2 - 1")), Sigma = Sigma),
#'     xlim = c(-w, w), ylim = c(-w, w)
#'   ) +
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

  if (is_uni) {
    # univariate/single-polynomial path
    n <- length(mpoly::vars(poly))
    if (is.data.frame(x)) {
      x <- as.matrix(x)
    }
    if (is.vector(x)) {
      if (n == 1) {
        x <- matrix(x, ncol = 1)
      } else {
        if (length(x) != n) {
          stop("x must be length n or a matrix with n columns.")
        }
        x <- matrix(x, nrow = 1)
      }
    } else if (!(is.matrix(x) && ncol(x) == n)) {
      stop("x must be length n or a matrix with n columns.")
    }
    if (!is.numeric(x) || any(!is.finite(x))) {
      stop("'x' must be finite numeric.")
    }

    g_func <- suppressMessages(as.function(poly))

    if (homo) {
      grad_g_obj <- mpoly::gradient(poly)
      if (mean(is.constant(grad_g_obj)) == 1) {
        const <- unname(unlist(grad_g_obj))

        # force() ensures const is evaluated in this iteration, not lazily
        force(const)
        grad_g <- function(...) const
      } else {
        grad_g <- suppressMessages(as.function(grad_g_obj))
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
          return(Inf)
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
          return(Inf)
        }
        g_val <- g_val / grad_g_val
      }
      z <- g_val / sd
      -(0.5 * z^2 + base::log(sd) + 0.5 * base::log(2 * base::pi))
    })
    return(if (log) log_density else base::exp(log_density))
  }

  vars <- mpoly::vars(poly)
  n <- length(vars)
  m <- length(poly)

  # multivariate/polynomial-list path
  X <- if (is.null(dim(x))) matrix(as.numeric(x), nrow = 1) else as.matrix(x)
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

  g_fns <- suppressMessages(as.function(poly))
  g_vals_mat <- t(apply(X, 1, g_fns))
  if (m == 1) {
    g_vals_mat <- matrix(g_vals_mat, ncol = 1)
  }

  grad_fun <- vector("list", m)
  for (i in seq_len(m)) {
    grad_fun[[i]] <- deriv(poly[[i]], var = vars)
    if (mean(is.constant(grad_fun[[i]])) == 1) {
      const <- unname(unlist(grad_fun[[i]]))
      grad_fun[[i]] <- local({
        const_i <- const
        force(const_i)
        function(...) const_i
      })
    } else {
      grad_fun[[i]] <- suppressMessages(
        as.function(deriv(poly[[i]], var = vars), varorder = vars)
      )
    }
  }

  # cholesky factor for efficient Mahalanobis distance
  L <- tryCatch(
    chol(Sigma),
    error = function(e) {
      stop(dispersion_label, " must be positive definite.", call. = FALSE)
    }
  )
  log_det_sigma <- 2 * sum(base::log(diag(L)))

  out_log <- numeric(nrow(X))
  for (i in seq_len(nrow(X))) {
    xi <- as.numeric(X[i, ])
    g_vals <- as.numeric(g_vals_mat[i, ])

    if (homo) {
      # jacobian of g(x): rows are equations, columns are variables
      J <- matrix(NA_real_, nrow = m, ncol = n)
      for (j in seq_len(m)) {
        J[j, ] <- grad_fun[[j]](xi)
      }

      sv <- svd(J)
      tol <- max(dim(J)) *
        .Machine$double.eps *
        ifelse(length(sv$d) > 0, sv$d[1], 0)
      r <- sum(sv$d > tol)

      if (n > m && r == m) {
        # full row rank (underdetermined): right pseudoinverse
        Jp <- t(J) %*% solve(J %*% t(J))
      } else if (m > n && r == n) {
        # full column rank (overdetermined): left pseudoinverse
        Jp <- solve(t(J) %*% J) %*% t(J)
      } else if (m == n && r == n) {
        # square and full rank: exact inverse
        Jp <- solve(J)
      } else {
        # rank-deficient fallback: SVD pseudoinverse with tolerance cutoff
        dplus <- ifelse(sv$d > tol, 1 / sv$d, 0)
        Jp <- sv$v %*% (dplus * t(sv$u))
      }

      v <- Jp %*% g_vals
      q <- n
    } else {
      v <- g_vals
      q <- m
    }

    quad <- sum(backsolve(L, v, transpose = TRUE)^2)
    out_log[i] <- -0.5 * (q * base::log(2 * base::pi) + log_det_sigma + quad)
  }

  if (log) out_log else base::exp(out_log)
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
