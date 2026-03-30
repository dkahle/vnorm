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
#' @param sigma For the single-polynomial case, a positive scalar standard
#'   deviation. For the multi-polynomial case, a scalar, vector, or matrix. If
#'   `homo = TRUE`, `sigma` must conform to the number of variables; if
#'   `homo = FALSE`, it must conform to the number of polynomials.
#' @param homo Logical; default is `TRUE`. If `TRUE`, compute the homoskedastic
#'   variety normal pseudo-density. If `FALSE`, compute the heteroskedastic
#'   pseudo-density.
#' @param log Logical. If `TRUE`, returns the log of the density.
#' @return A numeric scalar or vector containing the pseudo-density evaluated at
#'   `x`.
#'
#' @examples
#'
#' library("mpoly")
#'
#' ## Single polynomial in one variable
#' p1 <- mp("x")
#' pdvnorm(c(-1, 0, 1), p1, sigma = 1)
#' pdvnorm(0, p1, sigma = 2, log = TRUE)
#'
#' ## Multivariate (square) system: two polynomials in two variables
#' p2 <- mp(c("x", "y"))
#' x2 <- rbind(c(0, 0), c(1, 2), c(-1, 3))
#'
#' ## Different sigma forms
#' pdvnorm(x2, p2, sigma = 1)
#' pdvnorm(x2, p2, sigma = c(1, 2))
#' pdvnorm(x2, p2, sigma = diag(c(1, 4)))
#'
#' ## Multivariate (underdetermined): one polynomial in two variables
#' p3 <- mp("x + y")
#' x3 <- rbind(c(1, 1), c(2, -1), c(0, 3))
#' pdvnorm(x3, p3, sigma = 1)
#' pdvnorm(as.data.frame(x3), p3, sigma = 1)
#'
#' ## Multivariate (overdetermined): three polynomials in two variables
#' p4 <- mp(c("x", "y", "x + y"))
#' x4 <- rbind(c(1, 2), c(0, -1), c(2, 2))
#' pdvnorm(x4, p4, sigma = diag(2),    homo = TRUE)
#' pdvnorm(x4, p4, sigma = c(1, 2, 3), homo = FALSE)
#'
#' @export
pdvnorm <- function(x, poly, sigma, homo = TRUE, log = FALSE) {
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

  if (is_uni) {
    # univariate/single-polynomial path
    if (!is.numeric(sigma) || length(sigma) != 1 || sigma <= 0) {
      stop(
        "For the single-polynomial case, 'sigma' must be a single positive numeric (sd)."
      )
    }
    sd <- as.numeric(sigma)

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
  if (!is.numeric(sigma) || any(!is.finite(sigma))) {
    stop("'sigma' must be finite numeric.")
  }

  if (length(sigma) == 1) {
    if (sigma <= 0) {
      stop("'sigma' must be positive.")
    }
    Sigma <- if (homo) diag(sigma, n) else diag(sigma, m)
  } else if (is.vector(sigma)) {
    if (any(sigma <= 0)) {
      stop("'sigma' vector entries must be positive.")
    }
    if (homo && length(sigma) != n) {
      stop("When homo=TRUE, length(sigma) must be n.")
    }
    if (!homo && length(sigma) != m) {
      stop("When homo=FALSE, length(sigma) must be m.")
    }
    Sigma <- diag(sigma)
  } else {
    Sigma <- as.matrix(sigma)
    if (homo) {
      if (!all(dim(Sigma) == c(n, n))) {
        stop("'sigma' matrix must be n x n when homo=TRUE.")
      }
    } else {
      if (!all(dim(Sigma) == c(m, m))) {
        stop("'sigma' matrix must be m x m when homo=FALSE.")
      }
    }
  }

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

      # force() captures the current iteration's constant, not the final one
      force(const)
      grad_fun[[i]] <- function(...) {
        const
      }
    } else {
      grad_fun[[i]] <- suppressMessages(
        as.function(deriv(poly[[i]], var = vars), varorder = vars)
      )
    }
  }

  # Cholesky factor for efficient Mahalanobis distance
  L <- tryCatch(
    chol(Sigma),
    error = function(e) {
      stop("'sigma' must be positive definite.", call. = FALSE)
    }
  )
  log_det_sigma <- 2 * sum(base::log(diag(L)))

  out_log <- numeric(nrow(X))
  for (i in seq_len(nrow(X))) {
    xi <- as.numeric(X[i, ])
    g_vals <- as.numeric(g_vals_mat[i, ])

    if (homo) {
      # Jacobian of g(x): rows are equations, columns are variables
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
