#' Rejection Sampler for the Variety Normal Distribution
#'
#' Perform rejection sampling to generate draws from a
#' variety normal distribution.
#'
#' @param n The number of accepted draws to return.
#' @param poly An `mpoly` object or an `mpolyList` object.
#' @param sd The "standard deviation" component of the normal kernel.
#' @param vars A character vector of the indeterminates in the distribution.
#' @param w Proposal box constraints. If a single number, a box window
#'   `(-w, w)` is applied to all variables. If length 2, the same interval is
#'   used for all variables. A named list can be used to specify bounds for
#'   each variable.
#' @param output Either `"simple"` or `"tibble"` output format.
#' @param dist Either `"norm"` (normal) or `"unif"` (uniform).
#' @param homo If `TRUE`, sampling is done from a homoskedastic variety normal
#'   distribution.
#' @param correct_p_coefficients If `TRUE`, normalize polynomial coefficients.
#' @param correct_dp_coefficients If `TRUE`, normalize derivative coefficients.
#' @param message If `TRUE`, print progress messages showing remaining samples.
#' @param max_batches Maximum number of proposal batches before stopping with
#'   an informative error.
#'
#' @return A matrix or tibble containing the accepted samples.
#'
#' @examples
#' \dontrun{
#' library("mpoly")
#'
#' # single polynomial (circle)
#' p1 <- mp("x^2 + y^2 - 1")
#' set.seed(1)
#' rejection_sampler(100, p1, sd = 0.05, w = 1.5)
#'
#'
#'
#' # uniform band proposal around the variety, returning a tibble
#' rejection_sampler(
#'   100, p1, sd = 0.05, w = c(-1.5, 1.5),
#'   dist = "unif", output = "tibble"
#' )
#'
#'
#'
#' # two-polynomial system (upper/lower acceptance geometry differs by `homo`)
#' p2 <- mp(c("x^2 + y^2 - 1", "y"))
#' rejection_sampler(50, p2, sd = 0.05, w = 1.5, homo = TRUE)
#' rejection_sampler(50, p2, sd = c(0.05, 0.05), w = 1.5, homo = FALSE)
#' }
#'
#' @export
rejection_sampler <- function(n,
                              poly,
                              sd = .01,
                              vars = sort(mpoly::vars(poly)),
                              w = 1.25,
                              output = "simple",
                              dist = c("norm", "unif"),
                              homo = TRUE,
                              correct_p_coefficients = FALSE,
                              correct_dp_coefficients = FALSE,
                              message = FALSE,
                              max_batches = 1000L) {
  if (
    !is.numeric(n) || length(n) != 1L || !is.finite(n) ||
      n < 1 || n != as.integer(n)
  ) {
    stop("`n` must be a positive integer.", call. = FALSE)
  }
  n <- as.integer(n)
  if (!(is.mpoly(poly) || is.mpolyList(poly))) {
    stop("`poly` should be a `mpoly` or `mpolyList` object.", call. = FALSE)
  }
  if (length(mpoly::vars(poly)) == 0L) {
    stop("`poly` must contain at least one variable.", call. = FALSE)
  }
  if (!is.character(vars) || any(!nzchar(vars)) || anyDuplicated(vars)) {
    stop("`vars` must be a character vector of unique variable names.", call. = FALSE)
  }
  output <- match.arg(output, c("simple", "tibble"))
  dist <- match.arg(dist)
  if (!is.logical(homo) || length(homo) != 1L || is.na(homo)) {
    stop("`homo` must be TRUE or FALSE.", call. = FALSE)
  }
  if (
    !is.logical(correct_p_coefficients) || length(correct_p_coefficients) != 1L ||
      is.na(correct_p_coefficients)
  ) {
    stop("`correct_p_coefficients` must be TRUE or FALSE.", call. = FALSE)
  }
  if (
    !is.logical(correct_dp_coefficients) || length(correct_dp_coefficients) != 1L ||
      is.na(correct_dp_coefficients)
  ) {
    stop("`correct_dp_coefficients` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(message) || length(message) != 1L || is.na(message)) {
    stop("`message` must be TRUE or FALSE.", call. = FALSE)
  }
  if (
    !is.numeric(max_batches) || length(max_batches) != 1L ||
      !is.finite(max_batches) || max_batches < 1 ||
      max_batches != as.integer(max_batches)
  ) {
    stop("`max_batches` must be a positive integer.", call. = FALSE)
  }
  max_batches <- as.integer(max_batches)

  n_vars <- length(vars)
  w <- rejection_normalize_window(w, vars)
  sd <- rejection_validate_scale(sd, poly, homo, dist, n_vars)
  if (correct_p_coefficients) poly <- normalize_coefficients(poly)

  if (is.mpolyList(poly)) {
    # mpolyList case: g(x) is vector-valued and uses jacobian-based scaling
    n_polys <- length(poly)
    pf <- poly_function(poly, varorder = vars)
    dp <- dpfs <- vector(mode = "list", length = n_polys)
    for (i in seq_len(n_polys)) {
      dp[[i]] <- gradient_mpolylist(poly[[i]], vars)
      if (correct_dp_coefficients) dp[[i]] <- normalize_coefficients(dp[[i]])
      dpfs[[i]] <- poly_function(dp[[i]], varorder = vars)
    }

    # assemble jacobian row-by-row from per-polynomial gradients
    dpf <- function(x) {
      mat <- matrix(NA_real_, nrow = n_polys, ncol = n_vars)
      for (i in seq_len(n_polys)) mat[i, ] <- dpfs[[i]](x)
      mat
    }

    # J^+ g, the first-order displacement to the variety measured by every
    # homoskedastic kernel, with pdvnorm()'s conventions: least squares for
    # overdetermined systems, and no density (an infinite displacement) where
    # a square or underdetermined system's J loses rank and g leaves its
    # column space
    newton_residual <- function(x) {
      J <- dpf(x)
      g <- as.numeric(pf(x))
      J_inv <- pdvnorm_pseudoinverse(J)
      if (
        n_polys <= n_vars && attr(J_inv, "rank") < n_polys &&
          pdvnorm_outside_column_space(J, matrix(g, nrow = 1))
      ) {
        return(rep(Inf, n_vars))
      }
      as.numeric(J_inv %*% g)
    }
    # for dependent equations, whether g stays in J's column space is the
    # same at generic points, so the kernel's own rule at a few points of the
    # window shows whether any draw could be accepted
    if (homo && n_polys <= n_vars && dependent_equations(poly, vars)) {
      u <- matrix((seq_len(3L * n_vars) * 0.6180339887498949) %% 1, nrow = 3L)
      probes <- vapply(
        seq_len(n_vars), function(j) w[[j]][1] + u[, j] * diff(w[[j]]), numeric(3L)
      )
      vanishes <- vapply(seq_len(3L), function(i) {
        is.infinite(newton_residual(probes[i, ])[1L])
      }, logical(1))
      if (all(vanishes)) {
        stop(
          "The density is 0 off the variety: the Jacobian of `poly` loses rank ",
          "everywhere (its equations are dependent) and the residual leaves its ",
          "column space. Remove the dependent equations.",
          call. = FALSE
        )
      }
    }

    # build the acceptance kernel based on sd shape and homo setting
    if (is.vector(sd)) {
      # scalar/diagonal scale case
      if (homo) {
        if (length(sd) == 1) {
          pbar <- function(x) sum(newton_residual(x)^2)
          log_ptilde <- function(x) -pbar(x) / (2 * sd^2)
        } else {
          if (length(sd) != n_vars) {
            stop(
              "When `poly` is an `mpolyList` and `homo = TRUE`, ",
              "length(`sd`) must be 1, `length(vars)`, or a covariance matrix.",
              call. = FALSE
            )
          }
          pbar <- function(x) sum((newton_residual(x) / sd)^2)
          log_ptilde <- function(x) -pbar(x) / 2
        }
      } else {
        if (length(sd) == 1) {
          pbar <- function(x) {
            g <- as.numeric(pf(x))
            sum(g^2)
          }
          log_ptilde <- function(x) -pbar(x) / (2 * sd^2)
        } else {
          if (length(sd) != n_polys) {
            stop(
              "When `poly` is an `mpolyList` and `homo = FALSE`, ",
              "length(`sd`) must be 1, `length(poly)`, or a covariance matrix.",
              call. = FALSE
            )
          }
          pbar <- function(x) {
            g <- as.numeric(pf(x))
            sum((g / sd)^2)
          }
          log_ptilde <- function(x) -pbar(x) / 2
        }
      }
    } else if (is.matrix(sd)) {
      # full covariance case via spectral decomposition
      eig_sd <- eigen(sd, symmetric = TRUE)
      la <- eig_sd$values
      P <- eig_sd$vectors
      la_inv <- 1 / la
      # nrow keeps diag() from reading a 1 x 1 scale as an identity size
      sqrt_sd_inv <- diag(sqrt(la_inv), nrow = length(la_inv)) %*% t(P)

      if (homo) {
        pbar <- function(x) {
          v <- newton_residual(x)
          if (any(is.infinite(v))) return(Inf)
          as.numeric(crossprod(sqrt_sd_inv %*% v))
        }
      } else {
        pbar <- function(x) {
          g <- as.numeric(pf(x))
          as.numeric(crossprod(sqrt_sd_inv %*% g))
        }
      }

      log_ptilde <- function(x) -pbar(x) / 2
    }
  } else if (is.mpoly(poly)) {
    # single-polynomial case
    pf <- poly_function(poly, varorder = vars)
    dp <- deriv(poly, vars)
    if (correct_dp_coefficients) dp <- normalize_coefficients(dp)

    # sum of squared partial derivatives (gradient norm squared)
    if (n_vars > 1) {
      ssdp <- mp("0")
      for (i in seq_len(n_vars)) {
        ssdp <- ssdp + dp[[i]]^2
      }
    } else {
      ssdp <- dp^2
    }

    ssdpf <- if (is.constant(ssdp)) {
      function(x) ssdp[[1]][["coef"]]
    } else {
      poly_function(ssdp, varorder = vars)
    }

    if (homo) {
      # normalize by gradient magnitude for approximate arc-length scaling
      pbar <- function(x) pf(x) / sqrt(ssdpf(x))

      # log-space arithmetic avoids overflow of pbar(x)^2
      log_ptilde <- function(x) {
        -exp((2 * log(abs(pf(x))) - log(ssdpf(x))) - log(2) - 2 * log(sd))
      }
    } else {
      pbar <- pf
      log_ptilde <- function(x) {
        -exp((2 * log(abs(pf(x)))) - log(2) - 2 * log(sd))
      }
    }
  } else {
    stop("`poly` should be a `mpoly` or `mpolyList` object.", call. = FALSE)
  }

  mat <- matrix(nrow = 0, ncol = n_vars, dimnames = list(NULL, vars))
  n_remaining <- n
  n_proposed <- 0
  n_accepted <- 0
  batch_size <- n
  batch <- 0L

  # rejection loop: propose uniformly in box, accept by kernel
  while (n_remaining > 0) {
    batch <- batch + 1L
    if (batch > max_batches) {
      if (message) cat("\r", strrep(" ", 80), "\r", sep = "")
      stop(
        "rejection_sampler() accepted ",
        n - n_remaining,
        " of ",
        n,
        " requested draws after ",
        max_batches,
        " proposal batches; try enlarging `w`, increasing `sd`, ",
        "or increasing `max_batches`.",
        call. = FALSE
      )
    }

    if (message) cat("\r", strrep(" ", 80))
    if (message) {
      cat(
        "\r",
        scales::number_format(big.mark = ",")(n_remaining),
        " remaining...",
        sep = ""
      )
    }

    u <- matrix(nrow = batch_size, ncol = n_vars, dimnames = list(NULL, vars))
    for (i in seq_len(n_vars)) {
      u[, i] <- runif(batch_size, w[[i]][1], w[[i]][2])
    }

    if (dist == "norm") {
      gbars <- apply(u, 1, log_ptilde)
      accept_reject_ndcs <- which(log(runif(batch_size)) <= gbars)
    } else if (dist == "unif") {
      gbars <- apply(u, 1, pbar)
      accept_reject_ndcs <- which(abs(gbars) <= sd)
    }
    n_accept <- length(accept_reject_ndcs)
    n_proposed <- n_proposed + batch_size
    n_accepted <- n_accepted + n_accept

    if (n_accept > 0) {
      # accepted proposals are iid draws, so keeping the first ones still
      # needed does not bias the sample
      keep <- accept_reject_ndcs[seq_len(min(n_accept, n_remaining))]
      mat <- rbind(mat, u[keep, , drop = FALSE])
      n_remaining <- n_remaining - length(keep)
    }

    # size the next batch from the running acceptance rate; shrinking batches
    # to n_remaining would need about log(n) / rate batches at low rates
    batch_size <- if (n_accepted > 0) {
      max(n_remaining, min(ceiling(1.2 * n_remaining * n_proposed / n_accepted), 1e5))
    } else {
      n_remaining
    }
  }
  if (message) cat("\r", strrep(" ", 80))
  if (message) cat("\r")
  row.names(mat) <- NULL
  out <- if (output == "tibble") tibble::as_tibble(mat) else mat
  out

}

rejection_normalize_window <- function(w, vars) {
  check_bounds <- function(bounds, label) {
    if (
      !is.numeric(bounds) || length(bounds) != 2L ||
        any(!is.finite(bounds)) || bounds[1] >= bounds[2]
    ) {
      stop(label, " must be a finite numeric interval `c(lower, upper)`.", call. = FALSE)
    }
    as.numeric(bounds)
  }

  if (is.numeric(w) && length(w) == 1L) {
    if (!is.finite(w) || w <= 0) {
      stop("Scalar `w` must be positive and finite.", call. = FALSE)
    }
    out <- replicate(length(vars), c(-w, w), simplify = FALSE)
    names(out) <- vars
    return(out)
  }
  if (is.numeric(w) && length(w) == 2L) {
    bounds <- check_bounds(w, "`w`")
    out <- replicate(length(vars), bounds, simplify = FALSE)
    names(out) <- vars
    return(out)
  }
  if (!is.list(w) || is.null(names(w))) {
    stop("`w` must be a positive scalar, a length-2 interval, or a named list.", call. = FALSE)
  }
  missing_vars <- setdiff(vars, names(w))
  if (length(missing_vars) > 0L) {
    stop(
      "`w` is missing bounds for: ",
      paste(missing_vars, collapse = ", "),
      call. = FALSE
    )
  }
  out <- lapply(vars, function(var) check_bounds(w[[var]], paste0("`w$", var, "`")))
  names(out) <- vars
  out
}

rejection_validate_scale <- function(sd, poly, homo, dist, n_vars) {
  if (!is.numeric(sd) || any(!is.finite(sd))) {
    stop("`sd` must be finite numeric.", call. = FALSE)
  }
  if (dist == "unif" && !(is.null(dim(sd)) && length(sd) == 1L)) {
    stop("`dist = \"unif\"` requires scalar `sd`.", call. = FALSE)
  }

  if (is.mpoly(poly)) {
    if (!is.null(dim(sd)) || length(sd) != 1L) {
      stop("When `poly` is an `mpoly`, `sd` must be a positive scalar.", call. = FALSE)
    }
    if (sd <= 0) stop("`sd` must be positive.", call. = FALSE)
    return(as.numeric(sd))
  }

  n_polys <- length(poly)
  target_dim <- if (homo) n_vars else n_polys
  target_label <- if (homo) "`length(vars)`" else "`length(poly)`"

  if (is.matrix(sd)) {
    if (!all(dim(sd) == c(target_dim, target_dim))) {
      stop(
        "When `poly` is an `mpolyList`, matrix `sd` must be ",
        target_label,
        " by ",
        target_label,
        ".",
        call. = FALSE
      )
    }
    tryCatch(
      chol(sd),
      error = function(e) {
        stop("Matrix `sd` must be positive definite.", call. = FALSE)
      }
    )
    return(sd)
  }

  if (any(sd <= 0)) stop("`sd` must be positive.", call. = FALSE)
  if (length(sd) == 1L || length(sd) == target_dim) {
    return(as.numeric(sd))
  }

  stop(
    "When `poly` is an `mpolyList` and `homo = ",
    homo,
    "`, length(`sd`) must be 1 or ",
    target_label,
    ".",
    call. = FALSE
  )
}
