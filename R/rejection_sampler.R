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
#'
#' @return A matrix or tibble containing the accepted samples.
#'
#' @examples
#' \dontrun{
#' library("mpoly")
#'
#' # Single polynomial (circle)
#' p1 <- mp("x^2 + y^2 - 1")
#' set.seed(1)
#' rejection_sampler(100, p1, sd = 0.05, w = 1.5)
#'
#' # Uniform band proposal around the variety, returning a tibble
#' rejection_sampler(
#'   100, p1, sd = 0.05, w = c(-1.5, 1.5),
#'   dist = "unif", output = "tibble"
#' )
#'
#' # Two-polynomial system (upper/lower acceptance geometry differs by `homo`)
#' p2 <- mp(c("x^2 + y^2 - 1", "y"))
#' rejection_sampler(50, p2, sd = 0.05, w = 1.5, homo = TRUE)
#' rejection_sampler(50, p2, sd = c(0.05, 0.05), w = 1.5, homo = FALSE)
#' }



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
                              message = FALSE) {
  # expand scalar/length-2 window specs to per-variable bounds
  n_vars <- length(vars)
  if (is.numeric(w) && length(w) == 1) {
    w <- replicate(length(vars), c(-w, w), simplify = FALSE)
    names(w) <- vars
  } else if (is.numeric(w) && length(w) == 2) {
    w <- replicate(length(vars), w, simplify = FALSE)
    names(w) <- vars
  }
  dist <- match.arg(dist)
  if (correct_p_coefficients) poly <- normalize_coefficients(poly)

  if (is.mpolyList(poly)) {
    # mpolyList case: g(x) is vector-valued and uses Jacobian-based scaling
    n_polys <- length(poly)
    pf <- as.function(poly, varorder = vars, silent = TRUE)
    dp <- dpfs <- vector(mode = "list", length = n_polys)
    for (i in seq_len(n_polys)) {
      dp[[i]] <- deriv(poly[[i]], vars)
      if (correct_dp_coefficients) dp[[i]] <- normalize_coefficients(dp[[i]])
      dpfs[[i]] <- as.function(dp[[i]], varorder = vars, silent = TRUE)
    }

    # assemble Jacobian row-by-row from per-polynomial gradients
    dpf <- function(x) {
      mat <- matrix(NA_real_, nrow = n_polys, ncol = n_vars)
      for (i in seq_len(n_polys)) mat[i, ] <- dpfs[[i]](x)
      mat
    }

    # build the acceptance kernel based on sd shape and homo setting
    if (is.vector(sd)) {
      # scalar/diagonal scale case
      if (homo) {
        if (length(sd) == 1) {
          pbar <- function(x) {
            g <- pf(x)
            J <- dpf(x)
            as.numeric(g %*% solve(tcrossprod(J), g))
          }
          log_ptilde <- function(x) -pbar(x) / (2 * sd^2)
          ptilde <- function(x) exp(-pbar(x) / (2 * sd^2))
        } else {
          if (length(sd) != n_vars) {
            stop(
              "When `poly` is an `mpolyList` and `homo = TRUE`, ",
              "length(`sd`) must be 1, `length(vars)`, or a covariance matrix.",
              call. = FALSE
            )
          }
          pbar <- function(x) {
            g <- pf(x)
            J <- dpf(x)
            if (n_vars == n_polys) {
              J_inv <- solve(J)
            } else {
              J_inv <- MASS::ginv(J)
            }
            v <- (J_inv %*% g) / sd
            as.numeric(crossprod(v))
          }
          log_ptilde <- function(x) -pbar(x) / 2
          ptilde <- function(x) exp(-pbar(x) / 2)
        }
      } else {
        if (length(sd) == 1) {
          pbar <- function(x) {
            g <- as.numeric(pf(x))
            sum(g^2)
          }
          log_ptilde <- function(x) -pbar(x) / (2 * sd^2)
          ptilde <- function(x) exp(-pbar(x) / (2 * sd^2))
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
          ptilde <- function(x) exp(-pbar(x) / 2)
        }
      }
    } else if (is.matrix(sd)) {
      # full covariance case via spectral decomposition
      eig_sd <- eigen(sd, symmetric = TRUE)
      la <- eig_sd$values
      P <- eig_sd$vectors
      la_inv <- 1 / la
      sqrt_sd_inv <- diag(sqrt(la_inv)) %*% t(P)

      if (homo) {
        pbar <- function(x) {
          g <- pf(x)
          J <- dpf(x)
          if (n_vars == n_polys) {
            J_inv <- solve(J)
          } else {
            J_inv <- MASS::ginv(J)
          }
          as.numeric(crossprod(sqrt_sd_inv %*% J_inv %*% g))
        }
      } else {
        pbar <- function(x) {
          g <- as.numeric(pf(x))
          as.numeric(crossprod(sqrt_sd_inv %*% g))
        }
      }

      log_ptilde <- function(x) -pbar(x) / 2
      ptilde <- function(x) exp(-pbar(x) / 2)
    }
  } else if (is.mpoly(poly)) {
    # single-polynomial case
    pf <- as.function(poly, varorder = vars, silent = TRUE)
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
      as.function(ssdp, varorder = vars, silent = TRUE)
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
    ptilde <- function(x) exp(-pbar(x)^2 / (2 * sd^2))
  } else {
    stop("`poly` should be a `mpoly` or `mpolyList` object.", call. = FALSE)
  }

  mat <- matrix(nrow = 0, ncol = n_vars, dimnames = list(NULL, vars))
  n_remaining <- n

  # rejection loop: propose uniformly in box, accept by kernel
  while (n_remaining > 0) {
    if (message) cat("\r", strrep(" ", 80))
    if (message) {
      cat(
        "\r",
        scales::number_format(big.mark = ",")(n_remaining),
        " remaining...",
        sep = ""
      )
    }

    u <- matrix(nrow = n_remaining, ncol = n_vars, dimnames = list(NULL, vars))
    for (i in seq_len(n_vars)) {
      u[, i] <- runif(n_remaining, w[[i]][1], w[[i]][2])
    }

    if (dist == "norm") {
      gbars <- apply(u, 1, log_ptilde)
      accept_reject_ndcs <- which(log(runif(n_remaining)) <= gbars)
    } else if (dist == "unif") {
      gbars <- apply(u, 1, pbar)
      accept_reject_ndcs <- which(abs(gbars) <= sd)
    }
    n_accept <- length(accept_reject_ndcs)

    if (n_accept > 0) {
      samp <- u[accept_reject_ndcs, , drop = FALSE]
      mat <- rbind(mat, samp)
      n_remaining <- n_remaining - n_accept
    }
  }
  if (message) cat("\r", strrep(" ", 80))
  if (message) cat("\r")
  row.names(mat) <- NULL
  out <- if (output == "tibble") tibble::as_tibble(mat) else mat
  out

}
