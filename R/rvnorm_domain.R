# domain-noise draws for rvnorm(noise = "domain"): points approximately
# uniform on the variety, plus multivariate normal noise in the space of the
# variables. Split out from rvnorm.R, which dispatches here

rvnorm_domain <- function(
    n, poly, sd, Sigma, pilot_sd, output, rejection, seed, vars, w,
    code_only, pilot_args
  ) {
  if (!(is.mpoly(poly) || is.mpolyList(poly))) {
    stop(
      "`poly` should be either a character vector, mpoly, or mpolyList.",
      call. = FALSE
    )
  }
  if (length(mpoly::vars(poly)) == 0L) {
    stop("`poly` must contain at least one variable.", call. = FALSE)
  }
  sample_vars <- rvnorm_resolve_vars(vars, poly)
  noise_cov <- rvnorm_domain_covariance(sd, Sigma, sample_vars)
  if (is.null(pilot_sd)) {
    # the projected pilot draws are uniform on the variety up to terms of
    # order pilot_sd^2; a tenth of the noise scale makes those negligible
    pilot_sd <- sqrt(mean(diag(noise_cov))) / 10
  }
  proj_vars <- paste0(sample_vars, "_proj")
  if (output == "tibble" && any(proj_vars %in% sample_vars)) {
    stop(
      "Variable names clash with the `<var>_proj` columns of the points on ",
      "the variety: ", paste(intersect(proj_vars, sample_vars), collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  # the window W bounds the points on the variety. The pilot runs on W
  # widened by 6 pilot_sd, so that the normal lines of points near W's edge
  # are not cut short, which would thin them out; points that land outside W
  # are dropped. rvnorm(rejection = TRUE) uses w = 1.25 when none is given
  if (is.null(w) && rejection) w <- 1.25
  window <- if (!is.null(w)) rvnorm_normalize_window(w, sample_vars)
  if (!is.null(window)) {
    pilot_args$w <- rvnorm_widen_window(
      window, sample_vars, mpoly::vars(poly), 6 * pilot_sd
    )
  }

  if (code_only && !rejection) {
    return(rvnorm_domain_pilot(n, poly, pilot_sd, rejection, seed, pilot_args))
  }

  # the latent points must lie on the variety to within 1% of the noise's
  # smallest sd, or they would add noise of their own
  noise_sd_min <- sqrt(min(eigen(noise_cov, symmetric = TRUE, only.values = TRUE)$values))
  draws <- rvnorm_with_seed(seed, {
    latent <- rvnorm_domain_latent(
      n, poly, sample_vars, pilot_sd, rejection, seed, pilot_args, window,
      noise_sd_min
    )
    latent$x <- latent$points + rvnorm_gaussian_noise(n, noise_cov)
    latent
  })

  x <- draws$x
  colnames(x) <- sample_vars
  if (output == "simple") {
    return(as.data.frame(x))
  }
  latent <- draws$points
  colnames(latent) <- proj_vars
  out <- cbind(as.data.frame(x), as.data.frame(latent))
  if (!is.null(draws$index)) {
    out$.chain <- as.integer(draws$index[, ".chain"])
    out$.iteration <- as.integer(draws$index[, ".iteration"])
    out$.draw <- seq_len(n)
  }
  tibble::as_tibble(out)
}

rvnorm_domain_latent <- function(
    n, poly, vars, pilot_sd, rejection, seed, pilot_args, window,
    noise_sd_min = Inf
  ) {
  # n points nearly uniform on the variety, within the window if there is
  # one: pilot draws of the homogeneous variety normal, each moved to the
  # closest point on the variety. Each point must be located to within 1%
  # of noise_sd_min. Returns the points, the chain and iteration of their
  # pilot draws (NULL for rejection pilots), and how many projections were
  # tried and failed
  points <- matrix(numeric(0), nrow = 0, ncol = length(vars))
  index <- NULL
  chains_used <- 0L
  degree <- max(mpoly::totaldeg(poly))
  accuracy <- 0.01 * noise_sd_min
  n_projected <- 0L
  n_failed <- 0L
  n_off <- 0L
  n_far <- 0L
  n_accepted <- 0L
  n_stretched <- 0L
  n_imprecise <- 0L
  precision_missed <- numeric(0)
  # conditions under which the latent points cannot be trusted, checked
  # after every pilot run so as not to waste further runs
  stop_if_untrusted <- function() {
    if (n_stretched > 0.01 * n_accepted) {
      stop(
        "At ", n_stretched, " of ", n_accepted, " pilot draws, the ",
        "residual |J^+ g| was under two thirds of the distance to the ",
        "variety, so the pilot would weight those parts of the variety ",
        "more. That happens where an equation of `poly` has a repeated ",
        "factor, such as (x - 1)^2, or where `pilot_sd` is large next to ",
        "the variety's curvature. Use equations without repeated factors, ",
        "or lower `pilot_sd`.",
        call. = FALSE
      )
    }
    if (n_imprecise > 0.01 * (n_accepted + n_imprecise)) {
      stop(
        "At ", n_imprecise, " of ", n_accepted + n_imprecise, " pilot ",
        "draws, the projection locates the variety only to about ",
        signif(stats::median(precision_missed), 2), ", not within 1% of the ",
        "noise's smallest sd, ", signif(noise_sd_min, 2), ". Near a ",
        "repeated factor, such as (x - 1)^3, the expanded polynomial loses ",
        "precision. Use the polynomial without the repeated factor, or ",
        "larger noise.",
        call. = FALSE
      )
    }
    if (n_failed > 0.05 * n_projected) {
      stop(
        n_failed, " of ", n_projected, " projections onto the variety did ",
        "not converge: Gauss-Newton ran out of iterations, or the Jacobian ",
        "of `poly` lost rank. Parts of the variety would be missing. Check ",
        "it for singular components, or lower `pilot_sd`.",
        call. = FALSE
      )
    }
  }
  # spare draws for projections that fail or leave the window
  n_draw <- n + ceiling(if (is.null(window)) 0.02 * n else 0.1 * n) + 10L
  for (round in seq_len(5L)) {
    pilot <- as.data.frame(rvnorm_domain_pilot(
      n_draw, poly, pilot_sd, rejection,
      if (!is.null(seed)) seed + round - 1L, pilot_args
    ))
    x0 <- as.matrix(pilot[, vars, drop = FALSE])
    proj <- project_rows_closest(x0, poly, vars, scale = pilot_sd)
    # the pilot's residual |J^+ g| is the distance to the variety to first
    # order where poly is square-free, and about the distance / k on a
    # component where it has a factor of multiplicity k, like (x - 1)^2.
    # There the pilot's tube is k times wider, and its points weigh k
    # times as much. Crossings of two smooth branches keep the ratio below
    # sqrt(2)
    distance <- sqrt(rowSums((proj$points - x0)^2))
    multiplicity <- ifelse(
      proj$ok & proj$residual > 0, distance / proj$residual, 1
    )
    if (is.mpoly(poly)) {
      # a factor's multiplicity is at most the degree. Draws from a
      # component's tube lie within about 5 pilot_sd k of it; a projection
      # farther from its draw, or with a larger ratio, jumped to another
      # branch
      k <- pmax(1, round(multiplicity))
      near <- k <= degree & distance <= 10 * pilot_sd * k
    } else {
      # a system's correction is a product over the normal directions,
      # which one draw cannot estimate
      k <- rep(1, nrow(x0))
      near <- distance <= 10 * pilot_sd
      n_stretched <- n_stretched + sum(proj$ok & near & multiplicity > 1.5)
    }
    precise <- proj$ok & near & proj$precision <= accuracy
    n_imprecise <- n_imprecise + sum(proj$ok & near & !precise)
    precision_missed <- c(precision_missed, proj$precision[proj$ok & near & !precise])
    ok <- precise
    keep <- ok & rvnorm_in_window(proj$points, window, slack = proj$precision)
    # keep a draw on a component of multiplicity k with probability 1 / k,
    # which restores each component's share
    if (any(keep & k > 1)) keep <- keep & stats::runif(length(k)) < 1 / k
    # projections that did not converge lose points of the variety, which
    # must not pass silently; those that converged off the variety (at an
    # overdetermined system's least-squares points) or far from their draw
    # came from pilot draws that were not near it
    n_projected <- n_projected + nrow(x0)
    n_failed <- n_failed + sum(!proj$converged)
    n_off <- n_off + sum(proj$converged & !proj$ok)
    n_far <- n_far + sum(proj$ok & !near)
    n_accepted <- n_accepted + sum(ok)
    points <- rbind(points, proj$points[keep, , drop = FALSE])
    if (all(c(".chain", ".iteration") %in% names(pilot))) {
      # each run's chains follow the previous runs', so that a chain and
      # iteration identify one pilot draw
      run_index <- as.matrix(pilot[keep, c(".chain", ".iteration")])
      run_index[, ".chain"] <- run_index[, ".chain"] + chains_used
      chains_used <- chains_used + max(pilot$.chain)
      index <- rbind(index, run_index)
    }
    stop_if_untrusted()
    if (nrow(points) >= n) break
    rate <- max(mean(keep), 0.05)
    n_draw <- ceiling(1.5 * (n - nrow(points)) / rate) + 10L
  }
  if (nrow(points) < n) {
    stop(
      "Only ", nrow(points), " of ", n, " pilot draws projected onto the ",
      "variety", if (!is.null(window)) " inside the window `w`", ". Check ",
      "that the variety has points there and that `poly` has no repeated ",
      "factors, or lower `pilot_sd`.",
      call. = FALSE
    )
  }
  n_dropped <- n_failed + n_off + n_far + n_imprecise
  if (n_dropped > 0.01 * n_projected) {
    warning(
      n_dropped, " of ", n_projected, " pilot draws were dropped: ",
      n_failed, " projections did not converge, ", n_off, " converged off ",
      "the variety (as at an overdetermined system's least-squares points), ",
      n_far, " landed far from their draw, and ", n_imprecise, " were not ",
      "precise enough.",
      call. = FALSE
    )
  }
  # an evenly spaced subset keeps each chain's share of the draws
  rows <- ceiling(seq_len(n) * nrow(points) / n)
  list(
    points = points[rows, , drop = FALSE],
    index = if (!is.null(index)) index[rows, , drop = FALSE],
    n_projected = n_projected,
    n_failed = n_failed
  )
}

rvnorm_gaussian_noise <- function(n, Sigma) {
  # n draws of N(0, Sigma), one per row: rows of z %*% chol(Sigma) have
  # covariance t(chol(Sigma)) %*% chol(Sigma)
  matrix(stats::rnorm(n * ncol(Sigma)), n, ncol(Sigma)) %*% chol(Sigma)
}

rvnorm_domain_pilot <- function(n, poly, pilot_sd, rejection, seed, pilot_args) {
  # draws of the homogeneous variety normal with scale pilot_sd, whose
  # projections onto the variety are its nearly uniform points
  do.call(rvnorm, c(
    list(
      n = n, poly = poly, sd = pilot_sd, output = "tibble",
      rejection = rejection, homo = TRUE, noise = "codomain", seed = seed
    ),
    pilot_args
  ))
}

rvnorm_domain_covariance <- function(sd, Sigma, vars) {
  # the covariance of the noise in the space of the sampled variables, in the
  # order of vars; named entries are matched to variables by name
  n_vars <- length(vars)
  by_name <- function(v, label) {
    if (is.null(names(v))) return(v)
    if (anyDuplicated(names(v)) || !setequal(names(v), vars)) {
      stop(
        "Names of `", label, "` must match the variables: ",
        paste(vars, collapse = ", "), ".",
        call. = FALSE
      )
    }
    v[vars]
  }
  if (!is.null(Sigma)) {
    if (!is.numeric(Sigma) || any(!is.finite(Sigma))) {
      stop("`Sigma` must be finite numeric.", call. = FALSE)
    }
    if (is.null(dim(Sigma))) {
      rvnorm_check_positive_numeric(Sigma, "Sigma")
      if (length(Sigma) == 1L) Sigma <- rep(as.numeric(Sigma), n_vars)
      if (length(Sigma) != n_vars) {
        stop("Length of `Sigma` must match the number of variables.", call. = FALSE)
      }
      return(diag(as.numeric(by_name(Sigma, "Sigma")), nrow = n_vars))
    }
    Sigma <- as.matrix(Sigma)
    if (!all(dim(Sigma) == c(n_vars, n_vars))) {
      stop(
        "`Sigma` matrix dimensions must match the number of variables.",
        call. = FALSE
      )
    }
    rn <- rownames(Sigma)
    cn <- colnames(Sigma)
    if (!is.null(rn) || !is.null(cn)) {
      if (is.null(rn) || is.null(cn) || !setequal(rn, vars) || !setequal(cn, vars)) {
        stop(
          "Row and column names of `Sigma` must both match the variables: ",
          paste(vars, collapse = ", "), ".",
          call. = FALSE
        )
      }
      Sigma <- Sigma[vars, vars]
    }
    Sigma <- unname(Sigma)
    # chol() reads only the upper triangle
    if (!isSymmetric(Sigma)) {
      stop("`Sigma` must be symmetric.", call. = FALSE)
    }
    rvnorm_check_positive_definite(Sigma, "Sigma")
    return(Sigma)
  }
  if (is.null(sd)) {
    stop("`sd` must be supplied unless `Sigma` is supplied.", call. = FALSE)
  }
  if (is.matrix(sd)) {
    stop("Use `Sigma` to supply a covariance matrix.", call. = FALSE)
  }
  rvnorm_check_positive_numeric(sd, "sd")
  if (length(sd) == 1L) sd <- rep(as.numeric(sd), n_vars)
  if (length(sd) != n_vars) {
    stop("Length of `sd` must match the number of variables.", call. = FALSE)
  }
  diag(as.numeric(by_name(sd, "sd"))^2, nrow = n_vars)
}

rvnorm_widen_window <- function(window, vars, poly_vars, margin) {
  # a window from rvnorm_normalize_window(), widened by margin in the
  # variables poly involves. Projection does not move the others, so they
  # keep their bounds, which turns a scalar window into per-variable bounds;
  # without such variables it stays scalar, for the templates
  if (window$scalar) {
    if (all(vars %in% poly_vars)) return(window$value + margin)
    bounds <- rep(list(c(-window$value, window$value)), length(vars))
    window <- list(value = stats::setNames(bounds, vars), scalar = FALSE)
  }
  out <- lapply(names(window$value), function(var) {
    bounds <- window$value[[var]]
    if (var %in% poly_vars) bounds + c(-margin, margin) else bounds
  })
  names(out) <- names(window$value)
  out
}

rvnorm_in_window <- function(x, window, slack = 0) {
  # whether each row of x lies in the box rvnorm_normalize_window()
  # describes, to within slack (each row's precision) plus the rounding of
  # the bounds. A projection stops within its precision of the variety, so
  # a point of it on the box's edge, as the root 1 of x^2 - x is for w = 1,
  # can land just outside (1 + 2e-16), and an exact test dropped it
  if (is.null(window)) return(rep(TRUE, nrow(x)))
  if (window$scalar) {
    margin <- slack + 8 * .Machine$double.eps * window$value
    return(rowSums(abs(x) > window$value + margin) == 0)
  }
  inside <- rep(TRUE, nrow(x))
  for (var in names(window$value)) {
    bounds <- window$value[[var]]
    margin <- slack + 8 * .Machine$double.eps * max(abs(bounds))
    inside <- inside & x[, var] >= bounds[1] - margin &
      x[, var] <= bounds[2] + margin
  }
  inside
}

rvnorm_with_seed <- function(seed, code) {
  # evaluate code with R's RNG seeded by seed, then restore the caller's RNG
  # state; without a seed, evaluate it as is
  if (is.null(seed)) return(code)
  env <- globalenv()
  had_state <- exists(".Random.seed", envir = env, inherits = FALSE)
  if (had_state) old_state <- get(".Random.seed", envir = env, inherits = FALSE)
  on.exit(
    if (had_state) {
      assign(".Random.seed", old_state, envir = env)
    } else if (exists(".Random.seed", envir = env, inherits = FALSE)) {
      rm(".Random.seed", envir = env)
    },
    add = TRUE
  )
  set.seed(as.integer(seed %% .Machine$integer.max))
  code
}
