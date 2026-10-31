make_coefficients_data <- function(
    poly, num_of_vars, deg, basis = c("x", "y", "z")
  ) {
  # enumerate all monomials up to degree; zero-fill any not in poly
  vars <- basis[seq_len(num_of_vars)]
  monos <- mpoly::basis_monomials(vars, deg)
  stan_names <- stan_coef_names(monos, varorder = basis)

  required_coefs <- as.list(stats::setNames(rep(0L, length(stan_names)), stan_names))
  available_coef <- stan_coefficients(poly, varorder = basis)
  unknown_coefs <- setdiff(names(available_coef), stan_names)
  if (length(unknown_coefs) > 0L) {
    stop(
      "Coefficient name(s) do not match the precompiled Stan template: ",
      paste(unknown_coefs, collapse = ", "),
      call. = FALSE
    )
  }
  required_coefs[names(available_coef)] <- available_coef
  required_coefs
}

get_coefficients_data <- function(poly) {
  # collect named coefficients for either a single polynomial
  # or polynomial list
  if (is.mpoly(poly)) {
    data <- stan_coefficients(poly, varorder = sort(mpoly::vars(poly)))
  } else if (is.mpolyList(poly)) {
    poly <- canonicalize_mpolylist(poly)
    poly <- sort_mpolylist_lexicographically(poly)
    coefs <- list()
    for (i in seq_along(poly)) {
      coefs[[i]] <- stan_coefficients(
        poly[[i]],
        varorder = sort(mpoly::vars(poly[[i]])),
        suffix = i
      )
    }
    coefs <- unlist(coefs)
    data <- as.list(coefs)
  }
  data
}

write_stan_data <- function(data) {
  # cmdstanr writes data with 15 significant digits, which rounds
  # coefficients that need more; Stan reads 17 back exactly
  path <- tempfile(pattern = "vnorm-data-", fileext = ".json")
  jsonlite::write_json(data, path, auto_unbox = TRUE, digits = I(17), pretty = TRUE)
  path
}

check_and_replace_vars <- function(p, vars = mpoly::vars(p)) {
  # map arbitrary variable names to x/y/z for precompiled Stan templates;
  # `vars` may add variables that p does not involve, which take the
  # remaining template names
  current_vars <- vars
  num_vars <- length(current_vars)
  target_vars <- list(
    c("x"),           # for 1 indeterminate
    c("x", "y"),      # for 2 indeterminates
    c("x", "y", "z")  # for 3 indeterminates
  )

  if (num_vars > 3) {
    stop("The polynomial has more than 3 indeterminates.", call. = FALSE)
  }

  expected_vars <- target_vars[[num_vars]]
  if (setequal(current_vars, expected_vars)) {
    return(list(polynomial = p, mapping = list()))
  }

  # two-pass swap via temporary names to avoid collisions (e.g. x->y, y->x);
  # swap() only renames variables p involves
  temp_vars <- paste0("tmp", seq_along(current_vars))
  in_p <- which(current_vars %in% mpoly::vars(p))

  for (i in in_p) {
    p <- swap(p, current_vars[i], temp_vars[i])
  }

  for (i in in_p) {
    p <- swap(p, temp_vars[i], expected_vars[i])
  }

  list(
    polynomial = p,
    mapping = stats::setNames(as.list(current_vars), expected_vars)
  )
}

rvnorm_resolve_vars <- function(vars, poly) {
  # the variables to sample, sorted: those of poly unless `vars` is given,
  # in which case it must include them
  poly_vars <- mpoly::vars(poly)
  if (is.null(vars)) {
    return(sort(poly_vars))
  }
  if (
    !is.character(vars) || length(vars) == 0L || anyNA(vars) ||
      any(!nzchar(vars)) || anyDuplicated(vars)
  ) {
    stop("`vars` must be a character vector of unique variable names.", call. = FALSE)
  }
  missing_vars <- setdiff(poly_vars, vars)
  if (length(missing_vars) > 0L) {
    stop(
      "`vars` must include every variable of `poly`; it lacks ",
      paste(missing_vars, collapse = ", "), ".",
      call. = FALSE
    )
  }
  sort(vars)
}

rvnorm_normalize_window <- function(w, vars) {
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
    return(list(value = as.numeric(w), scalar = TRUE))
  }

  if (is.numeric(w) && length(w) == 2L) {
    bounds <- check_bounds(w, "`w`")
    out <- replicate(length(vars), bounds, simplify = FALSE)
    names(out) <- vars
    return(rvnorm_collapse_window(out, vars))
  }

  if (!is.list(w) || is.null(names(w)) || any(!nzchar(names(w)))) {
    stop("`w` must be a positive scalar, a length-2 interval, or a named list.", call. = FALSE)
  }

  unknown_vars <- setdiff(names(w), vars)
  if (length(unknown_vars) > 0L) {
    stop(
      "`w` contains bounds for unknown variable(s): ",
      paste(unknown_vars, collapse = ", "),
      call. = FALSE
    )
  }

  out <- lapply(names(w), function(var) check_bounds(w[[var]], paste0("`w$", var, "`")))
  names(out) <- names(w)
  rvnorm_collapse_window(out, vars)
}

rvnorm_collapse_window <- function(w, vars) {
  # the same box (-a, a) on every variable is the scalar window a, which the
  # pre-compiled templates take
  b <- w[[1]]
  same_box <- setequal(names(w), vars) &&
    all(vapply(w, identical, logical(1), b)) && b[1] == -b[2]
  if (same_box) list(value = b[2], scalar = TRUE) else list(value = w, scalar = FALSE)
}

rename_output_df <- function(df, replacement_list) {
  # restore original variable names in output after x/y/z remapping
  names(df) <- sapply(names(df), function(col) {
    if (col %in% names(replacement_list)) {
      replacement_list[[col]]
    } else {
      col
    }
  })
  df
}
