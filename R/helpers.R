stan_coef_name_from_term <- function(term, varorder) {
  if (is.mpoly(term)) {
    if (length(term) != 1L) {
      stop("`term` must be a single monomial.", call. = FALSE)
    }
    term <- unclass(term)[[1]]
  }

  powers <- term[names(term) != "coef"]
  if (length(powers) == 0L) return("b1")

  missing_vars <- setdiff(names(powers), varorder)
  if (length(missing_vars) > 0L) {
    stop(
      "Coefficient term contains variable(s) outside `varorder`: ",
      paste(missing_vars, collapse = ", "),
      call. = FALSE
    )
  }

  ordered_vars <- varorder[varorder %in% names(powers)]
  pieces <- vapply(ordered_vars, function(var) {
    power <- powers[[var]]
    if (identical(unname(power), 1)) var else paste0(var, power)
  }, character(1))

  paste0("b", paste0(pieces, collapse = ""))
}

stan_coef_names <- function(poly, varorder = sort(mpoly::vars(poly)), suffix = NULL) {
  if (is.mpoly(poly)) {
    terms <- unclass(poly)
  } else if (is.mpolyList(poly)) {
    terms <- unlist(lapply(poly, unclass), recursive = FALSE)
  } else {
    stop("`poly` must be an mpoly or mpolyList object.", call. = FALSE)
  }

  out <- vapply(terms, stan_coef_name_from_term, character(1), varorder = varorder)
  if (!is.null(suffix)) out <- paste0(out, "_", suffix)
  out
}

stan_coefficients <- function(poly, varorder = sort(mpoly::vars(poly)), suffix = NULL) {
  if (!is.mpoly(poly)) {
    stop("`poly` must be an mpoly object.", call. = FALSE)
  }

  coefs <- vapply(unclass(poly), function(term) unname(term[["coef"]]), numeric(1))
  names(coefs) <- stan_coef_names(poly, varorder = varorder, suffix = suffix)
  as.list(coefs)
}

stan_coef_lift <- function(poly, varorder = sort(mpoly::vars(poly)), suffix = NULL) {
  if (!is.mpoly(poly)) {
    stop("`poly` must be an mpoly object.", call. = FALSE)
  }

  coef_names <- stan_coef_names(poly, varorder = varorder, suffix = suffix)
  lifted <- unclass(poly)
  for (i in seq_along(lifted)) {
    lifted[[i]][["coef"]] <- 1
    lifted[[i]] <- structure(
      c(1, lifted[[i]]),
      names = c(coef_names[[i]], names(lifted[[i]]))
    )
  }
  structure(lifted, class = "mpoly")
}

get_listed_coefficients <- function(coefs) {
  # Backward-compatible printed-name conversion for legacy internals/tests.
  convert_names <- function(term) {
    term <- gsub("\\s+", "", term)
    if (term == "1") return("b1")
    term <- gsub("\\^", "", term)
    paste0("b", term)
  }
  names(coefs) <- sapply(names(coefs), convert_names)
  as.list(coefs)
}

poly_function <- function(x, varorder) {
  # evaluate an mpoly (or mpolyList) at points, coordinates in varorder order:
  # a matrix holds one point per row, and a vector is one point (or, with one
  # variable, one point per entry). mpoly's as.function() parses the printed
  # polynomial, so its coefficients follow getOption("digits"), and R's
  # parser misreads some 17-digit decimals besides; build the function from
  # the terms instead, with the coefficients themselves inlined
  is_list <- is.mpolyList(x)
  polys <- if (is_list) unclass(x) else list(x)
  unknown_vars <- setdiff(unlist(lapply(polys, mpoly::vars)), varorder)
  if (length(unknown_vars) > 0L) {
    stop(
      "`varorder` is missing variable(s): ",
      paste(unknown_vars, collapse = ", "),
      call. = FALSE
    )
  }
  n <- length(varorder)
  cols <- lapply(paste0(".x", seq_len(n)), as.name)

  term_call <- function(term) {
    powers <- term[names(term) != "coef"]
    out <- unname(term[["coef"]])
    for (var in names(powers)) {
      col <- cols[[match(var, varorder)]]
      power <- unname(powers[[var]])
      out <- call("*", out, if (power == 1) col else call("^", col, power))
    }
    out
  }
  poly_call <- function(p) {
    # rep_len() gives constant polynomials one value per point
    terms <- lapply(unclass(p), term_call)
    call("rep_len", Reduce(function(a, b) call("+", a, b), terms), quote(.k))
  }
  values <- if (is_list) {
    # one column per polynomial; a single point given as a vector gives a
    # plain vector
    substitute(
      {
        .out <- VALUES
        if (.k == 1L && !is.matrix(.)) .out[1L, ] else .out
      },
      list(VALUES = as.call(c(as.name("cbind"), lapply(polys, poly_call))))
    )
  } else {
    poly_call(polys[[1L]])
  }

  from_matrix <- lapply(seq_len(n), function(j) {
    substitute(COL <- .[, J], list(COL = cols[[j]], J = j))
  })
  from_vector <- if (n == 1L) {
    list(quote(.x1 <- .))
  } else {
    lapply(seq_len(n), function(j) {
      substitute(COL <- .[J], list(COL = cols[[j]], J = j))
    })
  }

  f <- function(.) NULL
  body(f) <- substitute(
    {
      if (is.matrix(.)) {
        if (ncol(.) != N) {
          stop("A matrix of points must have one column per variable.", call. = FALSE)
        }
        FROM_MATRIX
        .k <- nrow(.)
      } else {
        if (N > 1L && length(.) != N) {
          stop("A point must have one entry per variable.", call. = FALSE)
        }
        FROM_VECTOR
        .k <- length(.x1)
      }
      VALUES
    },
    list(
      N = n,
      FROM_MATRIX = as.call(c(as.name("{"), from_matrix)),
      FROM_VECTOR = as.call(c(as.name("{"), from_vector)),
      VALUES = values
    )
  )
  environment(f) <- baseenv()
  f
}

gradient_mpolylist <- function(poly, varorder) {
  # partial derivatives in varorder order, always as an mpolyList: deriv()
  # returns a bare mpoly for a single variable
  dp <- deriv(poly, var = varorder)
  if (is.mpoly(dp)) dp <- structure(list(dp), class = "mpolyList")
  dp
}

gradient_function <- function(poly, varorder) {
  # gradient of poly as a function of a point, entries in varorder order
  poly_function(gradient_mpolylist(poly, varorder), varorder = varorder)
}

hessian_function <- function(poly, varorder) {
  # hessian of poly as a function of a point, rows/columns in varorder order
  n <- length(varorder)
  row_fns <- lapply(
    gradient_mpolylist(poly, varorder),
    gradient_function,
    varorder = varorder
  )
  function(x) matrix(vapply(row_fns, function(f) f(x), numeric(n)), n, n)
}

jacobian_function <- function(poly, varorder) {
  # jacobian of a polynomial system as a function of a point: one row per
  # polynomial, columns in varorder order
  row_fns <- lapply(poly, gradient_function, varorder = varorder)
  n <- length(varorder)
  function(x) {
    matrix(
      vapply(row_fns, function(f) f(x), numeric(n)),
      nrow = length(row_fns), byrow = TRUE
    )
  }
}

dependent_equations <- function(poly, varorder) {
  # whether a polynomial system's jacobian has rank below min(m, n) at every
  # point. Clear full rank at one of a few fixed points settles that it does
  # not, once each row is scaled to unit size (x^100 - 1 has gradient
  # 100 x^99, about 1e-20 at x = 0.6). Rank deficiency there proves nothing,
  # since badly scaled or nearly dependent equations show it too, so the rank
  # is then computed exactly, in integers. When that cannot be done exactly
  # or quickly, the check is inconclusive and reports independence
  m <- length(poly)
  n <- length(varorder)
  k <- min(m, n)
  J_fn <- jacobian_function(poly, varorder)
  pts <- matrix(2 * ((seq_len(3L * n) * 0.6180339887498949) %% 1) - 1, nrow = 3L)
  for (i in seq_len(nrow(pts))) {
    J <- J_fn(pts[i, ])
    if (!all(is.finite(J))) next
    J <- J / pmax(apply(abs(J), 1L, max), .Machine$double.xmin)
    d <- svd(J, nu = 0L, nv = 0L)$d
    if (d[k] > sqrt(.Machine$double.eps) * d[1L]) return(FALSE)
  }
  rows <- integer_jacobian_rows(poly, varorder)
  !is.null(rows) && isTRUE(polynomial_matrix_rank(rows, max_work = 1e5) < k)
}

integer_jacobian_rows <- function(poly, varorder) {
  # the jacobian's rows (lists of mpoly) with integer coefficients, for
  # polynomial_matrix_rank(). Each equation is divided by the greatest common
  # divisor of its coefficients, an odd integer times a power of two, which
  # changes no rank and gives 0.1 x + 0.2 y (0.2 is exactly twice 0.1) the
  # coefficients 1 and 2. NULL once a coefficient reaches 2^53, past which
  # doubles do not hold every integer, as unrelated decimals' 53-bit
  # expansions do
  rows <- lapply(poly, function(p) {
    terms <- Filter(function(term) length(term) > 1L, unclass(p))
    if (length(terms) == 0L) return(rep(list(mp("0")), length(varorder)))
    parts <- binary_parts(vapply(terms, function(term) term[["coef"]], numeric(1)))
    gcd <- function(a, b) {
      # floor(a / b) is exact for integers below 2^53, where %% warns of a
      # loss of accuracy it does not have
      while (b > 0) {
        r <- a - floor(a / b) * b
        a <- b
        b <- r
      }
      a
    }
    coefs <- parts$m / Reduce(gcd, abs(parts$m)) * 2^(parts$e - min(parts$e))
    lapply(varorder, function(v) {
      # differentiate term by term: mpoly's deriv() rounds coef * power
      d_terms <- list()
      for (i in seq_along(terms)) {
        powers <- terms[[i]][names(terms[[i]]) != "coef"]
        if (!v %in% names(powers)) next
        power <- powers[[v]]
        powers[[v]] <- power - 1
        d_terms[[length(d_terms) + 1L]] <- c(powers[powers > 0], coef = coefs[i] * power)
      }
      if (length(d_terms) == 0L) mp("0") else mpoly::mpoly(d_terms)
    })
  })
  coefs <- unlist(lapply(rows, function(row) lapply(row, function(p) {
    vapply(unclass(p), function(term) term[["coef"]], numeric(1))
  })))
  if (!all(is.finite(coefs)) || any(abs(coefs) >= 2^53)) return(NULL)
  rows
}

binary_parts <- function(x) {
  # x = m * 2^e exactly for finite nonzero x, with m an odd integer: dividing
  # by a power of two is exact, and log2() can round across one
  e <- floor(log2(abs(x)))
  e <- e - (abs(x) < 2^e) + (abs(x) >= 2^(e + 1))
  e <- pmax(e - 52, -1074)
  m <- x / 2^e
  repeat {
    even <- m / 2 == floor(m / 2)
    if (!any(even)) break
    m[even] <- m[even] / 2
    e[even] <- e[even] + 1
  }
  list(m = m, e = e)
}

polynomial_matrix_rank <- function(rows, max_work = Inf) {
  # rank of a matrix of integer polynomials (a list of rows of mpoly entries)
  # at generic points, i.e. over the rational functions, by fraction-free
  # elimination with exact zero tests on the coefficients. Doubles hold
  # integers exactly below 2^53, so this gives NA once an update could reach
  # that (each coefficient of p a - f b is at most |p||a| + |f||b|, with |.|
  # the sum of absolute coefficients), and once the term products it would
  # multiply out pass max_work (about 30 microseconds each)
  is_zero <- function(p) all(vapply(unclass(p), function(term) term[["coef"]] == 0, logical(1)))
  l1 <- function(p) sum(abs(vapply(unclass(p), function(term) term[["coef"]], numeric(1))))
  rank <- 0L
  work <- 0
  repeat {
    if (length(rows) == 0L || length(rows[[1L]]) == 0L) return(rank)
    size <- do.call(rbind, lapply(rows, function(row) {
      vapply(row, function(p) if (is_zero(p)) Inf else length(p), numeric(1))
    }))
    if (all(is.infinite(size))) return(rank)
    # the sparsest pivot keeps the entries small
    pivot <- arrayInd(which.min(size), dim(size))
    pivot_row <- rows[[pivot[1L]]]
    j <- pivot[2L]
    p <- pivot_row[[j]]
    rows <- rows[-pivot[1L]]
    to_clear <- !vapply(rows, function(row) is_zero(row[[j]]), logical(1))
    work <- work + sum(vapply(rows[to_clear], function(row) {
      length(p) * sum(lengths(row[-j])) + length(row[[j]]) * sum(lengths(pivot_row[-j]))
    }, numeric(1)))
    if (work > max_work) return(NA_integer_)
    # 2^52 leaves room for rounding in the bound's own sums
    for (row in rows[to_clear]) {
      for (col in seq_along(pivot_row)[-j]) {
        bound <- l1(p) * l1(row[[col]]) + l1(row[[j]]) * l1(pivot_row[[col]])
        if (bound >= 2^52) return(NA_integer_)
      }
    }
    rows <- lapply(seq_along(rows), function(i) {
      row <- rows[[i]]
      if (!to_clear[i]) return(row[-j])
      f <- row[[j]]
      Map(function(a, b) p * a - f * b, row[-j], pivot_row[-j])
    })
    rank <- rank + 1L
  }
}

mpoly_to_stan <- function(mpoly) {
  # print polynomial in Stan-compatible syntax and parseable legend labels
  result <- print(
    mpoly,
    stars = TRUE,
    silent = TRUE,
    plus_pad = 0L,
    times_pad = 0L
  )
  result <- gsub("[*]{2}", "^", result)
  result
}

mpolyList_to_stan <- function(mpolyList) {
  # print a polynomial list as a Stan vector expression payload
  result <- print(
    mpolyList,
    silent = TRUE,
    stars = TRUE,
    plus_pad = 0,
    times_pad = 0
  )
  result <- gsub("\\*\\*", "^", result)
  result <- paste(result, collapse = ", ")
  result
}

get_derivative <- function(
    var, num_of_vars, deg, basis = c("x", "y", "z"), suffix = NULL
  ) {
  # construct symbolic derivative in Stan syntax for generated template models
  num_coef <- mpoly::basis_monomials(basis[seq_len(num_of_vars)], deg)
  num_coef <- lapply(num_coef, reorder, varorder = basis)
  num_coef <- lapply(num_coef, deriv, var = var)
  num_coef <- lapply(num_coef, coef)
  num_coef <- unlist(num_coef)
  num_coef <- unname(num_coef)

  indeterminates <- mpoly::basis_monomials(basis[seq_len(num_of_vars)], deg)
  indeterminates <- lapply(indeterminates, reorder, varorder = basis)
  indeterminates <- lapply(indeterminates, deriv, var = var)
  indeterminates <- lapply(indeterminates, function(item) {
    item[[1]][names(item[[1]]) == "coef"] <- 1
    item
  })
  indeterminates <- lapply(indeterminates, mpoly_to_stan)
  indeterminates <- unlist(indeterminates)
  indeterminates <- c(indeterminates)

  sym_coef <- mpoly::basis_monomials(basis[seq_len(num_of_vars)], deg)
  sym_coef <- stan_coef_names(sym_coef, varorder = basis, suffix = suffix)

  # pair numeric coefficients, symbolic names, and monomial terms for Stan output
  df_for_der <- data.frame(
    num_coef = num_coef,
    indeterminates = indeterminates,
    sym_coef = sym_coef
  )

  df_for_der <- dplyr::filter(df_for_der, num_coef != 0)

  # strip trivial 1* factors from the derivative expression
  gsub(
    "1\\*|\\*1",
    "",
    paste0(
      df_for_der$num_coef,
      "*",
      df_for_der$sym_coef,
      "*",
      df_for_der$indeterminates,
      collapse = "+"
    )
  )
}
