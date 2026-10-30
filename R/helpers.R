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

generic_jacobian_check <- function(poly, varorder) {
  # rank of a polynomial system's jacobian at a few fixed points with
  # irrational-looking coordinates, which miss the proper algebraic subsets
  # (singular loci) where an ordinary system's J loses rank, so rank
  # deficiency at all of them means dependent equations. Also whether g
  # leaves J's column space at all of them, where (for m <= n) pdvnorm()'s
  # density is 0
  m <- length(poly)
  n <- length(varorder)
  J_fn <- jacobian_function(poly, varorder)
  g_fn <- poly_function(poly, varorder)
  pts <- matrix(2 * ((seq_len(3L * n) * 0.6180339887498949) %% 1) - 1, nrow = 3L)
  checks <- lapply(seq_len(nrow(pts)), function(i) {
    J <- J_fn(pts[i, ])
    rank <- attr(pdvnorm_pseudoinverse(J), "rank")
    outside <- rank < m &&
      pdvnorm_outside_column_space(J, matrix(g_fn(pts[i, ]), nrow = 1))
    c(rank = rank, outside = outside)
  })
  checks <- do.call(rbind, checks)
  list(
    full_rank = max(checks[, "rank"]) == min(m, n),
    density_vanishes = m <= n && all(checks[, "outside"] == 1)
  )
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
