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

gradient_mpolylist <- function(poly, varorder) {
  # partial derivatives in varorder order, always as an mpolyList: deriv()
  # returns a bare mpoly for a single variable, and as.function() rejects a
  # bare constant mpoly when given a varorder
  dp <- deriv(poly, var = varorder)
  if (is.mpoly(dp)) dp <- structure(list(dp), class = "mpolyList")
  dp
}

gradient_function <- function(poly, varorder) {
  # gradient of poly as a function of a point, entries in varorder order
  as.function(gradient_mpolylist(poly, varorder), varorder = varorder, silent = TRUE)
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
