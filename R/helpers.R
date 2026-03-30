get_listed_coeficients <- function(coefs) {
  # convert monomial names to Stan data names (e.g., x^2 y -> bx2y)
  convert_names <- function(term) {
    term <- gsub("\\s+", "", term)
    if (term == "1") return("b1")
    term <- gsub("\\^", "", term)
    paste0("b", term)
  }
  names(coefs) <- sapply(names(coefs), convert_names)
  as.list(coefs)
}

mpoly_to_stan <- function(mpoly) {
  # print polynomial in Stan-compatible syntax and parseable legend labels
  p <- get("print.mpoly", asNamespace("mpoly"))
  result <- p(mpoly, stars = TRUE, silent = TRUE, plus_pad = 0L, times_pad = 0L)
  result <- gsub("[*]{2}", "^", result)
  result
}

mpolyList_to_stan <- function(mpolyList) {
  # print a polynomial list as a Stan vector expression payload
  p <- get("print.mpolyList", asNamespace("mpoly"))
  result <- p(
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

get_derivative <- function(var, num_of_vars, deg, basis = c("x", "y", "z")) {
  # construct symbolic derivative in Stan syntax for generated template models
  d <- get("deriv.mpoly", asNamespace("mpoly"))

  num_coef <- mpoly::basis_monomials(basis[seq_len(num_of_vars)], deg)
  num_coef <- lapply(num_coef, reorder, varorder = basis)
  num_coef <- lapply(num_coef, d, var = var)
  num_coef <- lapply(num_coef, mpoly:::coef.mpoly)
  num_coef <- unlist(num_coef)
  num_coef <- unname(num_coef)

  indeterminates <- mpoly::basis_monomials(basis[seq_len(num_of_vars)], deg)
  indeterminates <- lapply(indeterminates, reorder, varorder = basis)
  indeterminates <- lapply(indeterminates, d, var = var)
  indeterminates <- lapply(indeterminates, function(item) {
    item[[1]][names(item[[1]]) == "coef"] <- 1
    item
  })
  indeterminates <- lapply(indeterminates, mpoly_to_stan)
  indeterminates <- unlist(indeterminates)
  indeterminates <- c(indeterminates)

  sym_coef <- mpoly::basis_monomials(basis[seq_len(num_of_vars)], deg)
  sym_coef <- lapply(sym_coef, reorder, varorder = basis)
  sym_coef <- lapply(sym_coef, mpoly:::coef.mpoly)
  sym_coef <- unlist(sym_coef)
  sym_coef <- get_listed_coeficients(sym_coef)
  sym_coef <- names(sym_coef)

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
