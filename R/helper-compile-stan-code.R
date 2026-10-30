canonicalize_mpoly <- function(poly) {
  # normalize term ordering so structurally equivalent polynomials hash equally
  terms <- unclass(poly)

  # order variables within each term too: mp("y x") and mp("x y") store the
  # same monomial with differently ordered entries
  terms <- lapply(terms, function(term) {
    is_var <- names(term) != "coef"
    c(term[is_var][order(names(term)[is_var], method = "radix")], term[!is_var])
  })

  degrees <- sapply(
    terms,
    function(term) sum(term[names(term) != "coef"], na.rm = TRUE)
  )
  lex_keys <- sapply(
    terms,
    function(term) {
      paste(sort(names(term)[names(term) != "coef"]), collapse = "")
    }
  )
  # full monomial keys break ties between same-degree terms in the same
  # variables (e.g. x^2 y and x y^2), which otherwise kept input order
  mono_keys <- vapply(terms, monomial_key, character(1))

  # sort: highest total degree first, then lexicographic variable order
  sorted_terms <- terms[order(-degrees, lex_keys, mono_keys, method = "radix")]

  structure(sorted_terms, class = "mpoly")
}

monomial_key <- function(term) {
  # coefficient-free label of a term, e.g. "x2 y1"
  is_var <- names(term) != "coef"
  paste0(names(term)[is_var], term[is_var], collapse = " ")
}

canonicalize_mpolylist <- function(poly) {
  # canonicalize each polynomial in a list
  stopifnot(inherits(poly, "mpolyList"))

  canonicalized_list <- lapply(poly, canonicalize_mpoly)
  structure(canonicalized_list, class = "mpolyList")
}

sort_mpolylist_lexicographically <- function(poly) {
  # enforce deterministic polynomial ordering before hashing/codegen; the
  # keys are the canonical monomials alone, so changing a coefficient cannot
  # reorder the equations and change the model key, code, or coefficient data
  stopifnot(inherits(poly, "mpolyList"))

  structure_keys <- vapply(
    poly,
    function(p) {
      paste(vapply(unclass(p), monomial_key, character(1)), collapse = " + ")
    },
    character(1)
  )
  sorted_indices <- order(structure_keys, method = "radix")

  sorted_list <- poly[sorted_indices]
  structure(sorted_list, class = "mpolyList")
}

helper_for_derivative_for_mpoly_stan_code <- function(var, poly, varorder) {
  # return derivative expression in Stan syntax for one variable
  lifted_poly <- stan_coef_lift(poly, varorder = varorder)
  derivative <- deriv(lifted_poly, var)
  mpoly_to_stan(derivative)
}

helper_for_coef_lift_for_mpolylist <- function(p, i, varorder = sort(mpoly::vars(p))) {
  # replace numeric coefficients with symbolic names (e.g. bx2_1) for Stan
  stan_coef_lift(p, varorder = varorder, suffix = i)
}

helper_for_derivative_for_mpolylist_stan_code <- function(
    var, poly, i, varorder = sort(mpoly::vars(poly))
  ) {
  # derivative helper for mpolyList code generation
  lifted_poly <- helper_for_coef_lift_for_mpolylist(poly, i, varorder = varorder)
  derivative <- deriv(lifted_poly, var)
  mpoly_to_stan(derivative)
}

coef_lift_mpolylist_for_generating_names <- function(poly) {
  # lift all polynomials so model naming depends on structure, not coefficients
  lifted_mpolylist <- lapply(poly, coef_lift)
  class(lifted_mpolylist) <- "mpolyList"
  lifted_mpolylist
}

generate_model_name <- function(poly, windowed = FALSE, homo = TRUE) {
  # build deterministic cache key from canonicalized polynomial structure
  if (is.mpoly(poly)) {
    g <- canonicalize_mpoly(poly)
    g <- coef_lift(g)
    g <- digest::digest(g, algo = "md5")
  } else if (is.mpolyList(poly)) {
    g <- canonicalize_mpolylist(poly)
    g <- sort_mpolylist_lexicographically(g)
    g <- coef_lift_mpolylist_for_generating_names(g)
    g <- digest::digest(g, algo = "md5")
  }

  sprintf(
    "%s_%s%s.stan",
    g,
    if (homo) "vn" else "hvn",
    if (windowed) "_w" else ""
  )
}
