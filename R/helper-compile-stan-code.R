canonicalize_mpoly <- function(poly) {
  # normalize term ordering so structurally equivalent polynomials hash equally
  terms <- unclass(poly)

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

  # sort: highest total degree first, then lexicographic variable order
  sorted_terms <- terms[order(-degrees, lex_keys, method = "radix")]

  structure(sorted_terms, class = "mpoly")
}

canonicalize_mpolylist <- function(poly) {
  # canonicalize each polynomial in a list
  stopifnot(inherits(poly, "mpolyList"))

  canonicalized_list <- lapply(poly, canonicalize_mpoly)
  structure(canonicalized_list, class = "mpolyList")
}

sort_mpolylist_lexicographically <- function(poly) {
  # enforce deterministic polynomial ordering before hashing/codegen
  stopifnot(inherits(poly, "mpolyList"))

  poly_strings <- vapply(
    poly,
    function(p) paste(as.character(p), collapse = ""),
    character(1)
  )
  sorted_indices <- order(poly_strings, method = "radix")

  sorted_list <- poly[sorted_indices]
  structure(sorted_list, class = "mpolyList")
}

helper_for_derivative_for_mpoly_stan_code <- function(var, poly) {
  # return derivative expression in Stan syntax for one variable
  lifted_poly <- mpoly::coef_lift(poly)
  derivative <- deriv(lifted_poly, var)
  mpoly_to_stan(derivative)
}

helper_for_coef_lift_for_mpolylist <- function(p, i) {
  # replace numeric coefficients with symbolic names (e.g. bx2_1) for Stan
  monos <- monomials(p, unit = TRUE)
  printed_monos <- print(monos, silent = TRUE)
  printed_monos <- gsub("\\^", "", printed_monos)
  printed_monos <- gsub(" ", "", printed_monos)
  coefs_to_add <- paste0("b", printed_monos, "_", i)
  for (j in seq_along(p)) {
    p[[j]]["coef"] <- 1
    p[[j]] <- structure(c(1, p[[j]]), names = c(coefs_to_add[j], names(p[[j]])))
  }
  p
}

helper_for_derivative_for_mpolylist_stan_code <- function(var, poly, i) {
  # derivative helper for mpolyList code generation
  lifted_poly <- helper_for_coef_lift_for_mpolylist(poly, i)
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
