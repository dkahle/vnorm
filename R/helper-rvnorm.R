make_coefficients_data <- function(
    poly, num_of_vars, deg, basis = c("x", "y", "z")
  ) {
  # enumerate all monomials up to degree; zero-fill any not in poly
  monos <- mpoly::basis_monomials(basis[seq_len(num_of_vars)], deg)

  # single pass: reorder, extract coef name, convert to Stan name
  stan_names <- vapply(monos, function(m) {
    nm <- names(coef(reorder(m, varorder = basis)))
    nm <- gsub("\\s+", "", nm)
    if (nm == "1") return("b1")
    paste0("b", gsub("\\^", "", nm))
  }, character(1))

  required_coefs <- as.list(stats::setNames(rep(0L, length(stan_names)), stan_names))
  available_coef <- get_listed_coeficients(coef(poly))
  required_coefs[names(available_coef)] <- available_coef
  required_coefs
}

get_coefficients_data <- function(poly) {
  # collect named coefficients for either a single polynomial
  # or polynomial list
  if (is.mpoly(poly)) {
    data <- get_listed_coeficients(coef(poly))
  } else if (is.mpolyList(poly)) {
    poly <- canonicalize_mpolylist(poly)
    poly <- sort_mpolylist_lexicographically(poly)
    # suffix coefficient names with polynomial index (e.g. bx2_1, by_2)
    convert_names <- function(term, i) {
      term <- gsub("\\s+", "", term)
      if (term == "1") return(paste0("b1_", i))
      term <- gsub("\\^", "", term)
      paste0("b", term, "_", i)
    }
    coefs <- list()
    for (i in seq_along(poly)) {
      coefs[[i]] <- coef(poly[[i]])
      names(coefs[[i]]) <- sapply(names(coefs[[i]]), convert_names, i = i)
    }
    coefs <- unlist(coefs)
    data <- as.list(coefs)
  }
  data
}

check_and_replace_vars <- function(p) {
  # map arbitrary variable names to x/y/z for precompiled Stan templates
  current_vars <- mpoly::vars(p)
  num_vars <- length(current_vars)
  target_vars <- list(
    c("x"),           # for 1 indeterminate
    c("x", "y"),      # for 2 indeterminates
    c("x", "y", "z")  # for 3 indeterminates
  )

  if (num_vars > 3) {
    stop("The polynomial has more than 3 indeterminates.")
  }

  expected_vars <- target_vars[[num_vars]]
  if (setequal(current_vars, expected_vars)) {
    return(list(polynomial = p, mapping = list()))
  }

  # two-pass swap via temporary names to avoid collisions (e.g. x->y, y->x)
  var_mapping <- list()
  temp_vars <- paste0("tmp", seq_along(current_vars))

  for (i in seq_along(current_vars)) {
    p <- swap(p, current_vars[i], temp_vars[i])
    var_mapping[[temp_vars[i]]] <- current_vars[i]
  }

  for (i in seq_along(expected_vars)) {
    p <- swap(p, temp_vars[i], expected_vars[i])
    var_mapping[[expected_vars[i]]] <- var_mapping[[temp_vars[i]]]
    var_mapping[[temp_vars[i]]] <- NULL
  }

  list(polynomial = p, mapping = var_mapping)
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
