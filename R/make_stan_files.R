# internal helper to generate Stan files for precompiled univariate models

make_stan_files <- function(
  num_of_vars,
  totaldeg,
  homo = TRUE,
  windowed = TRUE,
  basis = c("x", "y", "z")
) {
  # enumerate parameter variable names used by this template
  vars <- basis[seq_len(num_of_vars)]

  # build the Stan data block containing polynomial coefficients
  var_for_data_block <- mpoly::basis_monomials(
    basis[seq_len(num_of_vars)],
    totaldeg
  )
  var_for_data_block <- lapply(var_for_data_block, reorder, varorder = basis)
  var_for_data_block <- lapply(var_for_data_block, coef)
  var_for_data_block <- unlist(var_for_data_block)
  var_for_data_block <- get_listed_coeficients(var_for_data_block)
  var_for_data_block <- names(var_for_data_block)

  data_block <- paste(
    sapply(var_for_data_block, function(x) paste("  real", x)),
    collapse = "; "
  )
  data_block <- paste0(data_block, ";")

  # build the Stan parameter block (optionally box-constrained by w)
  if (windowed) {
    data_block <- paste0(data_block, "\n  real w;")
  }

  data_block <- paste0("data {\n  real si;\n", data_block, "\n}\n")

  if (windowed) {
    params_block <- paste(sapply(vars, function(var) {
      paste0("  real<lower=-", "w", ", upper=", "w", "> ", var, ";")
    }), collapse = "\n")
  } else {
    params_block <- paste(sapply(vars, function(var) {
      paste0("  real ", var, ";")
    }), collapse = "\n")
  }

  params_block <- paste0("parameters {\n", params_block, "\n}\n")

  # build the symbolic polynomial g and its derivatives for transformed params
  g_coef <- mpoly::basis_monomials(basis[seq_len(num_of_vars)], totaldeg)
  g_coef <- lapply(g_coef, reorder, varorder = basis)
  g_coef <- lapply(g_coef, coef)
  g_coef <- unlist(g_coef)
  g_coef <- get_listed_coeficients(g_coef)
  g_coef <- names(g_coef)

  g_terms <- mpoly::basis_monomials(basis[seq_len(num_of_vars)], totaldeg)
  g_terms <- lapply(g_terms, reorder, varorder = basis)
  g_terms <- lapply(g_terms, mpoly_to_stan)
  g_terms <- unlist(g_terms)
  g_terms <- c(g_terms)

  g <- paste0(g_coef, "*", g_terms, collapse = "+")
  g <- gsub("1\\*|\\*1", "", g)

  derivatives <- lapply(
    vars,
    get_derivative,
    num_of_vars = num_of_vars,
    deg = totaldeg,
    basis = vars
  )

  derivative_names <- sapply(seq_along(vars), function(i) {
    paste0("dg", vars[i])
  })

  dg <- sapply(seq_along(vars), function(i) {
    paste(
      derivative_names[i],
      paste(derivatives[[i]], collapse = ""),
      sep = " = "
    )
  })

  dg <- paste0("  real ", dg)
  dg <- paste(dg, collapse = ";\n")

  if (homo) {
    ndg <- paste0("  real ndg = sqrt(",
                  paste0(derivative_names, "^2", collapse = " + "),
                  ");")
  } else {
    ndg <- "  real ndg = 1;"
  }

  # emit transformed-parameters and model blocks
  trans_block <- paste0(
    "transformed parameters {\n  real g = ", g, ";\n", dg, ";\n", ndg, "\n}\n"
  )

  model_block <- paste0(
    "model {\n  target += normal_lpdf(0.00 | g/ndg, si); \n}"
  )

  # assemble complete Stan program text
  stan_code <- paste0(
    data_block,
    params_block,
    trans_block,
    model_block,
    sep = ""
  )

  filename <- sprintf(
    "%i_%i_%s%s.stan",
    length(vars),
    totaldeg,
    if (homo) "vn" else "hvn",
    if (windowed) "_w" else ""
  )

  # write template source into src/Stan for package compilation
  path <- file.path("src", "stan", filename)
  writeLines(stan_code, con = path)
}
