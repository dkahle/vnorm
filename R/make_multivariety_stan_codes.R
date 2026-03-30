# internal helper to generate Stan code for multivariate variety models

make_multivariety_stan_codes <- function(num_of_vars,
                                         totaldeg,
                                         num_of_poly,
                                         homo = TRUE,
                                         windowed = TRUE,
                                         basis = c("x", "y", "z")) {
  # track variable sets for each polynomial in the system
  vars <- list()
  for (i in seq_len(num_of_poly)) {
    vars[[i]] <- basis[seq_len(num_of_vars[i])]
  }

  # build coefficient data declarations for every polynomial
  var_for_data_block <- list()
  for (i in seq_len(num_of_poly)) {
    var_for_data_block[[i]] <- mpoly::basis_monomials(
      basis[seq_len(num_of_vars[i])],
      totaldeg[i]
    )
    var_for_data_block[[i]] <- lapply(
      var_for_data_block[[i]],
      reorder,
      varorder = basis
    )
    var_for_data_block[[i]] <- lapply(var_for_data_block[[i]], coef)
    var_for_data_block[[i]] <- unlist(var_for_data_block[[i]])
    var_for_data_block[[i]] <- get_listed_coeficients(var_for_data_block[[i]])
    var_for_data_block[[i]] <- names(var_for_data_block[[i]])
    var_for_data_block[[i]] <- paste0(var_for_data_block[[i]], "_", i)
  }
  var_for_data_block <- unlist(var_for_data_block)
  data_block <- paste(
    sapply(var_for_data_block, function(x) paste("  real", x)),
    collapse = "; "
  )
  data_block <- paste0(data_block, ";")

  if (windowed) {
    data_block <- paste0(data_block, "\n  real w;")
  }

  data_block <- paste0("data {\n  real si;\n", data_block, "\n}\n")

  # parameter declarations for the union of variables across equations
  vars_for_params <- unique(unlist(vars))
  if (windowed) {
    params_block <- paste(sapply(vars_for_params, function(var) {
      paste0("  real<lower=-", "w", ", upper=", "w", "> ", var, ";")
    }), collapse = "\n")
  } else {
    params_block <- paste(sapply(vars_for_params, function(var) {
      paste0("  real ", var, ";")
    }), collapse = "\n")
  }
  params_block <- paste0("\nparameters {\n", params_block, "\n}\n")

  # build vector g and Jacobian entries used by the model
  g <- list()
  g_coef <- list()
  g_terms <- list()
  derivatives_pre <- list()
  derivatives <- list()

  for (i in seq_len(num_of_poly)) {
    g_coef[[i]] <- mpoly::basis_monomials(
      basis[seq_len(num_of_vars[i])],
      totaldeg[i]
    )
    g_coef[[i]] <- lapply(g_coef[[i]], reorder, varorder = basis)
    g_coef[[i]] <- lapply(g_coef[[i]], coef)
    g_coef[[i]] <- unlist(g_coef[[i]])
    g_coef[[i]] <- get_listed_coeficients(g_coef[[i]])
    g_coef[[i]] <- names(g_coef[[i]])

    g_terms[[i]] <- mpoly::basis_monomials(
      basis[seq_len(num_of_vars[i])],
      totaldeg[i]
    )
    g_terms[[i]] <- lapply(g_terms[[i]], reorder, varorder = basis)
    g_terms[[i]] <- lapply(g_terms[[i]], mpoly_to_stan)
    g_terms[[i]] <- unlist(g_terms[[i]])
    g_terms[[i]] <- c(g_terms[[i]])

    g[[i]] <- paste0(g_coef[[i]], "*", g_terms[[i]], collapse = "+")
    g[[i]] <- gsub("1\\*|\\*1", "", g[[i]])

    derivatives_pre[[i]] <- lapply(
      unique(unlist(vars)),
      get_derivative,
      num_of_vars = num_of_vars[i],
      deg = totaldeg[i],
      basis = unique(unlist(vars))
    )
  }
  g <- unlist(g)
  g <- paste0(
    "  vector[",
    length(g),
    "] g = [",
    paste(g, collapse = ","),
    "]';"
  )
  if (homo) {
    # homoskedastic case uses symbolic Jacobian entries
    for (i in seq_len(num_of_poly)) {
      derivatives[[i]] <- unlist(derivatives_pre[[i]])
    }
    # empty derivative expressions (from constant terms) become zero
    derivatives <- lapply(derivatives, function(v) {
      v[v == "**"] <- 0
      v
    })
    jac <- paste(
      sapply(
        derivatives,
        function(v) paste0("      [", paste(v, collapse = ","), "]")
      ),
      collapse = ",\n"
    )
  } else {
    # heteroskedastic path uses identity Jacobian surrogate
    n_eqs <- num_of_poly
    n_vars <- max(num_of_vars)
    jac <- array("", dim = c(n_eqs, n_vars))
    for (i in seq_len(n_eqs)) {
      for (j in seq_len(n_vars)) {
        jac[i, j] <- if (i == j) "1" else "0"
      }
    }
    jac <- apply(jac, 1L, paste, collapse = ", ")
    jac <- paste("      [", jac, "]", collapse = ", \n")
  }
  dg <- paste0(
    "  matrix[",
    num_of_poly,
    ",",
    max(num_of_vars),
    "] J = [ \n",
    jac,
    "\n    ];"
  )
  trans_block <- paste0("\ntransformed parameters {\n", g, "\n", dg, "\n}\n")

  # assemble full Stan program text
  model_block <- paste0(
    "\nmodel {\ntarget += normal_lpdf(0.00 | J' * ((J*J') \\ g), si);\n}"
  )
  stan_code <- paste0(data_block, params_block, trans_block, model_block)
  stan_code
}
