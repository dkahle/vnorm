# internal Stan code generation helper used by rvnorm()
# split out from rvnorm.R so coverage for rvnorm() can be assessed separately

create_stan_code <- function(poly, sd, n_eqs, w, homo, vars) {
  d <- get("deriv.mpoly", asNamespace("mpoly"))

  if (n_eqs == 1L) {
    if (!is.mpoly(poly)) poly <- mp(poly)
    if (missing(vars)) vars <- mpoly::vars(poly)
    n_vars <- length(vars)
    reorder.mpoly <- get("reorder.mpoly", asNamespace("mpoly"))
    poly <- reorder.mpoly(poly, varorder = sort(vars))

    g_string <- mpoly_to_stan(poly)
    if (homo) {
      grad <- if (n_vars > 1) {
        deriv(poly, var = mpoly::vars(poly))
      } else {
        gradient(poly)
      }

      # squared gradient norm for homoskedastic normalization
      ndg_sq <- if (n_vars >1)  Reduce(`+`, grad ^ 2) else grad^2
      ndg_string <- glue::glue("sqrt({mpoly_to_stan(ndg_sq)})")
    } else {
      ndg_string <- "1"
    }

    if (missing(w)) {
      parms <- glue::glue("real {vars};")
      parms <- paste(parms, collapse = "\n  ")
    } else if (is.numeric(w) && length(w) == 1L) {
      parms <- glue::glue("real<lower=-{w},upper={w}> {vars};")
      parms <- paste(parms, collapse = "\n  ")
    } else {
      parms <- sapply(seq_along(vars), function(i) {
        if (vars[i] %in% names(w)) {
          glue::glue(
            "real<lower={w[[vars[i]]][1]},upper={w[[vars[i]]][2]}> {vars[i]};"
          )
        } else {
          glue::glue("real {vars[i]};")
        }
      })
      parms <- paste(parms, collapse = "\n  ")
    }

    stan_code <- glue::glue(
      "
data {{
  real<lower=0> si;
}}

parameters {{
  {parms}
}}

transformed parameters {{
  real g = {g_string};
  real ndg = {ndg_string};
}}

model {{
  target += normal_lpdf(0.00 | g/ndg, si);
}}
    ")

  } else {
    vars <- mpoly::vars(poly)
    n_vars <- length(vars)
    printed_polys <- mpolyList_to_stan(poly)

    target_dim <- if (homo) n_vars else n_eqs
    data_string <- if (length(sd) == 1) {
      "real<lower=0> si"
    } else {
      paste0("cov_matrix[", target_dim, "] si")
    }
    model_string <- if (length(sd) == 1) {
      "normal_lpdf("
    } else {
      "multi_normal_lpdf("
    }
    mu_string <- if (length(sd) == 1) {
      "0.00"
    } else {
      paste0("[", paste(rep(0.00, target_dim), collapse = ","), "]'")
    }
    gbar_string <- "g"

    if (homo) {
      # Jacobian setup for homoskedastic pseudoinverse normalization
      printed_jac <- array("", dim = c(n_eqs, n_vars))
      for (i in seq_len(n_eqs)) {
        for (j in seq_len(n_vars)) {
          printed_jac[i, j] <- mpoly_to_stan(d(poly[[i]], vars[j]))
        }
      }

      # format Jacobian rows as Stan matrix literal
      printed_jac <- apply(printed_jac, 1, paste, collapse = ", ")
      printed_jac <- paste("      [", printed_jac, "]", collapse = ", \n")
      printed_jac <- paste0("[\n", printed_jac, "\n    ]")

      # pseudoinverse form depends on system shape (square, over, under)
      gbar_string <- if (n_vars == n_eqs) {
        "J \\ g"
      } else if (n_vars > n_eqs) {
        "J' * ((J*J') \\ g)"
      } else {
        "(J'*J) \\ (J'*g)"
      }
    }

    if (missing(w)) {
      parms <- paste(sprintf("real %s;", vars), collapse = "\n  ")
    } else if (is.numeric(w) && length(w) == 1L) {
      parms <- glue::glue("real<lower=-{w},upper={w}> {vars};")
      parms <- paste(parms, collapse = "\n  ")
    } else {
      parms <- sapply(seq_along(vars), function(i) {
        if (vars[i] %in% names(w)) {
          sprintf(
            "real<lower=%s,upper=%s> %s;",
            w[[vars[i]]][1],
            w[[vars[i]]][2],
            vars[i]
          )
        } else {
          sprintf("real %s;", vars[i])
        }
      })
      parms <- paste(parms, collapse = "\n  ")
    }
    transformed_parameters <- if (homo) {
      glue::glue(
        "  vector[{n_eqs}] g = [{printed_polys}]';
  matrix[{n_eqs},{n_vars}] J = {printed_jac};"
      )
    } else {
      glue::glue("  vector[{n_eqs}] g = [{printed_polys}]';")
    }

    stan_code <- glue::glue(
      "
data {{
  {data_string};
}}

parameters {{
  {parms}
}}

transformed parameters {{
{transformed_parameters}
}}

model {{
  target += {model_string}{mu_string} | {gbar_string}, si);
}}
      ")
  }
  stan_code
}
