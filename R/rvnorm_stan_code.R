# internal Stan code generation helper used by rvnorm()
# split out from rvnorm.R so coverage for rvnorm() can be assessed separately

create_stan_code <- function(poly, sd, n_eqs, w, homo, vars) {
  if (n_eqs == 1L) {
    if (!is.mpoly(poly)) poly <- mp(poly)
    if (missing(vars)) vars <- mpoly::vars(poly)
    # the gradient runs over the variables poly involves; any others in
    # `vars` add nothing to it
    n_vars <- length(mpoly::vars(poly))
    poly <- reorder(poly, varorder = sort(vars))

    g_string <- stan_poly_expr(poly)
    if (homo) {
      grad <- if (n_vars > 1) {
        deriv(poly, var = mpoly::vars(poly))
      } else {
        gradient(poly)
      }

      # squared gradient norm for homoskedastic normalization
      ndg_sq <- if (n_vars >1)  Reduce(`+`, grad ^ 2) else grad^2
      ndg_string <- glue::glue("sqrt({stan_poly_expr(ndg_sq)})")
    } else {
      ndg_string <- "1"
    }

    if (missing(w)) {
      parms <- glue::glue("real {vars};")
      parms <- paste(parms, collapse = "\n  ")
    } else if (is.numeric(w) && length(w) == 1L) {
      parms <- glue::glue("real<lower={stan_real(-w)},upper={stan_real(w)}> {vars};")
      parms <- paste(parms, collapse = "\n  ")
    } else {
      parms <- sapply(seq_along(vars), function(i) {
        if (vars[i] %in% names(w)) {
          glue::glue(
            "real<lower={stan_real(w[[vars[i]]][1])},upper={stan_real(w[[vars[i]]][2])}> {vars[i]};"
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
    # sorted order fixes the parameter vector and Jacobian columns, so a
    # vector sd / Sigma over variables (homo = TRUE) means what it does in
    # pdvnorm() regardless of the order variables appear in the terms
    poly_vars <- sort(mpoly::vars(poly))
    n_vars <- length(poly_vars)
    # parameters for every sampled variable; those the equations do not
    # involve stay out of J
    if (missing(vars)) vars <- poly_vars
    printed_polys <- paste(vapply(poly, stan_poly_expr, character(1)), collapse = ", ")

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
      # jacobian setup for homoskedastic pseudoinverse normalization
      printed_jac <- array("", dim = c(n_eqs, n_vars))
      for (i in seq_len(n_eqs)) {
        for (j in seq_len(n_vars)) {
          printed_jac[i, j] <- stan_poly_expr(deriv(poly[[i]], poly_vars[j]))
        }
      }

      # format jacobian rows as Stan matrix literal
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
      parms <- glue::glue("real<lower={stan_real(-w)},upper={stan_real(w)}> {vars};")
      parms <- paste(parms, collapse = "\n  ")
    } else {
      parms <- sapply(seq_along(vars), function(i) {
        if (vars[i] %in% names(w)) {
          sprintf(
            "real<lower=%s,upper=%s> %s;",
            stan_real(w[[vars[i]]][1]),
            stan_real(w[[vars[i]]][2]),
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
  stan_real_literals(stan_code)
}

stan_real_literals <- function(code) {
  # Stan integers are 32-bit, so integer-valued coefficients printed without a
  # decimal point fail to compile past 2^31 - 1; write long ones as reals
  gsub(
    "(?<![[:alnum:]_.])([0-9]{10,})(?![[:alnum:]_.])", "\\1.0", code,
    perl = TRUE
  )
}

stan_poly_expr <- function(poly) {
  # print a polynomial as a Stan expression with its exact coefficients;
  # mpoly's printer rounds them to getOption("digits") and prints any within
  # sqrt(eps) of 1 as 1. Same layout as mpoly_to_stan() otherwise
  terms <- vapply(unclass(poly), function(term) {
    powers <- term[names(term) != "coef"]
    coef <- unname(term[["coef"]])
    if (length(powers) == 0L) return(stan_real(coef))
    factors <- paste0(names(powers), ifelse(powers == 1, "", paste0("^", powers)))
    lead <- if (coef == 1) "" else if (coef == -1) "-" else paste0(stan_real(coef), "*")
    paste0(lead, paste(factors, collapse = "*"))
  }, character(1))
  gsub("+-", "-", paste(terms, collapse = "+"), fixed = TRUE)
}

stan_real <- function(x) {
  # shortest decimal that stanc and C++ read back as x: 15 significant digits
  # when they do, else 17, which always do. Both round correctly, as C's
  # strtod() does under jsonlite's parser; R's own parser does not on every
  # platform (on arm64 it reads -0.000000000261404920558413 as
  # -0x1.1f6af1a526408p-32, one ulp from the correctly rounded value)
  out <- format(x, digits = 15, scientific = FALSE)
  if (jsonlite::parse_json(out) != x) out <- format(x, digits = 17, scientific = FALSE)
  out
}
