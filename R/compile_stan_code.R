#' Compile a Stan Model for User-Defined Polynomials
#'
#' Compile and cache a Stan model template for a polynomial shape so repeated
#' calls to [rvnorm()] can reuse the compiled binary.
#'
#' @param poly An `mpoly` or `mpolyList` object. An `mpolyList` of one
#'   polynomial is treated as that polynomial, as in [rvnorm()].
#' @param custom_stan_code If `TRUE`, a custom model is compiled even if the
#'   general case of the polynomial is already included during package
#'   installation. Defaults to `FALSE`.
#' @param windowed If `TRUE`, the compiled Stan model includes a window (box
#'   constraint) parameter `w` in its data block. Defaults to `FALSE`.
#' @param homo If `TRUE`, sampling is done from a homoskedastic variety
#'   normal distribution. Defaults to `TRUE`.
#' @details The compiled model metadata is stored in an internal package cache
#'   used by [rvnorm()] with `user_compiled = TRUE`.
#' @return The compiled `CmdStanModel`.
#'
#' @export
#' @examples
#'
#' \dontrun{
#' # compile a model that looks like b0 + bx6 x^6 + by6 y^6 for later input
#' p <- mp("x^6 + y^6 - 1") # template polynomial
#' samps <- rvnorm(1000, p, sd = .05)
#' head(samps)
#' compile_stan_code(p) # allows to change coefficients
#' p <- mp("x^6 + 8 y^6 - 1")
#' rvnorm(1e4, p, .05, user_compiled = TRUE)
#'
#' }
#'
compile_stan_code <- function(
    poly, custom_stan_code = FALSE, windowed = FALSE, homo = TRUE
  ) {
  if (!(is.mpoly(poly) || is.mpolyList(poly))) {
    stop("`poly` should be an mpoly or mpolyList object.", call. = FALSE)
  }
  # rvnorm() samples an mpolyList of one polynomial as that polynomial, so its
  # model must be cached under the polynomial's name
  if (is.mpolyList(poly) && length(poly) == 1L) {
    poly <- poly[[1]]
  }

  if (!custom_stan_code && is.mpoly(poly)) {
    # reuse shipped templates when available unless user forces a custom
    # compile
    if (length(mpoly::vars(poly)) < 4 && base::max(mpoly::totaldeg(poly)) < 4) {
      stop(
        paste0(
          "Pre-compiled model for the general case already ",
          "exists from installation. "
        ),
        "Use custom_stan_code = TRUE to use a custom model anyway.",
        call. = FALSE
      )
    }
  }

  stan_code <- get_custom_stan_code(poly = poly, windowed = windowed, homo = homo)
  model_name <- generate_model_name(poly = poly, windowed = windowed, homo = homo)
  # write source to a temp .Stan file and cache the mapping for rvnorm()
  model_path <- cmdstanr::write_stan_file(stan_code, dir = tempdir())
  info_before <- get_compiled_stan_info()
  add_compiled_stan_info(name = model_name, path = model_path)
  info_after <- get_compiled_stan_info()
  if (nrow(info_before) == 0) {
    message("compiled_stan_info cache created in package namespace")
  } else if (nrow(info_after) > nrow(info_before)) {
    message("compiled_stan_info cache updated in package namespace")
  } else {
    message("compiled_stan_info cache already contained this model")
  }

  cmdstan_model(model_path)
}

get_custom_stan_code <- function(poly, windowed = FALSE, homo = TRUE) {
  if (is.mpoly(poly)) {
    # single-polynomial program: scalar g and scalar normalized distance
    poly <- canonicalize_mpoly(poly)
    vars <- sort(mpoly::vars(poly))
    num_of_vars <- length(vars)

    # data block: lifted coefficients (+ optional box width w)
    var_for_data_block <- names(stan_coefficients(poly, varorder = vars))

    data_block <- paste(
      sapply(var_for_data_block, function(x) paste("  real", x)),
      collapse = "; "
    )
    data_block <- paste0(data_block, ";")
    if (windowed) {
      data_block <- paste0(data_block, "\n  real w;")
    }
    data_block <- paste0("data {\n  real<lower=0> si;\n", data_block, "\n}\n")

    # parameter block: unconstrained or box-constrained coordinates
    if (windowed) {
      params_block <- paste(sapply(vars, function(var) {
        paste0("  real<lower=-", "w", ", upper=", "w", "> ", var, ";")
      }), collapse = "\n")
    } else {
      params_block <- paste(sapply(vars, function(var) {
        paste0("  real ", var, ";")
      }), collapse = "\n")
    }
    params_block <- paste0("parameters {\n", params_block, "\n }\n")

    g <- stan_coef_lift(poly, varorder = vars)
    g <- mpoly_to_stan(g)

    derivatives <- lapply(
      vars,
      helper_for_derivative_for_mpoly_stan_code,
      poly = poly,
      varorder = vars
    )

    derivative_names <- sapply(seq_along(vars), function(i) {
      paste0("dg", vars[i])
    })

    ndg <- if (homo) {
      paste0(
        "  real ndg = sqrt(",
        paste0(derivative_names, "^2", collapse = " + "),
        ");"
      )
    } else {
      "  real ndg = 1;"
    }

    dg <- sapply(seq_along(vars), function(i) {
      paste(
        derivative_names[i],
        paste(derivatives[[i]], collapse = ""),
        sep = " = "
      )
    })
    dg <- paste0("  real ", dg)
    dg <- paste(dg, collapse = ";")

    model_block <- paste0(
      "model {\n  real g = ", g, ";\n", dg, ";\n", ndg,
      "\n  target += normal_lpdf(0.00 | g/ndg, si); \n}"
    )

    stan_code <- paste0(data_block, params_block, model_block, sep = "")
  } else if (is.mpolyList(poly)) {
    # multi-polynomial program: vector g and matrix jacobian J
    poly <- canonicalize_mpolylist(poly)
    poly <- sort_mpolylist_lexicographically(poly)
    n_eqs <- length(poly)
    n_vars <- length(mpoly::vars(poly))

    vars_by_poly <- vector("list", length(poly))
    for (i in seq_along(poly)) {
      vars_by_poly[[i]] <- sort(mpoly::vars(poly[[i]]))
    }

    # extract and suffix coefficients per polynomial for Stan data block
    var_for_data_block <- vector("list", length(poly))
    for (i in seq_along(poly)) {
      var_for_data_block[[i]] <- names(stan_coefficients(
        poly[[i]],
        varorder = vars_by_poly[[i]],
        suffix = i
      ))
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

    data_block <- paste0("data {\n  real<lower=0> si;\n", data_block, "\n}\n")
    vars_for_params <- sort(unique(unlist(vars_by_poly)))
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

    g <- vector("list", length(poly))
    derivatives_pre <- vector("list", length(poly))
    derivatives <- vector("list", length(poly))

    for (i in seq_along(poly)) {
      g[[i]] <- helper_for_coef_lift_for_mpolylist(
        poly[[i]],
        i,
        varorder = vars_by_poly[[i]]
      )
      g[[i]] <- mpoly_to_stan(g[[i]])
      derivatives_pre[[i]] <- lapply(
        sort(mpoly::vars(poly)),
        helper_for_derivative_for_mpolylist_stan_code,
        poly = poly[[i]],
        i = i,
        varorder = vars_by_poly[[i]]
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

    gbar_string <- "g"
    if (homo) {
      for (i in seq_along(poly)) {
        derivatives[[i]] <- unlist(derivatives_pre[[i]])
      }
      jac <- paste(
        sapply(
          derivatives,
          function(v) paste0("      [", paste(v, collapse = ","), "]")
        ),
        collapse = ",\n"
      )
      dg <- paste0("  matrix[", n_eqs, ",", n_vars, "] J = [ \n", jac, "\n    ];")
      trans_block <- paste0("\ntransformed parameters {\n", g, "\n", dg, "\n}\n")

      gbar_string <- if (n_vars == n_eqs) {
        "J \\ g"
      } else if (n_vars > n_eqs) {
        "J' * ((J*J') \\ g)"
      } else {
        "(J'*J) \\ (J'*g)"
      }
    } else {
      trans_block <- paste0("\ntransformed parameters {\n", g, "\n}\n")
    }

    model_block <- paste0(
      "\nmodel {\n  target += normal_lpdf(0.00 |",
      gbar_string,
      ", si);\n}"
    )
    stan_code <- paste0(data_block, params_block, trans_block, model_block)
  } else {
    stop(
      "`poly` should either be an mpoly or an mpolyList object.",
      call. = FALSE
    )
  }

  stan_code
}
