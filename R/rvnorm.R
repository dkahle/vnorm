#' The Variety Normal Distribution
#'
#' Unnormalized density evaluation and random generation for the variety normal
#' distribution with mean defined by `poly` and scale `sd`. See Details for
#' caveats.
#'
#' If the variety you are interested in is connected, this strategy should work
#' well out of the box. If it is not, you will likely need multiple chains, and
#' sampling may become biased toward one or more components while down-sampling
#' others.
#' Question: what is the relative likelihood of each component, or an equal unit
#' of length, on different components? How does this generalize to more
#' varieties of varying dimensions?
#'
#' @param n Number of draws desired from each chain after warmup.
#' @param poly An `mpoly` object.
#' @param sd Scale parameter for the normal kernel. If `Sigma` is supplied,
#'   `sd` is replaced by `Sigma`.
#' @param Sigma Full covariance matrix or a diagonal vector of covariance terms.
#' @param output One of `"simple"`, `"tibble"`, or `"stanfit"`.
#' @param rejection If `TRUE`, rejection sampling is used.
#' @param chains The number of chains to run for the random number generation,
#'   see `stan()`.
#' @param cores The number of CPU cores to distribute the chains across, see
#'   `stan()`.
#' @param warmup Number of warmup iterations in `stan()`.
#' @param inc_warmup If `TRUE`, the MCMC warmup steps are included in the
#'   output.
#' @param thin `stan()` `thin` parameter.
#' @param homo If `TRUE`, sampling is from a homoskedastic variety normal
#'   distribution.
#' @param verbose If `TRUE`, print additional progress messages.
#' @param vars Character vector of polynomial indeterminates.
#' @param numerator,denominator Character scalars containing printed numerator
#'   and denominator forms.
#' @param w A named list of box constraints for vectors to be passed to Stan,
#'   see examples. If a single number, a box window `(-w, w)` is applied to all
#'   variables.
#' @param refresh The `refresh` argument of `stan()`, which governs how
#'   much information is provided to the user while sampling.
#' @param pre_compiled Whether to use precompiled Stan models. Available
#'   for polynomials with three indeterminates and three degrees. Defaults to
#'   `TRUE`.
#' @param code_only If `TRUE`, only formulate and return Stan code.
#' @param ... Additional parameters passed to `stan()`.
#' @param user_compiled If `TRUE`, use a user-compiled Stan program produced by
#'   [compile_stan_code()]. Defaults to `FALSE`.
#' @param show_messages If `TRUE`, Stan sampler messages are shown.
#' @param seed Optional integer seed passed to CmdStan's sampler for
#'   reproducibility. Note that R's `set.seed()` does not control Stan's RNG.
#' @name rvnorm
#' @return Either (1) matrix whose rows are the individual draws from the
#'   distribution, (2) a [tbl_df-class] object with the draws along with
#'   additional information, or (3) an object of class `stanfit`.
#' @examples
#'
#'
#' library("tidyverse")
#' options("mc.cores" = parallel::detectCores() - 1)
#'
#' \dontrun{
#' ## basic usage
#' ########################################
#'
#' # single polynomial
#' p <- mp("x^2 + y^2 - 1")
#' samps <- rvnorm(1000, p, sd = .05)
#' head(samps)
#' str(samps)
#' plot(samps, asp = 1)
#'
#' # returning a data frame
#' (samps <- rvnorm(5000, p, sd = .05, w = 2, output = "tibble"))
#'
#' ggplot(samps, aes(x, y)) +
#'   geom_point(size = .5) +
#'   coord_equal()
#'
#' ggplot(samps, aes(x, y)) +
#'   geom_point(size = .5) +
#'   geom_variety(poly = p) +
#'   coord_equal()
#'
#' ggplot(samps, aes(x, y)) +
#'   geom_point(aes(color = factor(.chain)), size = .5) +
#'   geom_variety(poly = p, linewidth = 1) +
#'   coord_equal()
#'
#' ggplot(samps, aes(x, y)) +
#'   stat_density2d(
#'     aes(fill = after_stat(density)),
#'     geom = "raster", contour = FALSE
#'    ) +
#'   geom_variety(poly = p) +
#'   coord_equal()
#'
#' library("ggdensity")
#' ggplot(samps, aes(x, y)) +
#'   geom_hdr(xlim = c(-2,2), ylim = c(-2,2)) +
#'   geom_variety(poly = p) +
#'   coord_equal()
#'
#' # in three variables
#' (samps <- rvnorm(20, mp("x^2 + y^2 + z^2 - 1"), sd = .05, w = 2))
#' apply(samps, 1, function(v) sqrt(sum(v^2)))
#'
#' # more than one polynomial, # vars > # eqns, underdetermined system
#' p <- mp(c("x^2 + y^2 + z^2 - 1", "z"))
#' (samps <- rvnorm(500, p, sd = .1, output = "tibble"))
#' ggplot(samps, aes(x, y)) + geom_point(size = .5) + coord_equal()
#'
#' ggplot(samps, aes(x, y, color = `g[1]`)) + geom_point() +
#'   scale_color_gradient2(mid = "gray80") + coord_equal()
#'
#' ggplot(samps, aes(x, y, color = `g[2]`)) + geom_point() +
#'   scale_color_gradient2(mid = "gray80") + coord_equal()
#'
#' ggplot(samps, aes(x, z, color = `g[1]`)) + geom_point() +
#'   scale_color_gradient2(mid = "gray80") + coord_equal()
#'
#' # more than one polynomial, # vars < # eqns, overdetermined system
#' p <- mp(c("3 x", "3 y", "2 x + 2 y", "3 (x^2 + y)", "3 (x^2 - y)"))
#' (samps <- rvnorm(500, p, sd = .1, output = "tibble"))
#'
#' samps |>
#'   select(x, y, starts_with("g")) |>
#'   pivot_longer(
#'     starts_with("g"),
#'     names_to = "equation",
#'     values_to = "value"
#'   ) |>
#'   ggplot(aes(x, y, color = value)) + geom_point() +
#'     scale_color_gradient2(mid = "gray80") + coord_equal() +
#'     facet_wrap(~ equation)
#'
#' ## using refresh to get more info
#' ########################################
#'
#' rvnorm(2000, p, sd = .1, "tibble", verbose = TRUE)
#' rvnorm(2000, p, sd = .1, "tibble", refresh = 500)
#' rvnorm(2000, p, sd = .1, "tibble", refresh = 0) # default
#' rvnorm(2000, p, sd = .1, "tibble", refresh = -1)
#'
#' ## many chains in parallel
#' ########################################
#'
#' options(mc.cores = parallel::detectCores())
#' p <- mp("x^2 + (4 y)^2 - 1")
#' (samps <- rvnorm(1e4, p, sd = .01, "tibble", verbose = TRUE, chains = 8))
#' ggplot(samps, aes(x, y)) + geom_bin2d(binwidth = .01*c(1,1)) + coord_equal()
#' # decrease sd to get more uniform sampling
#'
#' ## windowing for unbounded varieties
#' ########################################
#' # windowing is needed for unbounded varieties
#' # in the following, look at the parameters block
#'
#' p <- mp("x y - 1") # unbounded variety, 1 poly
#' p <- mp(c("x y - 1", "y - x")) # 2 polys
#'
#' rvnorm(1e3, p, sd = .01, "tibble", code_only = TRUE)
#' rvnorm(1e3, p, sd = .01, "tibble",  w = 1.15)
#'
#' window <- list("x" = c(-1.5, 1.5))
#' rvnorm(1e3, p, sd = .01, "tibble",  w = window)
#'
#' ## the importance of normalizing
#' ########################################
#' # one of the effects of the normalizing is to stabilize variances, making
#' # them roughly equivalent globally over the variety.
#'
#' # lemniscate of bernoulli
#' p <- mp("(x^2 + y^2)^2 - 2 (x^2 - y^2)")
#'
#' # normalized, good
#' (samps <- rvnorm(2000, p, .05, "tibble"))
#' ggplot(samps, aes(x, y)) + geom_point(size = .5) + coord_equal()
#' ggplot(samps, aes(x, y)) + geom_bin2d(binwidth = .05*c(1,1)) + coord_equal()
#'
#' # unnormalized, bad
#' (samps <- rvnorm(2000, p, .05, "tibble", homo = FALSE))
#' ggplot(samps, aes(x, y)) + geom_point(size = .5) + coord_equal()
#' ggplot(samps, aes(x, y)) + geom_bin2d(binwidth = .05*c(1,1)) + coord_equal()
#'
#' ## semi-algebraic sets
#' ########################################
#' # inside the semialgebraic set x^2 + y^2 <= 1
#' # this is the same as x^2 + y^2 - 1 <= 0, so that
#' # x^2 + y^2 - 1 + s^2 == 0 for some slack variable s
#' # this is the projection of the sphere into the xy-plane.
#'
#' p <- mp("1 - (x^2 + y^2) - s^2")
#' samps <- rvnorm(1e4, p, sd = .1, "tibble", chains = 8, refresh = 1e3)
#' ggplot(samps, aes(x, y)) + geom_bin2d(binwidth = .05*c(1,1)) + coord_equal()
#'
#' ggplot(sample_n(samps, 2e3), aes(x, y, color = s)) +
#'   geom_point(size = .5) +
#'   scale_color_gradient2() +
#'   coord_equal()
#'
#' # alternative representation
#' # x^2 + y^2 - 1 <= 0 iff s^2 (x^2 + y^2 - 1) + 1 == 0
#' # note that it's gradient is more complicated.
#' p <- mp("s^2 (x^2 + y^2 - 1) + 1")
#' samps <- rvnorm(1e4, p, sd = .1, "tibble", chains = 8, w = 2, refresh = 1e3)
#' ggplot(samps, aes(x, y)) + geom_bin2d(binwidth = .05*c(1,1)) + coord_equal()
#'
#' ## keeping the warmup / the importance of multiple chains
#' ########################################
#'
#' p <- mp("((x + 1.5)^2 + y^2 - 1) ((x - 1.5)^2 + y^2 - 1)")
#' ggplot() +
#'   geom_variety(poly = p, xlim = c(-3,3)) +
#'   theme(legend.position = "top") +
#'   coord_equal()
#'
#' # notice the migration of chains initialized away from the distribution
#' # (it helps to make the graphic large on your screen)
#' samps <- rvnorm(500, p, sd = .05, "tibble", chains = 8, inc_warmup = TRUE)
#'
#' ## ideal-variety correspondence considerations
#' ########################################
#'
#' p <- mp("x^2 + y^2 - 1")
#'
#' samps_1 <- rvnorm(250, p^1, sd = .1, output = "tibble", chains = 8)
#' samps_2 <- rvnorm(250, p^2, sd = .1, output = "tibble", chains = 8)
#' samps_3 <- rvnorm(250, p^3, sd = .1, output = "tibble", chains = 8)
#' samps_4 <- rvnorm(250, p^4, sd = .1, output = "tibble", chains = 8)
#' samps <- bind_rows(mget(apropos("samps_[1-4]")))
#' samps$power <- rep(seq_along(apropos("samps_[1-4]")), each = 250)
#'
#' ggplot(samps, aes(x, y, color = g < 0)) +
#'   geom_point(size = .5) +
#'   coord_equal(xlim = c(-3,3), ylim = c(-3,3)) +
#'   facet_wrap(~ power)
#'
#' ## neat examples
#' ########################################
#' # an implicit Lissajous region, view in separate window large
#'
#' # x = cos(m t + p)
#' # y = sin(n t + q)
#' (p <- lissajous(3, 2,  -pi/2, 0))
#' (p <- lissajous(4, 3,  -pi/2, 0))
#' (p <- lissajous(5, 4,  -pi/2, 0))
#' (p <- lissajous(3, 3,  0, 0))
#' (p <- lissajous(5, 5,  0, 0))
#' (p <- lissajous(7, 7,  0, 0))
#' ggplot() +
#'   geom_variety(poly = p, xlim = c(-2, 2), ylim = c(-2, 2), n = 201,
#'   show.legend = FALSE) +
#'   coord_equal()
#'
#' p <- plug(p, "x", mp(".5 x"))
#' p <- plug(p, "y", mp(".5 y"))
#'
#' # algebraic set
#' samps <- rvnorm(5e3, p, sd = .01, "tibble", chains = 8, cores = 8)
#' ggplot(samps, aes(x, y, color = factor(.chain))) +
#'   geom_point(size = .5) +
#'   coord_equal()
#'
#' ggplot(samps, aes(x, y)) + geom_hdr() + coord_equal()
#'
#' ggplot(samps, aes(x, y, color = factor(.chain))) +
#'   geom_point(size = .5) +
#'   coord_equal() +
#'   facet_wrap(~ factor(.chain))
#'
#' # semi-algebraic set
#' samps_normd <- rvnorm(1e4, p + mp("s^2"), sd = .01, "tibble", chains = 8,
#'   cores = 8, homo = TRUE
#' )
#' samps_unormd <- rvnorm(1e4, p + mp("s^2"), sd = .01, "tibble", chains = 8,
#'   cores = 8, homo = FALSE
#' )
#'
#' bind_rows(
#'   samps_normd  |> mutate(normd = TRUE),
#'   samps_unormd |> mutate(normd = FALSE)
#' ) |>
#'   ggplot(aes(x, y)) +
#'     geom_point(size = .5) +
#'     # geom_bin2d(binwidth = .05*c(1,1)) +
#'     facet_grid(normd ~ .chain) +
#'     coord_equal()
#'
#' ggplot(samps_normd, aes(x, y)) +
#'   geom_bin2d(binwidth = .05*c(1,1)) + coord_equal()
#' }
#' @rdname rvnorm
#' @export
rvnorm <- function(
  n, poly, sd, output = "simple", Sigma = NULL, rejection = FALSE,
  chains = 4L, warmup = max(500, floor(n / 2)), inc_warmup = FALSE,
  thin = 1L, verbose = FALSE, cores = min(chains, getOption("mc.cores", 1L)),
  homo = TRUE, w, vars, numerator, denominator, refresh = 0L,
  code_only = FALSE, pre_compiled = TRUE, user_compiled = FALSE,
  show_messages = FALSE, seed = NULL, ...
) {

  if (!is.numeric(n) || length(n) != 1L || n < 1 || n != as.integer(n)) {
    stop("`n` must be a positive integer.", call. = FALSE)
  }
  n <- as.integer(n)

  if (rejection) {
    if (missing(w)) {
      w <- 1
    }
    if (!(is.numeric(sd) && is.vector(sd) && length(sd) == 1L)) {
      stop("This sd is not supported yet for rejection sampler.")
    }
    sample <- rejection_sampler(
      n = n,
      poly = poly,
      sd = sd,
      vars = sort(mpoly::vars(poly)),
      w = w,
      output = output,
      homo = homo,
      message = verbose
    )
    return(sample)
  }

  output_needs_rewriting <- FALSE

  if (is.character(poly)) {
    poly <- mp(poly)
  }

  if (is.mpoly(poly)) {
    n_eqs <- 1L
  } else if (is.mpolyList(poly)) {
    n_eqs <- length(poly)
  } else {
    stop(
      "`poly` should be either a character vector, mpoly, or mpolyList.",
      call. = FALSE
    )
  }
  # resolve sd/Sigma into a single value for Stan
  sigma_stan <- if (!is.null(Sigma)) Sigma else sd
  n_vars <- length(mpoly::vars(poly))
  if (is.numeric(sigma_stan) && length(sigma_stan) == 1) {
  } else if (is.vector(sigma_stan) && length(sigma_stan) == n_vars) {
    sigma_stan <- diag(sigma_stan)
  } else if (is.matrix(sigma_stan) &&
    nrow(sigma_stan) == n_vars && ncol(sigma_stan) == n_vars) {
  } else {
    stop(
      "`Sigma` should be a scalar, a vector of length equal to the number of ",
      "variables, or a square matrix with dimensions equal to the number of ",
      "variables.",
      call. = FALSE
    )
  }

  if (isTRUE(verbose) && refresh == 0L) {
    refresh <- max(ceiling(n / 10), 1L)
  }

  if (
    n_eqs > 1 ||
      length(mpoly::vars(poly)) > 3 ||
      base::max(mpoly::totaldeg(poly)) > 3
  ) {
    # precompiled binaries only cover a small catalog of low-degree templates
    pre_compiled <- FALSE
  }

  if (code_only) {
    stan_code <- create_stan_code(poly, sigma_stan, n_eqs, w, homo)
    return(stan_code)
  }

  if (user_compiled) {
    # look up a previously compiled user template from the internal cache
    model_name <- generate_model_name(
      poly = poly, homo = homo, windowed = !missing(w)
    )
    compiled_stan_info <- get_compiled_stan_info()
    if (nrow(compiled_stan_info) == 0) {
      stop(
        "No compiled model cache found. Run compile_stan_code() first.",
        call. = FALSE
      )
    }

    model_path <- compiled_stan_info$path[compiled_stan_info$name == model_name]
    if (length(model_path) == 0) {
      stop(
        "Requested compiled model not found in cache. ",
        "Run compile_stan_code() for this polynomial first.",
        call. = FALSE
      )
    }

    model <- cmdstan_model(model_path[1])
    stan_data <- get_coefficients_data(poly)
    stan_data <- if (missing(w)) {
      c(stan_data, "si" = sigma_stan)
    } else {
      c(stan_data, "w" = w, "si" = sigma_stan)
    }
  } else if (pre_compiled) {
    poly_original <- poly
    output_needs_rewriting_original <- output_needs_rewriting
    # remap variable names to x/y/z because shipped models are keyed that way
    list_for_transformation <- check_and_replace_vars(poly_original)
    if (!is.null(unlist(list_for_transformation$mapping))) {
      poly <- list_for_transformation$polynomial
      var_info <- list_for_transformation$mapping
      output_needs_rewriting <- TRUE
    }

    deg <- mpoly::totaldeg(poly)
    num_of_vars <- length(mpoly::vars(poly))
    stan_file_name <- paste(
      num_of_vars,
      deg,
      if (homo) "vn" else "hvn",
      sep = "_"
    )
    if (!missing(w)) {
      stan_file_name <- paste(stan_file_name, "w", sep = "_")
    }

    model <- tryCatch(
      instantiate::stan_package_model(
        name = stan_file_name,
        package = "vnorm"
      ),
      error = function(e) {
        message(
          "pre-compiled Stan model not available; ",
          "falling back to regular rvnorm sampling"
        )
        NULL
      }
    )
    if (!model_has_usable_executable(model)) {
      poly <- poly_original
      output_needs_rewriting <- output_needs_rewriting_original
      stan_code <- create_stan_code(poly, sigma_stan, n_eqs, w, homo, vars)
      stan_file <- write_stan_file(stan_code)
      model <- cmdstan_model(stan_file)
      stan_data <- list("si" = sigma_stan)
    } else {
      stan_data <- make_coefficients_data(poly, num_of_vars, deg)
      stan_data <- if (missing(w)) {
        c(stan_data, "si" = sigma_stan)
      } else {
        c(stan_data, "w" = w, "si" = sigma_stan)
      }
    }
  } else {
    # fall back to generating and compiling a temporary Stan model
    stan_code <- create_stan_code(poly, sigma_stan, n_eqs, w, homo, vars)
    stan_file <- write_stan_file(stan_code)
    model <- cmdstan_model(stan_file)
    stan_data <- list("si" = sigma_stan)
  }

  sample_args <- list(
    data = stan_data,
    refresh = refresh,
    iter_warmup = warmup,
    iter_sampling = ceiling(n / chains),
    chains = chains,
    parallel_chains = cores,
    adapt_delta = .999,
    max_treedepth = 20L,
    save_warmup = inc_warmup,
    show_messages = show_messages,
    ...
  )
  # seed is injected separately to avoid passing NULL to cmdStan
  if (!is.null(seed)) sample_args$seed <- seed
  samps <- do.call(model$sample, sample_args)

  if (output == "simple") {
    df <- as.data.frame(samps$draws(format = "df", inc_warmup = inc_warmup))
    if (!inc_warmup) df <- tail(df, n)
    df <- df[, mpoly::vars(poly), drop = FALSE]
    row.names(df) <- NULL
    if (output_needs_rewriting) {
      df <- rename_output_df(df, replacement_list = var_info)
    }
    return(df)
  }

  if (output == "tibble") {
    df <- as.data.frame(samps$draws(format = "df", inc_warmup = inc_warmup))
    if (!inc_warmup) df <- tail(df, n)
    row.names(df) <- NULL
    df <- cbind(df[, -1], df[, 1])
    if (output_needs_rewriting) {
      df <- rename_output_df(df, replacement_list = var_info)
    }
    return(tibble::as_tibble(df))
  }

  if (output == "stanfit") {
    if (output_needs_rewriting) {
      message(
        "The stanfit object has variable names as x,y,z instead of the ones ",
        "on the polynomial. Don't use pre-compiled model if you need output ",
        "names to match your polynomial."
      )
    }
    return(samps)
  }
}

model_has_usable_executable <- function(model) {
  if (is.null(model)) return(FALSE)

  exe_ref <- model$exe_file
  if (is.null(exe_ref)) {
    # Test doubles may only provide a sample() method; treat those as usable.
    return(is.function(model$sample))
  }

  exe_path <- tryCatch(
    if (is.function(exe_ref)) exe_ref() else exe_ref,
    error = function(e) NULL
  )

  is.character(exe_path) &&
    length(exe_path) >= 1 &&
    nzchar(exe_path[[1]]) &&
    file.exists(exe_path[[1]])
}
