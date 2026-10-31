#' The Variety Normal Distribution
#'
#' Random generation for the variety normal distribution with mean defined by
#' `poly` and scale `sd`; see [pdvnorm()] for density evaluation. See Details
#' for caveats.
#'
#' @param n Total number of post-warmup draws to return (after any thinning);
#'   they are split evenly across `chains`.
#' @param poly An `mpoly` or `mpolyList` object, or a character vector that
#'   [mpoly::mp()] can parse. An `mpolyList` of one polynomial is treated as
#'   that polynomial.
#' @param sd Scale parameter for the normal kernel. If `Sigma` is supplied,
#'   `sd` is replaced by `Sigma`. For single-polynomial sampling, `sd` must be
#'   a positive scalar. For polynomial-list sampling, `sd` may be a positive
#'   scalar or a positive vector matching the Stan target dimension. With
#'   `noise = "domain"`, `sd` gives the standard deviations of the noise in the
#'   space of the variables: a scalar, or one per variable in `vars`.
#' @param Sigma Covariance parameter. For single-polynomial sampling, only a
#'   positive scalar variance is supported. For polynomial-list sampling,
#'   `Sigma` may be a positive scalar variance, a positive diagonal vector, or a
#'   positive-definite covariance matrix matching the Stan target dimension.
#'   With `noise = "domain"`, `Sigma` is the covariance of the noise in the
#'   space of the variables: a scalar variance, a vector of variances, or a
#'   covariance matrix over the variables in `vars`, in sorted order unless
#'   its entries are named.
#' @param output One of `"simple"`, `"tibble"`, or `"stanfit"`.
#' @param rejection If `TRUE`, rejection sampling is used.
#' @param chains The number of Markov chains to run.
#' @param cores The number of chains to run in parallel (CmdStan's
#'   `parallel_chains`).
#' @param warmup Number of warmup iterations per chain. By default this is
#'   conservative and scales with the number of post-warmup draws per chain.
#' @param inc_warmup If `TRUE`, the MCMC warmup steps are included in the
#'   output.
#' @param thin Thinning interval passed to CmdStan. Sampling runs `thin` times
#'   as many iterations, so `n` draws are still returned.
#' @param adapt_delta,max_treedepth Conservative CmdStan sampling defaults.
#'   Pass smaller values only when you are comfortable with the resulting
#'   diagnostics.
#' @param homo If `TRUE`, sampling is from a homoskedastic variety normal
#'   distribution.
#' @param noise Where the noise lives. `"codomain"` (the default) gives the
#'   variety normal distribution, a normal kernel on the polynomial's values
#'   (see `homo`). `"domain"` gives points nearly uniform on the variety plus
#'   multivariate normal noise in the space of the variables, with covariance
#'   from `sd` or `Sigma`; see Details.
#' @param pilot_sd With `noise = "domain"`, the scale of the homogeneous
#'   variety normal whose draws are projected onto the variety to give its
#'   uniform points. Where the variety is smooth and `poly` has no repeated
#'   factors, their distribution differs from uniform by terms of order
#'   `pilot_sd^2` times the variety's squared curvature; the cost of
#'   sampling grows roughly like `1 / pilot_sd`. Defaults to a tenth of the
#'   noise scale, `sqrt(mean(diag(Sigma))) / 10`.
#' @param verbose If `TRUE`, print additional progress messages.
#' @param vars Character vector of the variables to sample. Defaults to the
#'   variables of `poly`, and must include all of them. The distribution is
#'   flat along variables that `poly` does not involve, so `w` must bound
#'   those, and they are sampled uniformly within it.
#' @param numerator,denominator Character scalars containing printed numerator
#'   and denominator forms.
#' @param w A named list of box constraints for vectors to be passed to Stan,
#'   see examples. If a single number, a box window `(-w, w)` is applied to all
#'   variables. With `rejection = TRUE`, a missing `w` means `w = 1.25`, the
#'   [rejection_sampler()] default.
#' @param refresh CmdStan's `refresh`, which governs how often sampling
#'   progress is reported.
#' @param pre_compiled If `TRUE` (the default), sample with a Stan model
#'   compiled when vnorm was installed whenever one covers the request: a
#'   single polynomial in at most three variables (counting `vars`) of total
#'   degree at most three, with no window or the same window `(-a, a)` on
#'   every variable. Other requests, and installations without these models,
#'   compile a model for `poly`.
#' @param code_only If `TRUE`, only formulate and return Stan code.
#' @param ... Additional arguments passed to the CmdStan model's `$sample()`
#'   method.
#' @param user_compiled If `TRUE`, use a user-compiled Stan program produced by
#'   [compile_stan_code()]. Defaults to `FALSE`.
#' @param show_messages If `TRUE`, Stan sampler messages are shown.
#' @param seed Optional integer seed passed to CmdStan's sampler for
#'   reproducibility. When `seed` is `NULL`, cmdstanr draws one from R's
#'   random number generator, so `set.seed()` also makes runs reproducible.
#' @details Variable columns in `"simple"` and `"tibble"` output appear in
#'   sorted variable order (`sort(vars)`), regardless of the order in which
#'   variables appear in the polynomial's terms.
#'
#'   When the variety has several disconnected components, HMC chains cannot
#'   cross the low-density regions between them, so the share of draws in
#'   each component reflects how many chains initialized near it rather than
#'   the components' relative mass. For such varieties, consider increasing
#'   `chains`, or use `rejection = TRUE`, which samples components in
#'   proportion to their mass.
#'
#'   With `homo = TRUE`, the Stan programs normalize a polynomial system's
#'   residual through its Jacobian, which must have full rank away from the
#'   system's singular points, so systems with dependent equations are
#'   rejected; drop the redundant equations or use `homo = FALSE`. With
#'   `rejection = TRUE`, draws follow the density [pdvnorm()] defines, and
#'   only systems where it is 0 everywhere off the variety are rejected.
#'
#'   With `noise = "domain"`, each draw is `x = u + e`, where `u` is nearly
#'   uniform on the variety (by arclength, area, and so on, and within the
#'   window `w` when one is given) and `e` is normal with mean 0 and
#'   covariance from `sd` or `Sigma`. The points `u` are draws of the
#'   homogeneous variety normal with scale `pilot_sd`, each moved to the
#'   closest point on the variety; draws whose projection fails, as near
#'   singular points, are replaced by further pilot runs, whose chains
#'   continue the numbering in `.chain`. The window bounds `u`, not `x`.
#'
#'   Uniform points need a variety of finite size: compact, or bounded by
#'   the window `w`. They also need equations without repeated factors. On
#'   a component where `poly` has a factor such as `(x - 1)^2`, the pilot's
#'   residual is about half the distance to the variety, so that component
#'   would get twice its share; `rvnorm()` detects this and stops. The
#'   approximation is of order `pilot_sd^2` where the variety is smooth, and
#'   isolated singular points, such as crossings, do not affect it. As with
#'   codomain noise, chains do not cross between disconnected components, so
#'   `rejection = TRUE` suits such varieties better. `"tibble"` output adds
#'   the points `u` as `<var>_proj` columns; `"stanfit"` output and
#'   `inc_warmup = TRUE` are not available. A given `seed` also seeds the
#'   noise, without changing R's random number stream.
#'
#' @name rvnorm
#' @return Either (1) a data frame whose rows are the individual draws from the
#'   distribution, (2) a [tbl_df-class] object with the draws along with
#'   additional information (with `noise = "domain"`, the points on the
#'   variety as `<var>_proj` columns), or (3) an object of class
#'   `CmdStanMCMC`.
#' @examples
#'
#'
#' library("tidyverse")
#' library("mpoly")
#' options("mc.cores" = parallel::detectCores() - 1)
#'
#' \dontrun{
#'
#' # single polynomial
#' p <- mp("x^2 + y^2 - 1")
#' samps <- rvnorm(1000, p, sd = .05)
#' head(samps)
#' str(samps)
#' plot(samps, asp = 1)
#'
#'
#'
#' # returning a data frame
#' (samps <- rvnorm(1e5, p, sd = .05, output = "tibble", chains = 8, cores = 8))
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
#' f <- function(x, y) pdvnorm(cbind(x, y), poly = p, sd = .05)
#'
#' ggplot(samps, aes(x, y)) +
#'   geom_hdr(xlim = c(-2,2), ylim = c(-2,2)) +
#'   # geom_hdr(xlim = c(-2,2), ylim = c(-2,2), method = method_kde(adjust = c(0.2, 0.2))) +
#'   geom_hdr_lines_fun(fun = f, xlim = c(-2,2), ylim = c(-2,2), color = "cyan", linewidth = .25) +
#'   # geom_point(aes(x, y), size = .25, data = samps) +
#'   geom_variety(poly = p) +
#'   coord_equal()
#'
#' ggplot() +
#'   geom_hdr_fun(fun = f, xlim = c(-2,2), ylim = c(-2,2)) +
#'   # geom_point(aes(x, y), size = .25, data = samps) +
#'   geom_variety(poly = p) +
#'   coord_equal()
#'
#'
#' # in three variables
#' (samps <- rvnorm(20, mp("x^2 + y^2 + z^2 - 1"), sd = .05, w = 2))
#' apply(samps, 1, function(v) sqrt(sum(v^2)))
#'
#'
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
#'
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
#'
#'
#' # using refresh to get more info
#'
#'
#' rvnorm(2000, p, sd = .1, "tibble", verbose = TRUE)
#' rvnorm(2000, p, sd = .1, "tibble", refresh = 500)
#' rvnorm(2000, p, sd = .1, "tibble", refresh = 0) # default
#' rvnorm(2000, p, sd = .1, "tibble", refresh = -1)
#'
#'
#'
#' # many chains in parallel
#'
#'
#' options(mc.cores = parallel::detectCores())
#' p <- mp("x^2 + (4 y)^2 - 1")
#' (samps <- rvnorm(1e4, p, sd = .01, "tibble", verbose = TRUE, chains = 8))
#' ggplot(samps, aes(x, y)) + geom_bin2d(binwidth = .01*c(1,1)) + coord_equal()
#' # decrease sd to get more uniform sampling
#'
#'
#'
#' # windowing for unbounded varieties
#'
#'
#'
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
#'
#'
#' # the importance of normalizing
#'
#' # one of the effects of the normalizing is to stabilize variances, making
#' # them roughly equivalent globally over the variety
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
#'
#'
#' # semi-algebraic sets
#'
#' # inside the semialgebraic set x^2 + y^2 <= 1
#' # this is the same as x^2 + y^2 - 1 <= 0, so that
#' # x^2 + y^2 - 1 + s^2 == 0 for some slack variable s
#' # this is the projection of the sphere into the xy-plane
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
#'
#'
#' # alternative representation
#' # x^2 + y^2 - 1 <= 0 iff s^2 (x^2 + y^2 - 1) + 1 == 0
#' # note that its gradient is more complicated
#' p <- mp("s^2 (x^2 + y^2 - 1) + 1")
#' samps <- rvnorm(1e4, p, sd = .1, "tibble", chains = 8, w = 2, refresh = 1e3)
#' ggplot(samps, aes(x, y)) + geom_bin2d(binwidth = .05*c(1,1)) + coord_equal()
#'
#'
#'
#' # keeping the warmup / the importance of multiple chains
#'
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
#'
#'
#' # ideal-variety correspondence considerations
#'
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
#'
#'
#' # neat examples
#'
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
#'
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
#'
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
#'
#'
#'
#' # noise in the domain: points uniform on an ellipse (red) plus correlated
#' # normal noise
#' p <- mp("0.25 x^2 + y^2 - 1")
#' samps <- rvnorm(
#'   2000, p, noise = "domain", Sigma = matrix(c(.02, .015, .015, .02), 2),
#'   output = "tibble"
#' )
#' ggplot(samps, aes(x, y)) +
#'   geom_point(size = .5) +
#'   geom_point(aes(x_proj, y_proj), size = .5, color = "red") +
#'   coord_equal()
#' }
#' @rdname rvnorm
#' @export
rvnorm <- function(
  n, poly, sd, output = "simple", Sigma = NULL, rejection = FALSE,
  chains = 4L, warmup = max(500, ceiling(n / chains / 2)), inc_warmup = FALSE,
  thin = 1L, adapt_delta = .999, max_treedepth = 20L,
  verbose = FALSE, cores = min(chains, getOption("mc.cores", 1L)),
  homo = TRUE, noise = c("codomain", "domain"), pilot_sd = NULL,
  w, vars, numerator, denominator, refresh = 0L,
  code_only = FALSE, pre_compiled = TRUE, user_compiled = FALSE,
  show_messages = FALSE, seed = NULL, ...
) {

  if (
    !is.numeric(n) || length(n) != 1L || !is.finite(n) ||
      n < 1 || n != as.integer(n)
  ) {
    stop("`n` must be a positive integer.", call. = FALSE)
  }
  n <- as.integer(n)
  output <- match.arg(output, c("simple", "tibble", "stanfit"))
  chains <- rvnorm_validate_positive_integer(chains, "chains")
  cores <- rvnorm_validate_positive_integer(cores, "cores")
  warmup <- rvnorm_validate_nonnegative_integer(warmup, "warmup")
  thin <- rvnorm_validate_positive_integer(thin, "thin")
  if (
    !is.numeric(adapt_delta) || length(adapt_delta) != 1L ||
      !is.finite(adapt_delta) || adapt_delta <= 0 || adapt_delta >= 1
  ) {
    stop("`adapt_delta` must be a finite numeric scalar between 0 and 1.", call. = FALSE)
  }
  max_treedepth <- rvnorm_validate_positive_integer(max_treedepth, "max_treedepth")
  if (!is.logical(inc_warmup) || length(inc_warmup) != 1L || is.na(inc_warmup)) {
    stop("`inc_warmup` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    stop("`verbose` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(homo) || length(homo) != 1L || is.na(homo)) {
    stop("`homo` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(code_only) || length(code_only) != 1L || is.na(code_only)) {
    stop("`code_only` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(pre_compiled) || length(pre_compiled) != 1L || is.na(pre_compiled)) {
    stop("`pre_compiled` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(user_compiled) || length(user_compiled) != 1L || is.na(user_compiled)) {
    stop("`user_compiled` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(show_messages) || length(show_messages) != 1L || is.na(show_messages)) {
    stop("`show_messages` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.numeric(refresh) || length(refresh) != 1L || !is.finite(refresh)) {
    stop("`refresh` must be a finite numeric scalar.", call. = FALSE)
  }
  if (!is.null(seed) && (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed))) {
    stop("`seed` must be NULL or a finite numeric scalar.", call. = FALSE)
  }
  noise <- match.arg(noise)
  if (
    !is.null(pilot_sd) && (
      !is.numeric(pilot_sd) || length(pilot_sd) != 1L ||
        !is.finite(pilot_sd) || pilot_sd <= 0
    )
  ) {
    stop("`pilot_sd` must be NULL or a positive finite numeric scalar.", call. = FALSE)
  }
  if (!missing(numerator) || !missing(denominator)) {
    warning(
      "`numerator` and `denominator` are deprecated and ignored.",
      call. = FALSE
    )
  }
  vars_arg <- if (missing(vars)) NULL else vars

  if (is.character(poly)) {
    poly <- mp(poly)
  }
  # an mpolyList of one polynomial is that polynomial, as in pdvnorm(), so
  # the single-polynomial programs and templates serve it
  if (is.mpolyList(poly) && length(poly) == 1L) {
    poly <- poly[[1]]
  }

  if (noise == "domain") {
    if (!homo) {
      stop(
        "`homo` describes codomain noise; `noise = \"domain\"` takes its ",
        "points on the variety from the homogeneous variety normal.",
        call. = FALSE
      )
    }
    if (inc_warmup) {
      stop(
        "`inc_warmup = TRUE` is not available with `noise = \"domain\"`: ",
        "warmup draws do not follow the pilot distribution, so their ",
        "projections are not uniform on the variety.",
        call. = FALSE
      )
    }
    if (output == "stanfit") {
      stop(
        "`output = \"stanfit\"` is not available with `noise = \"domain\"`; ",
        "\"tibble\" output includes the points on the variety.",
        call. = FALSE
      )
    }
    # sampler settings for the pilot draws; rvnorm_domain() sets their
    # scale, window, and seed
    pilot_args <- c(
      list(
        chains = chains, warmup = warmup, thin = thin,
        adapt_delta = adapt_delta, max_treedepth = max_treedepth,
        verbose = verbose, cores = cores, refresh = refresh,
        code_only = code_only, pre_compiled = pre_compiled,
        user_compiled = user_compiled, show_messages = show_messages
      ),
      if (!is.null(vars_arg)) list(vars = vars_arg),
      list(...)
    )
    return(rvnorm_domain(
      n = n, poly = poly, sd = if (missing(sd)) NULL else sd, Sigma = Sigma,
      pilot_sd = pilot_sd, output = output, rejection = rejection,
      seed = seed, vars = vars_arg, w = if (missing(w)) NULL else w,
      code_only = code_only, pilot_args = pilot_args
    ))
  }

  if (rejection) {
    if (!is.null(Sigma)) {
      stop("`rejection = TRUE` supports `sd`, not `Sigma`.", call. = FALSE)
    }
    if (missing(w)) {
      # rejection_sampler()'s default box; w = 1 clipped the band around
      # unit-scale varieties such as the unit circle
      w <- 1.25
    }
    if (missing(sd) || !(is.numeric(sd) && is.vector(sd) && length(sd) == 1L)) {
      stop("This sd is not supported yet for rejection sampler.")
    }
    sample <- rejection_sampler(
      n = n,
      poly = poly,
      sd = sd,
      vars = rvnorm_resolve_vars(vars_arg, poly),
      w = w,
      output = output,
      homo = homo,
      message = verbose
    )
    # rejection_sampler() returns a matrix; "simple" output is a data frame
    if (output == "simple") sample <- as.data.frame(sample)
    return(sample)
  }

  output_needs_rewriting <- FALSE

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
  n_vars <- length(mpoly::vars(poly))
  if (n_vars == 0L) {
    stop("`poly` must contain at least one variable.", call. = FALSE)
  }
  # variables to sample (output names, fixed before any x/y/z remapping
  # below), which `vars` can extend past those of poly
  sample_vars <- rvnorm_resolve_vars(vars_arg, poly)
  extra_vars <- setdiff(sample_vars, mpoly::vars(poly))
  # the model's names for sample_vars; the templates rename them
  stan_vars <- sample_vars
  w_supplied <- !missing(w)
  if (w_supplied) {
    w_info <- rvnorm_normalize_window(w, sample_vars)
    w <- w_info$value
    w_is_scalar <- w_info$scalar
  } else {
    w_is_scalar <- FALSE
  }
  unbounded_vars <- if (!w_supplied) {
    extra_vars
  } else if (w_is_scalar) {
    character()
  } else {
    setdiff(extra_vars, names(w))
  }
  if (length(unbounded_vars) > 0L) {
    stop(
      "The distribution is flat along variables in `vars` that `poly` does ",
      "not involve; bound ", paste(unbounded_vars, collapse = ", "),
      " with `w`.",
      call. = FALSE
    )
  }
  scale_info <- rvnorm_resolve_scale(
    poly = poly,
    sd = if (missing(sd)) NULL else sd,
    Sigma = Sigma,
    has_sd = !missing(sd),
    homo = homo,
    n_vars = n_vars,
    n_eqs = n_eqs
  )
  sigma_stan <- scale_info$si
  if (user_compiled && !rvnorm_is_scalar(sigma_stan)) {
    stop(
      "`user_compiled = TRUE` currently supports only scalar `sd` or scalar `Sigma`.",
      call. = FALSE
    )
  }

  if (w_supplied && !w_is_scalar && (pre_compiled || user_compiled)) {
    if (isTRUE(verbose)) {
      message("non-scalar `w` requires generated Stan code; falling back to regular rvnorm sampling")
    }
    pre_compiled <- FALSE
    user_compiled <- FALSE
  }
  if (length(extra_vars) > 0L && user_compiled) {
    if (isTRUE(verbose)) {
      message("compile_stan_code() models sample only the variables of `poly`; ignoring `user_compiled`")
    }
    user_compiled <- FALSE
  }

  if (isTRUE(verbose) && refresh == 0L) {
    refresh <- max(ceiling(n / 10), 1L)
  }

  if (
    n_eqs > 1 ||
      length(sample_vars) > 3 ||
      base::max(mpoly::totaldeg(poly)) > 3
  ) {
    # precompiled binaries only cover a small catalog of low-degree templates
    pre_compiled <- FALSE
  }

  if (code_only) {
    stan_code <- if (w_supplied) {
      create_stan_code(poly, sigma_stan, n_eqs, w, homo, sample_vars)
    } else {
      create_stan_code(poly, sigma_stan, n_eqs, homo = homo, vars = sample_vars)
    }
    return(stan_code)
  }

  # the programs solve with J (J \ g, J*J', or J'*J), which fails everywhere
  # when the equations are dependent
  if (homo && n_eqs > 1L && dependent_equations(poly, sort(mpoly::vars(poly)))) {
    stop(
      "With `homo = TRUE`, the Jacobian of `poly` must have full rank away ",
      "from its singular points, but it loses rank everywhere (the equations ",
      "are dependent). Remove the dependent equations or use `homo = FALSE`.",
      call. = FALSE
    )
  }

  if (user_compiled) {
    # look up a previously compiled user template from the internal cache
    model_name <- generate_model_name(
      poly = poly, homo = homo, windowed = w_supplied
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
    stan_data <- if (!w_supplied) {
      c(stan_data, "si" = sigma_stan)
    } else {
      c(stan_data, "w" = w, "si" = sigma_stan)
    }
  } else if (pre_compiled) {
    poly_original <- poly
    output_needs_rewriting_original <- output_needs_rewriting
    # remap variable names to x/y/z because shipped models are keyed that way
    list_for_transformation <- check_and_replace_vars(poly_original, sample_vars)
    if (!is.null(unlist(list_for_transformation$mapping))) {
      poly <- list_for_transformation$polynomial
      var_info <- list_for_transformation$mapping
      output_needs_rewriting <- TRUE
    }

    deg <- mpoly::totaldeg(poly)
    # the template for all sampled variables; the coefficients of terms in
    # variables poly does not involve are 0
    num_of_vars <- length(sample_vars)
    stan_file_name <- paste(
      num_of_vars,
      deg,
      if (homo) "vn" else "hvn",
      sep = "_"
    )
    if (w_supplied) {
      stan_file_name <- paste(stan_file_name, "w", sep = "_")
    }

    model <- tryCatch(
      instantiate::stan_package_model(
        name = stan_file_name,
        package = "vnorm"
      ),
      error = function(e) NULL
    )
    if (!model_has_usable_executable(model)) {
      rvnorm_note_missing_templates()
      poly <- poly_original
      output_needs_rewriting <- output_needs_rewriting_original
      stan_code <- if (w_supplied) {
        create_stan_code(poly, sigma_stan, n_eqs, w, homo, sample_vars)
      } else {
        create_stan_code(poly, sigma_stan, n_eqs, homo = homo, vars = sample_vars)
      }
      stan_file <- write_stan_file(stan_code)
      model <- cmdstan_model(stan_file)
      stan_data <- list("si" = sigma_stan)
    } else {
      stan_vars <- c("x", "y", "z")[seq_len(num_of_vars)]
      stan_data <- make_coefficients_data(poly, num_of_vars, deg)
      stan_data <- if (!w_supplied) {
        c(stan_data, "si" = sigma_stan)
      } else {
        c(stan_data, "w" = w, "si" = sigma_stan)
      }
    }
  } else {
    # fall back to generating and compiling a temporary Stan model
    stan_code <- if (w_supplied) {
      create_stan_code(poly, sigma_stan, n_eqs, w, homo, sample_vars)
    } else {
      create_stan_code(poly, sigma_stan, n_eqs, homo = homo, vars = sample_vars)
    }
    stan_file <- write_stan_file(stan_code)
    model <- cmdstan_model(stan_file)
    stan_data <- list("si" = sigma_stan)
  }

  sample_args <- list(
    data = write_stan_data(stan_data),
    refresh = refresh,
    iter_warmup = warmup,
    # thin = k keeps every k-th iteration, so run k times as many to return n
    iter_sampling = ceiling(n * thin / chains),
    chains = chains,
    parallel_chains = cores,
    adapt_delta = adapt_delta,
    max_treedepth = max_treedepth,
    thin = thin,
    save_warmup = inc_warmup,
    show_messages = show_messages,
    ...
  )
  # seed is injected separately to avoid passing NULL to cmdStan
  if (!is.null(seed)) sample_args$seed <- seed
  # cmdstanr names output files after the minute and a draw from R's RNG,
  # so runs after the same set.seed() shared them, and garbage collecting
  # the first run's fit deleted the second's files before they were read.
  # tempfile() names do not depend on R's RNG
  if (is.null(sample_args$output_basename)) {
    sample_args$output_basename <- basename(tempfile("vnorm-"))
  }
  samps <- do.call(model$sample, sample_args)

  if (output == "simple") {
    df <- as.data.frame(samps$draws(format = "df", inc_warmup = inc_warmup))
    if (!inc_warmup) df <- tail(df, n)
    df <- df[, stan_vars, drop = FALSE]
    row.names(df) <- NULL
    if (output_needs_rewriting) {
      df <- rename_output_df(df, replacement_list = var_info)
    }
    # match the package-wide sorted variable-order convention
    df <- df[, sort(names(df)), drop = FALSE]
    return(df)
  }

  if (output == "tibble") {
    df <- as.data.frame(samps$draws(format = "df", inc_warmup = inc_warmup))
    if (!inc_warmup) df <- tail(df, n)
    row.names(df) <- NULL
    # move the log-density diagnostic column to the end, by name
    if ("lp__" %in% names(df)) {
      df <- df[, c(setdiff(names(df), "lp__"), "lp__"), drop = FALSE]
    }
    if (output_needs_rewriting) {
      df <- rename_output_df(df, replacement_list = var_info)
    }
    # variable columns first in sorted order, then the model's other
    # quantities (g, ndg, ...) in Stan's order, then sampler diagnostics
    diag_cols <- intersect(names(df), c(".chain", ".iteration", ".draw", "lp__"))
    var_cols <- intersect(sample_vars, names(df))
    other_cols <- setdiff(names(df), c(var_cols, diag_cols))
    df <- df[, c(var_cols, other_cols, diag_cols), drop = FALSE]
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

rvnorm_validate_positive_integer <- function(x, label) {
  if (
    !is.numeric(x) || length(x) != 1L || !is.finite(x) ||
      x < 1 || x != as.integer(x)
  ) {
    stop("`", label, "` must be a positive integer.", call. = FALSE)
  }
  as.integer(x)
}

rvnorm_validate_nonnegative_integer <- function(x, label) {
  if (
    !is.numeric(x) || length(x) != 1L || !is.finite(x) ||
      x < 0 || x != as.integer(x)
  ) {
    stop("`", label, "` must be a non-negative integer.", call. = FALSE)
  }
  as.integer(x)
}

rvnorm_is_scalar <- function(x) {
  is.numeric(x) && is.null(dim(x)) && length(x) == 1L
}

rvnorm_check_positive_numeric <- function(x, label) {
  if (!is.numeric(x) || any(!is.finite(x))) {
    stop("`", label, "` must be finite numeric.", call. = FALSE)
  }
  if (any(x <= 0)) {
    stop("`", label, "` must be positive.", call. = FALSE)
  }
  invisible(x)
}

rvnorm_check_positive_definite <- function(x, label) {
  tryCatch(
    chol(x),
    error = function(e) {
      stop("`", label, "` must be positive definite.", call. = FALSE)
    }
  )
  invisible(x)
}

rvnorm_resolve_scale <- function(poly, sd, Sigma, has_sd, homo, n_vars, n_eqs) {
  if (is.null(Sigma) && !has_sd) {
    stop("`sd` must be supplied unless `Sigma` is supplied.", call. = FALSE)
  }

  if (is.mpoly(poly)) {
    if (!is.null(Sigma)) {
      if (!is.numeric(Sigma) || any(!is.finite(Sigma))) {
        stop("`Sigma` must be finite numeric.", call. = FALSE)
      }
      if (!rvnorm_is_scalar(Sigma)) {
        stop(
          "For single-polynomial sampling, `Sigma` must be a single positive variance.",
          call. = FALSE
        )
      }
      rvnorm_check_positive_numeric(Sigma, "Sigma")
      return(list(si = sqrt(as.numeric(Sigma)), target_dim = 1L))
    }

    rvnorm_check_positive_numeric(sd, "sd")
    if (!rvnorm_is_scalar(sd)) {
      stop(
        "For single-polynomial sampling, `sd` must be a single positive numeric.",
        call. = FALSE
      )
    }
    return(list(si = as.numeric(sd), target_dim = 1L))
  }

  target_dim <- if (homo) n_vars else n_eqs
  dim_label <- if (homo) "number of variables" else "number of polynomials"

  if (!is.null(Sigma)) {
    if (!is.numeric(Sigma) || any(!is.finite(Sigma))) {
      stop("`Sigma` must be finite numeric.", call. = FALSE)
    }
    if (rvnorm_is_scalar(Sigma)) {
      rvnorm_check_positive_numeric(Sigma, "Sigma")
      return(list(si = sqrt(as.numeric(Sigma)), target_dim = target_dim))
    }
    if (is.null(dim(Sigma)) && is.vector(Sigma)) {
      rvnorm_check_positive_numeric(Sigma, "Sigma")
      if (length(Sigma) != target_dim) {
        stop(
          "Length of `Sigma` must match the ", dim_label, ".",
          call. = FALSE
        )
      }
      return(list(si = diag(as.numeric(Sigma)), target_dim = target_dim))
    }
    Sigma <- as.matrix(Sigma)
    if (!all(dim(Sigma) == c(target_dim, target_dim))) {
      stop(
        "`Sigma` matrix dimensions must match the ", dim_label, ".",
        call. = FALSE
      )
    }
    rvnorm_check_positive_definite(Sigma, "Sigma")
    return(list(si = Sigma, target_dim = target_dim))
  }

  if (is.matrix(sd)) {
    stop("Use `Sigma` to supply a covariance matrix.", call. = FALSE)
  }
  rvnorm_check_positive_numeric(sd, "sd")
  if (rvnorm_is_scalar(sd)) {
    return(list(si = as.numeric(sd), target_dim = target_dim))
  }
  if (!is.vector(sd)) {
    stop("`sd` must be a scalar or vector.", call. = FALSE)
  }
  if (length(sd) != target_dim) {
    stop(
      "Length of `sd` must match the ", dim_label, ".",
      call. = FALSE
    )
  }
  list(si = diag(as.numeric(sd)^2), target_dim = target_dim)
}

rvnorm_note_missing_templates <- function() {
  # installations can lack the template binaries (R CMD check skips them, as
  # do VNORM_PRECOMPILE_STAN=false and a failed compile), and rvnorm() then
  # compiles every model; say so once per session
  if (isTRUE(.vnorm_state$noted_missing_templates)) return(invisible(FALSE))
  .vnorm_state$noted_missing_templates <- TRUE
  message(
    "vnorm's pre-compiled Stan models are not installed, so rvnorm() ",
    "compiles a model instead, which is slower. Reinstalling vnorm with ",
    "CmdStan available builds them. This message is shown once per session."
  )
  invisible(TRUE)
}

model_has_usable_executable <- function(model) {
  if (is.null(model)) return(FALSE)

  exe_ref <- model$exe_file
  if (is.null(exe_ref)) {
    # test doubles may only provide a sample() method; treat those as usable
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
