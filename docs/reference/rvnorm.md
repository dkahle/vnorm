# The Variety Normal Distribution

Unnormalized density evaluation and random generation for the variety
normal distribution with mean defined by `poly` and scale `sd`. See
Details for caveats.

## Usage

``` r
rvnorm(
  n,
  poly,
  sd,
  output = "simple",
  Sigma = NULL,
  rejection = FALSE,
  chains = 4L,
  warmup = max(500, floor(n/2)),
  inc_warmup = FALSE,
  thin = 1L,
  verbose = FALSE,
  cores = min(chains, getOption("mc.cores", 1L)),
  homo = TRUE,
  w,
  vars,
  numerator,
  denominator,
  refresh = 0L,
  code_only = FALSE,
  pre_compiled = TRUE,
  user_compiled = FALSE,
  show_messages = FALSE,
  seed = NULL,
  ...
)
```

## Arguments

- n:

  Number of draws desired from each chain after warmup.

- poly:

  An `mpoly` object.

- sd:

  Scale parameter for the normal kernel. If `Sigma` is supplied, `sd` is
  replaced by `Sigma`. For single-polynomial sampling, `sd` must be a
  positive scalar. For polynomial-list sampling, `sd` may be a positive
  scalar or a positive vector matching the Stan target dimension.

- output:

  One of `"simple"`, `"tibble"`, or `"stanfit"`.

- Sigma:

  Covariance parameter. For single-polynomial sampling, only a positive
  scalar variance is supported. For polynomial-list sampling, `Sigma`
  may be a positive scalar variance, a positive diagonal vector, or a
  positive-definite covariance matrix matching the Stan target
  dimension.

- rejection:

  If `TRUE`, rejection sampling is used.

- chains:

  The number of chains to run for the random number generation, see
  `stan()`.

- warmup:

  Number of warmup iterations in `stan()`.

- inc_warmup:

  If `TRUE`, the MCMC warmup steps are included in the output.

- thin:

  `stan()` `thin` parameter.

- verbose:

  If `TRUE`, print additional progress messages.

- cores:

  The number of CPU cores to distribute the chains across, see `stan()`.

- homo:

  If `TRUE`, sampling is from a homoskedastic variety normal
  distribution.

- w:

  A named list of box constraints for vectors to be passed to Stan, see
  examples. If a single number, a box window `(-w, w)` is applied to all
  variables.

- vars:

  Character vector of polynomial indeterminates.

- numerator, denominator:

  Character scalars containing printed numerator and denominator forms.

- refresh:

  The `refresh` argument of `stan()`, which governs how much information
  is provided to the user while sampling.

- code_only:

  If `TRUE`, only formulate and return Stan code.

- pre_compiled:

  Whether to use precompiled Stan models. Available for polynomials with
  three indeterminates and three degrees. Defaults to `TRUE`.

- user_compiled:

  If `TRUE`, use a user-compiled Stan program produced by
  [`compile_stan_code()`](https://www.kahle.io/vnorm/reference/compile_stan_code.md).
  Defaults to `FALSE`.

- show_messages:

  If `TRUE`, Stan sampler messages are shown.

- seed:

  Optional integer seed passed to CmdStan's sampler for reproducibility.
  Note that R's [`set.seed()`](https://rdrr.io/r/base/Random.html) does
  not control Stan's RNG.

- ...:

  Additional parameters passed to `stan()`.

## Value

Either (1) matrix whose rows are the individual draws from the
distribution, (2) a
[tbl_df](https://rdrr.io/pkg/tibble/man/tbl_df-class.html) object with
the draws along with additional information, or (3) an object of class
`stanfit`.

## Details

If the variety you are interested in is connected, this strategy should
work well out of the box. If it is not, you will likely need multiple
chains, and sampling may become biased toward one or more components
while down-sampling others. Question: what is the relative likelihood of
each component, or an equal unit of length, on different components? How
does this generalize to more varieties of varying dimensions?

## Examples

``` r


library("tidyverse")
#> ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
#> ✔ dplyr     1.2.1     ✔ readr     2.2.0
#> ✔ forcats   1.0.1     ✔ stringr   1.6.0
#> ✔ lubridate 1.9.5     ✔ tibble    3.3.1
#> ✔ purrr     1.2.2     ✔ tidyr     1.3.2
#> ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
#> ✖ dplyr::filter() masks stats::filter()
#> ✖ dplyr::lag()    masks stats::lag()
#> ✖ dplyr::vars()   masks mpoly::vars(), ggplot2::vars()
#> ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors
options("mc.cores" = parallel::detectCores() - 1)

if (FALSE) { # \dontrun{
# basic usage


# single polynomial
p <- mp("x^2 + y^2 - 1")
samps <- rvnorm(1000, p, sd = .05)
head(samps)
str(samps)
plot(samps, asp = 1)



# returning a data frame
(samps <- rvnorm(1e5, p, sd = .05, output = "tibble", chains = 8, cores = 8))

ggplot(samps, aes(x, y)) +
  geom_point(size = .5) +
  coord_equal()

ggplot(samps, aes(x, y)) +
  geom_point(size = .5) +
  geom_variety(poly = p) +
  coord_equal()

ggplot(samps, aes(x, y)) +
  geom_point(aes(color = factor(.chain)), size = .5) +
  geom_variety(poly = p, linewidth = 1) +
  coord_equal()

ggplot(samps, aes(x, y)) +
  stat_density2d(
    aes(fill = after_stat(density)),
    geom = "raster", contour = FALSE
   ) +
  geom_variety(poly = p) +
  coord_equal()

library("ggdensity")
f <- function(x, y) pdvnorm(cbind(x, y), poly = p, sd = .05)

ggplot(samps, aes(x, y)) +
  geom_hdr(xlim = c(-2,2), ylim = c(-2,2), method = method_kde(adjust = c(0.2, 0.2))) +
  geom_hdr_lines_fun(fun = f, xlim = c(-2,2), ylim = c(-2,2), color = "cyan", linewidth = .25) +
  # geom_point(aes(x, y), size = .25, data = samps) +
  geom_variety(poly = p) +
  coord_equal()

ggplot() +
  geom_hdr_fun(fun = f, xlim = c(-2,2), ylim = c(-2,2)) +
  # geom_point(aes(x, y), size = .25, data = samps) +
  geom_variety(poly = p) +
  coord_equal()


# in three variables
(samps <- rvnorm(20, mp("x^2 + y^2 + z^2 - 1"), sd = .05, w = 2))
apply(samps, 1, function(v) sqrt(sum(v^2)))



# more than one polynomial, # vars > # eqns, underdetermined system
p <- mp(c("x^2 + y^2 + z^2 - 1", "z"))
(samps <- rvnorm(500, p, sd = .1, output = "tibble"))
ggplot(samps, aes(x, y)) + geom_point(size = .5) + coord_equal()

ggplot(samps, aes(x, y, color = `g[1]`)) + geom_point() +
  scale_color_gradient2(mid = "gray80") + coord_equal()

ggplot(samps, aes(x, y, color = `g[2]`)) + geom_point() +
  scale_color_gradient2(mid = "gray80") + coord_equal()

ggplot(samps, aes(x, z, color = `g[1]`)) + geom_point() +
  scale_color_gradient2(mid = "gray80") + coord_equal()



# more than one polynomial, # vars < # eqns, overdetermined system
p <- mp(c("3 x", "3 y", "2 x + 2 y", "3 (x^2 + y)", "3 (x^2 - y)"))
(samps <- rvnorm(500, p, sd = .1, output = "tibble"))

samps |>
  select(x, y, starts_with("g")) |>
  pivot_longer(
    starts_with("g"),
    names_to = "equation",
    values_to = "value"
  ) |>
  ggplot(aes(x, y, color = value)) + geom_point() +
    scale_color_gradient2(mid = "gray80") + coord_equal() +
    facet_wrap(~ equation)



# using refresh to get more info


rvnorm(2000, p, sd = .1, "tibble", verbose = TRUE)
rvnorm(2000, p, sd = .1, "tibble", refresh = 500)
rvnorm(2000, p, sd = .1, "tibble", refresh = 0) # default
rvnorm(2000, p, sd = .1, "tibble", refresh = -1)



# many chains in parallel


options(mc.cores = parallel::detectCores())
p <- mp("x^2 + (4 y)^2 - 1")
(samps <- rvnorm(1e4, p, sd = .01, "tibble", verbose = TRUE, chains = 8))
ggplot(samps, aes(x, y)) + geom_bin2d(binwidth = .01*c(1,1)) + coord_equal()
# decrease sd to get more uniform sampling



# windowing for unbounded varieties



# windowing is needed for unbounded varieties
# in the following, look at the parameters block

p <- mp("x y - 1") # unbounded variety, 1 poly
p <- mp(c("x y - 1", "y - x")) # 2 polys

rvnorm(1e3, p, sd = .01, "tibble", code_only = TRUE)
rvnorm(1e3, p, sd = .01, "tibble",  w = 1.15)

window <- list("x" = c(-1.5, 1.5))
rvnorm(1e3, p, sd = .01, "tibble",  w = window)



# the importance of normalizing

# one of the effects of the normalizing is to stabilize variances, making
# them roughly equivalent globally over the variety

# lemniscate of bernoulli
p <- mp("(x^2 + y^2)^2 - 2 (x^2 - y^2)")

# normalized, good
(samps <- rvnorm(2000, p, .05, "tibble"))
ggplot(samps, aes(x, y)) + geom_point(size = .5) + coord_equal()
ggplot(samps, aes(x, y)) + geom_bin2d(binwidth = .05*c(1,1)) + coord_equal()

# unnormalized, bad
(samps <- rvnorm(2000, p, .05, "tibble", homo = FALSE))
ggplot(samps, aes(x, y)) + geom_point(size = .5) + coord_equal()
ggplot(samps, aes(x, y)) + geom_bin2d(binwidth = .05*c(1,1)) + coord_equal()



# semi-algebraic sets

# inside the semialgebraic set x^2 + y^2 <= 1
# this is the same as x^2 + y^2 - 1 <= 0, so that
# x^2 + y^2 - 1 + s^2 == 0 for some slack variable s
# this is the projection of the sphere into the xy-plane

p <- mp("1 - (x^2 + y^2) - s^2")
samps <- rvnorm(1e4, p, sd = .1, "tibble", chains = 8, refresh = 1e3)
ggplot(samps, aes(x, y)) + geom_bin2d(binwidth = .05*c(1,1)) + coord_equal()

ggplot(sample_n(samps, 2e3), aes(x, y, color = s)) +
  geom_point(size = .5) +
  scale_color_gradient2() +
  coord_equal()



# alternative representation
# x^2 + y^2 - 1 <= 0 iff s^2 (x^2 + y^2 - 1) + 1 == 0
# note that its gradient is more complicated
p <- mp("s^2 (x^2 + y^2 - 1) + 1")
samps <- rvnorm(1e4, p, sd = .1, "tibble", chains = 8, w = 2, refresh = 1e3)
ggplot(samps, aes(x, y)) + geom_bin2d(binwidth = .05*c(1,1)) + coord_equal()



# keeping the warmup / the importance of multiple chains


p <- mp("((x + 1.5)^2 + y^2 - 1) ((x - 1.5)^2 + y^2 - 1)")
ggplot() +
  geom_variety(poly = p, xlim = c(-3,3)) +
  theme(legend.position = "top") +
  coord_equal()

# notice the migration of chains initialized away from the distribution
# (it helps to make the graphic large on your screen)
samps <- rvnorm(500, p, sd = .05, "tibble", chains = 8, inc_warmup = TRUE)



# ideal-variety correspondence considerations


p <- mp("x^2 + y^2 - 1")

samps_1 <- rvnorm(250, p^1, sd = .1, output = "tibble", chains = 8)
samps_2 <- rvnorm(250, p^2, sd = .1, output = "tibble", chains = 8)
samps_3 <- rvnorm(250, p^3, sd = .1, output = "tibble", chains = 8)
samps_4 <- rvnorm(250, p^4, sd = .1, output = "tibble", chains = 8)
samps <- bind_rows(mget(apropos("samps_[1-4]")))
samps$power <- rep(seq_along(apropos("samps_[1-4]")), each = 250)

ggplot(samps, aes(x, y, color = g < 0)) +
  geom_point(size = .5) +
  coord_equal(xlim = c(-3,3), ylim = c(-3,3)) +
  facet_wrap(~ power)



# neat examples

# an implicit Lissajous region, view in separate window large

# x = cos(m t + p)
# y = sin(n t + q)
(p <- lissajous(3, 2,  -pi/2, 0))
(p <- lissajous(4, 3,  -pi/2, 0))
(p <- lissajous(5, 4,  -pi/2, 0))
(p <- lissajous(3, 3,  0, 0))
(p <- lissajous(5, 5,  0, 0))
(p <- lissajous(7, 7,  0, 0))
ggplot() +
  geom_variety(poly = p, xlim = c(-2, 2), ylim = c(-2, 2), n = 201,
  show.legend = FALSE) +
  coord_equal()

p <- plug(p, "x", mp(".5 x"))
p <- plug(p, "y", mp(".5 y"))



# algebraic set
samps <- rvnorm(5e3, p, sd = .01, "tibble", chains = 8, cores = 8)
ggplot(samps, aes(x, y, color = factor(.chain))) +
  geom_point(size = .5) +
  coord_equal()

ggplot(samps, aes(x, y)) + geom_hdr() + coord_equal()

ggplot(samps, aes(x, y, color = factor(.chain))) +
  geom_point(size = .5) +
  coord_equal() +
  facet_wrap(~ factor(.chain))



# semi-algebraic set
samps_normd <- rvnorm(1e4, p + mp("s^2"), sd = .01, "tibble", chains = 8,
  cores = 8, homo = TRUE
)
samps_unormd <- rvnorm(1e4, p + mp("s^2"), sd = .01, "tibble", chains = 8,
  cores = 8, homo = FALSE
)

bind_rows(
  samps_normd  |> mutate(normd = TRUE),
  samps_unormd |> mutate(normd = FALSE)
) |>
  ggplot(aes(x, y)) +
    geom_point(size = .5) +
    # geom_bin2d(binwidth = .05*c(1,1)) +
    facet_grid(normd ~ .chain) +
    coord_equal()

ggplot(samps_normd, aes(x, y)) +
  geom_bin2d(binwidth = .05*c(1,1)) + coord_equal()
} # }
```
