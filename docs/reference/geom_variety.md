# One-Dimensional Varieties in Two Dimensions

Plot implicit polynomial varieties with `ggplot2`.

## Usage

``` r
stat_variety(
  mapping = NULL,
  data = NULL,
  geom = GeomVariety,
  position = "identity",
  ...,
  poly = NULL,
  n = 201,
  nx = n,
  ny = n,
  xlim = NULL,
  ylim = NULL,
  shift = 0,
  projection = c("auto", "on", "off"),
  na.rm = FALSE,
  show.legend = NA,
  inherit.aes = TRUE
)

geom_variety(
  mapping = NULL,
  data = NULL,
  stat = "variety",
  position = "identity",
  ...,
  poly,
  vary_colour = FALSE,
  shift = 0,
  projection = c("auto", "on", "off"),
  na.rm = FALSE,
  show.legend = NA,
  inherit.aes = TRUE
)
```

## Arguments

- mapping:

  Aesthetic mappings created with
  [`ggplot2::aes()`](https://ggplot2.tidyverse.org/reference/aes.html).

- data:

  Layer data.

- geom:

  The geometric object used to display data; defaults to GeomVariety.

- position:

  Position adjustment for the layer.

- ...:

  Additional parameters forwarded to
  [`ggplot2::layer()`](https://ggplot2.tidyverse.org/reference/layer.html).

- poly:

  An `mpoly` or `mpolyList` object describing the implicit polynomial(s)
  to plot.

- n:

  Number of grid points used in both x and y directions when `nx` and
  `ny` are not supplied.

- nx, ny:

  Number of grid points in the x and y directions.

- xlim, ylim:

  Length-2 numeric vectors giving plotting limits. If not supplied,
  limits are taken from the plot scales.

- shift:

  A numeric constant added to the evaluated surface before contouring,
  i.e. the plotted level set is `poly + shift = 0`. This is mainly
  useful when the polynomial does not cross zero on the plotting grid
  (for example, `p^2`). If `shift = 0` and all sampled values have one
  sign, `geom_variety()` prints a message suggesting a shift value.

- projection:

  Whether to project contour points back onto the true variety after
  contour extraction. `"off"` returns the raw shifted level-0 contour.
  `"on"` always projects. `"auto"` projects when a shift is used, when
  the polynomial has no strict sign change on the plotting grid, or
  when the raw contour is noticeably off the zero set. For shifted
  repeated-factor cases, `geom_variety()` also prints a caution that the
  recovered contour may still miss branches.

- na.rm:

  If `FALSE`, the default, missing values are removed with a warning. If
  `TRUE`, missing values are silently removed.

- show.legend:

  logical. Should this layer be included in the legends? `NA`, the
  default, includes if any aesthetics are mapped. `FALSE` never
  includes, and `TRUE` always includes. It can also be a named logical
  vector to finely select the aesthetics to display. To include legend
  keys for all levels, even when no data exists, use `TRUE`. If `NA`,
  all levels are shown in legend, but unobserved levels are omitted.

- inherit.aes:

  If `FALSE`, overrides the default aesthetics, rather than combining
  with them. This is most useful for helper functions that define both
  data and aesthetics and shouldn't inherit behaviour from the default
  plot specification, e.g.
  [`annotation_borders()`](https://ggplot2.tidyverse.org/reference/annotation_borders.html).

- stat:

  The statistical transformation to use on the data for this layer. When
  using a `geom_*()` function to construct a layer, the `stat` argument
  can be used to override the default coupling between geoms and stats.
  The `stat` argument accepts the following:

  - A `Stat` ggproto subclass, for example `StatCount`.

  - A string naming the stat. To give the stat as a string, strip the
    function name of the `stat_` prefix. For example, to use
    [`stat_count()`](https://ggplot2.tidyverse.org/reference/geom_bar.html),
    give the stat as `"count"`.

  - For more information and other ways to specify the stat, see the
    [layer
    stat](https://ggplot2.tidyverse.org/reference/layer_stats.html)
    documentation.

- vary_colour:

  Logical. If `TRUE`, map colour to the polynomial label so users can
  control per-polynomial colours with `scale_colour_*()`. Defaults to
  `FALSE`, which keeps a constant line colour and varies only linetype
  across an `mpolyList`.

## Details

`geom_variety()` extracts contour paths from a grid evaluation of the
polynomial and, by default, projects those paths back onto `poly = 0`
when that helps recover the intended zero set.

## Aesthetics

`geom_variety()` understands the following aesthetics. `x` and `y` are
computed by the stat, so users typically do not map them manually:

- x

- y

- alpha

- color

- group

- linetype

- linewidth

- subgroup

## Computed variables

- Polynomial:

  A parseable label for the polynomial, useful for
  `after_stat(Polynomial)` mappings (for example, linetype or colour).

- group:

  Contour path group identifier used internally by the layer.

## See also

[`ggplot2::geom_path()`](https://ggplot2.tidyverse.org/reference/geom_path.html)

## Examples

``` r
library("ggplot2")

# 1) Ellipse
p1 <- mp("x^2 + 4 y^2 - 1")
ggplot() +
  geom_variety(poly = p1, xlim = c(-2, 2), ylim = c(-2, 2)) +
  coord_equal()


# Works with standard ggplot2 styling
ggplot() +
  geom_variety(
    poly = p1, xlim = c(-2, 2), ylim = c(-2, 2),
    colour = "steelblue", linewidth = 0.5
  ) +
  coord_equal() +
  theme_minimal()


# 2) Folium of Descartes (singular variety)
p2 <- mp("x^3 + y^3 - 3 x y")
ggplot() +
  geom_variety(poly = p2, xlim = c(-2, 3), ylim = c(-2, 3)) +
  coord_equal()


# 3) "Heart" curve (classic implicit heart)
p3 <- mp("(x^2 + y^2 - 1)^3 - x^2 y^3")
ggplot() +
  geom_variety(poly = p3, xlim = c(-2, 2), ylim = c(-2, 2)) +
  coord_equal() +
  theme(legend.position = "top")


# 4) A 2-polynomial system (mpolyList): circle and xy = 0.25
p4 <- mp(c("x^2 + y^2 - 1", "x y - 0.25"))
# By default, polynomials differ by linetype (not color).
ggplot() +
  geom_variety(poly = p4, xlim = c(-2, 2), ylim = c(-2, 2)) +
  coord_equal()


# With different colors (optional)
ggplot() +
  geom_variety(poly = p4, xlim = c(-2, 2), ylim = c(-2, 2), vary_colour = TRUE) +
  coord_equal() +
  scale_colour_manual(values = c("steelblue", "firebrick"))


# You can also customize linetypes and legend placement with ggplot2 scales/themes
ggplot() +
  geom_variety(poly = p4, xlim = c(-2, 2), ylim = c(-2, 2), vary_colour = TRUE) +
  coord_equal() +
  scale_colour_manual(values = c("steelblue", "firebrick")) +
  scale_linetype_manual(values = c("solid", "22"), guide = "none") +
  theme(legend.position = "top")
#> Scale for linetype is already present.
#> Adding another scale for linetype, which will replace the existing scale.


## common contouring situations
########################################

# 5) Squared polynomial (same zero set, but no sign change on the grid)
# geom_variety() will suggest a negative shift when no contour is found.
p5 <- mp("x^2 + y^2 - 1")
ggplot() +
  geom_variety(poly = p5^2, xlim = c(-2, 2), ylim = c(-2, 2)) +
  coord_equal()
#> All values are positive on the plotting grid; try shift = -0.00101684.
#> Zero contours were generated


# Use the suggested shift (your printed value may differ slightly).
ggplot() +
  geom_variety(poly = p5^2, xlim = c(-2, 2), ylim = c(-2, 2), shift = -0.001) +
  coord_equal()
#> Using shift = -0.001 to contour a nearby level set. For repeated-factor varieties, the projected result may still miss branches; if the unsquared polynomial is available, prefer plotting that directly.


# Inspect the raw shifted level set versus the default projected recovery.
p6 <- mp("y^2 - x^2")
ggplot() +
  geom_variety(
    poly = p6^2,
    xlim = c(-2, 2), ylim = c(-2, 2),
    shift = -0.004, projection = "off"
  ) +
  coord_equal()


ggplot() +
  geom_variety(poly = p6^2, xlim = c(-2, 2), ylim = c(-2, 2), shift = -0.004) +
  coord_equal()
#> Using shift = -0.004 to contour a nearby level set. For repeated-factor varieties, the projected result may still miss branches; if the unsquared polynomial is available, prefer plotting that directly.


```
