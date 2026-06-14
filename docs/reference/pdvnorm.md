# Pseudo-Density for the Variety Normal Distribution

Evaluate the variety normal pseudo-density in either the homoskedastic
or heteroskedastic setting.

## Usage

``` r
pdvnorm(x, poly, sd, homo = TRUE, log = FALSE, Sigma = NULL, ...)
```

## Arguments

- x:

  A numeric vector of length equal to the number of variables in `poly`,
  or a numeric matrix/data frame with that many columns (one row per
  evaluation point).

- poly:

  An `mpoly` object (single polynomial) or an `mpolyList` object
  (multiple polynomials).

- sd:

  Scale parameter for the normal kernel. For polynomial systems, a
  scalar is recycled across the relevant dimension, while a vector
  supplies component-wise standard deviations.

- homo:

  Logical; default is `TRUE`. If `TRUE`, compute the homoskedastic
  variety normal pseudo-density. If `FALSE`, compute the heteroskedastic
  pseudo-density.

- log:

  Logical. If `TRUE`, returns the log of the density.

- Sigma:

  Full covariance matrix, scalar covariance, or a diagonal vector of
  covariance terms. If supplied, `Sigma` replaces `sd`.

- ...:

  Deprecated. A named `sigma` argument is accepted for backward
  compatibility.

## Value

A numeric scalar or vector containing the pseudo-density evaluated at
`x`.

## Examples

``` r

# m = 1 polynomial in n = 1 variable, 0d variety
p <- mp("x")

pdvnorm(c(-1, 0, 1), poly = p, sd = 1)
#> [1] 0.2419707 0.3989423 0.2419707
pdvnorm(c(-1, 0, 1), poly = p, sd = 1, log = TRUE)
#> [1] -1.4189385 -0.9189385 -1.4189385



# m = 2 polynomials in n = 2 variables (square system), 0d variety
p <- mp(c("x", "y"))

X <- rbind(c(-.25, 0), c(1, 2), c(-1, 3))
pdvnorm(X, poly = p, sd = 1)
#> [1] 0.154258260 0.013064233 0.001072378



# m = 1 polynomial in n = 2 variables, 1d variety
p <- mp("x^2 + y^2 - 1")

X <- rbind(c(-.25, 0), c(1, 2), c(-1, 3))
pdvnorm(X, poly = p, sd = 1)
#> [1] 0.06878628 0.26741901 0.14493955



## Different dispersion forms
p <- mp(c("x", "y"))
X <- rbind(c(0, 0), c(1, 2), c(-1, 3))
pdvnorm(X, p, sd = 1)
#> [1] 0.159154943 0.013064233 0.001072378
pdvnorm(X, p, sd = c(1, 2))
#> [1] 0.07957747 0.02927492 0.01566973
pdvnorm(X, p, Sigma = diag(c(1, 4)))
#> [1] 0.07957747 0.02927492 0.01566973



## Multivariate (underdetermined): one polynomial in two variables
p <- mp("x^2 + y^2 - 1")
X <- rbind(c(1, 1), c(2, -1), c(0, 3))
pdvnorm(X, p, sd = 1)
#> [1] 0.3747716 0.2674190 0.1640101
pdvnorm(as.data.frame(X), p, sd = 1)
#> [1] 0.3747716 0.2674190 0.1640101
pdvnorm(c(1, 1), p, Sigma = diag(c(1, 4)))
#> [1] 0.2460836



## Multivariate (overdetermined): three polynomials in two variables
p <- mp(c("x", "y", "x + y"))
X <- rbind(c(1, 2), c(0, -1), c(2, 2))
pdvnorm(X, p, Sigma = diag(2), homo = TRUE)
#> [1] 0.013064233 0.096532353 0.002915024
pdvnorm(X, p, sd = c(1, 2, 3), homo = FALSE)
#> [1] 0.002361224 0.008834148 0.000357111



if (FALSE) { # \dontrun{

library("ggplot2")
library("ggfunction")



# m = 1 polynomial in n = 1 variable, 0d variety
p <- mp("x")

ggplot() +
  geom_pdf(fun = pdvnorm, args = list(poly = p, sd = 1), xlim = c(-5, 5))



# m = 2 polynomials in n = 2 variables (square system), 0d variety
p <- mp(c("x", "y"))

ggplot() +
  geom_pdf_2d(
    fun = pdvnorm, args = list(poly = p, sd = 1),
    xlim = c(-5, 5), ylim = c(-5, 5)
  ) +
  coord_equal()



# m = 1 polynomial in n = 2 variables, 1d variety
p <- mp("x^2 + y^2 - 1")

ggplot() +
  geom_pdf_2d(
    fun = pdvnorm, args = list(poly = p, sd = .1),
    xlim = c(-1.5, 1.5), ylim = c(-1.5, 1.5)
  ) +
  coord_equal()



# m = 1 polynomial in n = 2 variables, 1d variety
p <- mp("-1 x^2 (x + 1) + y^2", varorder = c("x", "y"))

ggplot() +
  geom_pdf_2d(
    fun = pdvnorm, args = list(poly = p, sd = .1),
    xlim = c(-1.5, 1.5), ylim = c(-1.5, 1.5)
  ) +
  coord_equal()

ggplot() +
  geom_pdf_2d(
    fun = pdvnorm, args = list(poly = p, sd = .1, homo = FALSE),
    xlim = c(-1.5, 1.5), ylim = c(-1.5, 1.5)
  ) +
  coord_equal()



# m = 1 polynomial in n = 2 variables, 1d varieties
si <- .025
w <- 1.5

p1 <- ggplot() +
  geom_pdf_2d(
    fun = pdvnorm, args = list(poly = mp("x^2 + (4 y)^2 - 1"), sd = si),
    xlim = c(-w, w), ylim = c(-w, w)
  ) +
  coord_equal(xlim = c(-w, w), ylim = c(-w, w))

p2 <- ggplot() +
  geom_pdf_2d(
    fun = pdvnorm, args = list(poly = mp("-1 (x - y) (x + y)"), sd = si),
    xlim = c(-w, w), ylim = c(-w, w)
  ) +
  coord_equal(xlim = c(-w, w), ylim = c(-w, w))

p3 <- ggplot() +
  geom_pdf_2d(
    fun = pdvnorm, args = list(poly = mp("(x^2 + y^2)^3 - 4 x^2 y^2"), sd = si),
    xlim = c(-w, w), ylim = c(-w, w)
  ) +
  coord_equal(xlim = c(-w, w), ylim = c(-w, w))

p4 <- ggplot() +
  geom_pdf_2d(
    fun = pdvnorm, args = list(poly = mp("(x^2 + y^2 - 1)^3 - x^2 y^3"), sd = si),
    xlim = c(-w, w), ylim = c(-w, w)
  ) +
  coord_equal(xlim = c(-w, w), ylim = c(-w, w))

library("patchwork")
p1 + p2 + p3 + p4 +
  plot_layout(nrow = 1, widths = 1, guides = "collect")



# m = 2 polynomials in n = 2 variables, 0d varieties
# different representations and "covariances"

Sigma <- .2^2
w <- 1.5

p1 <- ggplot() +
  geom_pdf_2d(
    fun = pdvnorm, args = list(poly = mp(c("x", "y")), Sigma = Sigma),
    xlim = c(-w, w), ylim = c(-w, w)
  ) +
  coord_equal(xlim = c(-w, w), ylim = c(-w, w)) +
  labs(title = latex2exp::TeX(r"(V(x, y))"))

p2 <- ggplot() +
  geom_pdf_2d(
    fun = pdvnorm, args = list(poly = mp(c("-1 (x - y)", "(x + y)")), Sigma = Sigma),
    xlim = c(-w, w), ylim = c(-w, w)
  ) +
  coord_equal(xlim = c(-w, w), ylim = c(-w, w)) +
  labs(title = latex2exp::TeX(r"(V(y - x, y + x))"))

p3 <- ggplot() +
  geom_pdf_2d(
    fun = pdvnorm, args = list(poly = mp(c("(-1 (x - y))^2", "(x + y)^2")), Sigma = Sigma),
    xlim = c(-w, w), ylim = c(-w, w)
  ) +
  coord_equal(xlim = c(-w, w), ylim = c(-w, w)) +
  labs(title = latex2exp::TeX(r"(V((y - x)^2, (y + x)^2))"))

p4 <- ggplot() +
  geom_pdf_2d(
    fun = pdvnorm, args = list(poly = mp(c("x y^3 - x^3 y", "x^2 + y^2 - 1")), Sigma = Sigma),
    xlim = c(-w, w), ylim = c(-w, w)
  ) +
  coord_equal(xlim = c(-w, w), ylim = c(-w, w)) +
  labs(title = latex2exp::TeX(r"(V(x y^3 - x^3 y, x^2 + y^2 - 1))"))

p1 + p2 + p3 + p4 +
  plot_layout(nrow = 1, widths = 1, guides = "collect")



S <- (.20 * diag(c(1, 1)))
R <- matrix(c(1, .90, .90, 1), nrow = 2)
Sigma <- S %*% R %*% S
w <- 1.5

p1 <- ggplot() +
  geom_pdf_2d(
    fun = pdvnorm, args = list(poly = mp(c("x", "y")), Sigma = Sigma),
    xlim = c(-w, w), ylim = c(-w, w)
  ) +
  coord_equal(xlim = c(-w, w), ylim = c(-w, w)) +
  labs(title = latex2exp::TeX(r"(V(x, y))"))

p2 <- ggplot() +
  geom_pdf_2d(
    fun = pdvnorm, args = list(poly = mp(c("-1 (x - y)", "(x + y)")), Sigma = Sigma),
    xlim = c(-w, w), ylim = c(-w, w)
  ) +
  coord_equal(xlim = c(-w, w), ylim = c(-w, w)) +
  labs(title = latex2exp::TeX(r"(V(y - x, y + x))"))

p3 <- ggplot() +
  geom_pdf_2d(
    fun = pdvnorm, args = list(poly = mp(c("(-1 (x - y))^2", "(x + y)^2")), Sigma = Sigma),
    xlim = c(-w, w), ylim = c(-w, w)
  ) +
  coord_equal(xlim = c(-w, w), ylim = c(-w, w)) +
  labs(title = latex2exp::TeX(r"(V((y - x)^2, (y + x)^2))"))

p4 <- ggplot() +
  geom_pdf_2d(
    fun = pdvnorm, args = list(poly = mp(c("x y^3 - x^3 y", "x^2 + y^2 - 1")), Sigma = Sigma),
    xlim = c(-w, w), ylim = c(-w, w)
  ) +
  coord_equal(xlim = c(-w, w), ylim = c(-w, w)) +
  labs(title = latex2exp::TeX(r"(V(x y^3 - x^3 y, x^2 + y^2 - 1))"))

p1 + p2 + p3 + p4 +
  plot_layout(nrow = 1, widths = 1, guides = "collect")



# geom_hdr_fun accepts functions f(x, y) with vectors x and y,
# but pdvnorm accepts the packed f(matrix), so a wrapper is needed.
library("ggdensity")

p <- mp("x^2 + y^2 - 1")
f <- function(x, y, ...) pdvnorm(cbind(x, y), ...)

f(1, 2, poly = p, sd = .05)
f(1:2, 2:3, poly = p, sd = .05)

ggplot() +
  geom_hdr_fun(
    fun = f, xlim = c(-1.25, 1.25), ylim = c(-1.25, 1.25),
    args = list(poly = p, sd = .10)
  ) +
  coord_equal()

S <- (.10 * diag(c(1, 1)))
R <- matrix(c(1, .95, .95, 1), nrow = 2)

ggplot() +
  geom_hdr_fun(
    fun = f, xlim = c(-1.25, 1.25), ylim = c(-1.25, 1.25),
    args = list(poly = p, Sigma = S %*% R %*% S), n = 250
  ) +
  coord_equal()
} # }
```
