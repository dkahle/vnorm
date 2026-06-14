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

## Single polynomial in one variable
p1 <- mp("x")
pdvnorm(c(-1, 0, 1), p1, sd = 1)
#> [1] 0.2419707 0.3989423 0.2419707
pdvnorm(0, p1, sd = 2, log = TRUE)
#> [1] -1.612086

## Single polynomial in two variables
p2 <- mp("x^2 + y^2 - 1")
x2 <- rbind(c(1, 0), c(0, 1), c(1, 1))
pdvnorm(x2, p2, sd = 0.1)
#> [1] 3.989422804 3.989422804 0.007701398
pdvnorm(as.data.frame(x2), p2, sd = 0.1)
#> [1] 3.989422804 3.989422804 0.007701398
pdvnorm(c(1, 1), p2, Sigma = diag(c(1, 4)))
#> [1] 0.2460836

## Polynomial systems
p3 <- mp(c("x", "y"))
x3 <- rbind(c(0, 0), c(1, 2), c(-1, 3))
pdvnorm(x3, p3, sd = 1)
#> [1] 0.159154943 0.013064233 0.001072378
pdvnorm(x3, p3, sd = c(1, 2))
#> [1] 0.07957747 0.02927492 0.01566973
pdvnorm(x3, p3, Sigma = diag(c(1, 4)))
#> [1] 0.07957747 0.02927492 0.01566973

## Overdetermined systems use the variable dimension when homo = TRUE
## and the equation dimension when homo = FALSE.
p4 <- mp(c("x", "y", "x + y"))
x4 <- rbind(c(1, 2), c(0, -1), c(2, 2))
pdvnorm(x4, p4, Sigma = diag(2), homo = TRUE)
#> [1] 0.013064233 0.096532353 0.002915024
pdvnorm(x4, p4, sd = c(1, 2, 3), homo = FALSE)
#> [1] 0.002361224 0.008834148 0.000357111
```
