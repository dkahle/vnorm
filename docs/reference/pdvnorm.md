# Pseudo-Density for the Variety Normal Distribution

Evaluate the variety normal pseudo-density in either the homoskedastic
or heteroskedastic setting.

## Usage

``` r
pdvnorm(x, poly, sigma, homo = TRUE, log = FALSE)
```

## Arguments

- x:

  A numeric vector of length equal to the number of variables in `poly`,
  or a numeric matrix/data frame with that many columns (one row per
  evaluation point).

- poly:

  An `mpoly` object (single polynomial) or an `mpolyList` object
  (multiple polynomials).

- sigma:

  For the single-polynomial case, a positive scalar standard deviation.
  For the multi-polynomial case, a scalar, vector, or matrix. If
  `homo = TRUE`, `sigma` must conform to the number of variables; if
  `homo = FALSE`, it must conform to the number of polynomials.

- homo:

  Logical; default is `TRUE`. If `TRUE`, compute the homoskedastic
  variety normal pseudo-density. If `FALSE`, compute the heteroskedastic
  pseudo-density.

- log:

  Logical. If `TRUE`, returns the log of the density.

## Value

A numeric scalar or vector containing the pseudo-density evaluated at
`x`.

## Examples

``` r

library("mpoly")
#> 
#> Attaching package: ‘mpoly’
#> The following object is masked from ‘package:ggplot2’:
#> 
#>     vars

## Single polynomial in one variable
p1 <- mp("x")
pdvnorm(c(-1, 0, 1), p1, sigma = 1)
#> [1] 0.2419707 0.3989423 0.2419707
pdvnorm(0, p1, sigma = 2, log = TRUE)
#> [1] -1.612086

## Multivariate (square) system: two polynomials in two variables
p2 <- mp(c("x", "y"))
x2 <- rbind(c(0, 0), c(1, 2), c(-1, 3))

## Different sigma forms
pdvnorm(x2, p2, sigma = 1)
#> [1] 0.15915494 0.05167004 0.09653235
pdvnorm(x2, p2, sigma = c(1, 2))
#> [1] 0.11253954 0.06412310 0.08764588
pdvnorm(x2, p2, sigma = diag(c(1, 4)))
#> [1] 0.07957747 0.06006823 0.07022687

## Multivariate (underdetermined): one polynomial in two variables
p3 <- mp("x + y")
x3 <- rbind(c(1, 1), c(2, -1), c(0, 3))
pdvnorm(x3, p3, sigma = 1)
#> [1] 0.14676266 0.31069656 0.04204821
pdvnorm(as.data.frame(x3), p3, sigma = 1)
#> [1] 0.14676266 0.31069656 0.04204821

## Multivariate (overdetermined): three polynomials in two variables
p4 <- mp(c("x", "y", "x + y"))
x4 <- rbind(c(1, 2), c(0, -1), c(2, 2))
pdvnorm(x4, p4, sigma = diag(2),    homo = TRUE)
#> [1] 0.05854983 0.14241810 0.02689930
pdvnorm(x4, p4, sigma = c(1, 2, 3), homo = FALSE)
#> [1] 0.0012905390 0.0170882873 0.0000896711
```
