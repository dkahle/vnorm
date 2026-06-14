# Compile a Stan Model for User-Defined Polynomials

Compile and cache a Stan model template for a polynomial shape so
repeated calls to
[`rvnorm()`](https://www.kahle.io/vnorm/reference/rvnorm.md) can reuse
the compiled binary.

## Usage

``` r
compile_stan_code(
  poly,
  custom_stan_code = FALSE,
  windowed = FALSE,
  homo = TRUE
)
```

## Arguments

- poly:

  An `mpoly` or `mpolyList` object.

- custom_stan_code:

  If `TRUE`, a custom model is compiled even if the general case of the
  polynomial is already included during package installation. Defaults
  to `FALSE`.

- windowed:

  If `TRUE`, the compiled Stan model includes a window (box constraint)
  parameter `w` in its data block. Defaults to `FALSE`.

- homo:

  If `TRUE`, sampling is done from a homoskedastic variety normal
  distribution. Defaults to `TRUE`.

  The compiled model metadata is stored in an internal package cache
  used by [`rvnorm()`](https://www.kahle.io/vnorm/reference/rvnorm.md)
  with `user_compiled = TRUE`.

## Examples

``` r

if (FALSE) { # \dontrun{
# compile a model that looks like b0 + bx6 x^6 + by6 y^6 for later input
p <- mp("x^6 + y^6 - 1") # template polynomial
samps <- rvnorm(1000, p, sd = .05)
head(samps)
compile_stan_code(p) # allows to change coefficients
p <- mp("x^6 + 8 y^6 - 1")
rvnorm(1e4, p, .05, user_compiled = TRUE)

} # }
```
