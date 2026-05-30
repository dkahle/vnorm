# Solve a Zero-Dimensional Variety

Estimates the solution of a zero-dimensional variety by sampling with
[`rvnorm()`](https://www.kahle.io/vnorm/reference/rvnorm.md) and
averaging the draws. For well-behaved varieties, this method works
reliably. When multiple isolated components exist, the sampler may favor
some over others, and the posterior mean may fall between true
solutions.

## Usage

``` r
variety_solve(
  polylist,
  sd = 0.01,
  n = 1e+05,
  sig_digit = 3,
  vars = mpoly::vars(polylist),
  inc_warmup = FALSE,
  show_message = FALSE,
  stanfit = FALSE,
  ...
)
```

## Arguments

- polylist:

  An `mpolyList` or `mpoly` object containing the polynomials.

- sd:

  Numeric scalar standard deviation used by
  [`rvnorm()`](https://www.kahle.io/vnorm/reference/rvnorm.md)
  (default: `0.01`).

- n:

  Integer, number of draws used for averaging (default: `1e5`).

- sig_digit:

  Integer, significant digits used to round the solution (default: `3`).

- vars:

  Character vector of variable names to average; defaults to
  `mpoly::vars(polylist)`.

- inc_warmup:

  Logical; include warmup draws when extracting (default: `FALSE`).

- show_message:

  Logical; show Stan sampling messages (default: `FALSE`).

- stanfit:

  Logical; if `TRUE`, return both the Stan fit and the solution
  (default: `FALSE`).

- ...:

  Additional arguments passed to
  [`rvnorm()`](https://www.kahle.io/vnorm/reference/rvnorm.md).

## Value

If `stanfit = FALSE`, a named numeric vector of rounded posterior means.
If `stanfit = TRUE`, a list with `stanfit` (Stan fit) and `results` (the
vector).

## Examples

``` r
library(mpoly)
if (FALSE) { # \dontrun{
polylist <- mp(c("x^2 - y", "x^2 + y"))

# Posterior-mean solution only
variety_solve(polylist, n = 2e4, sd = 0.01, sig_digit = 3)

# Return Stan fit and solution
variety_solve(polylist, n = 1e4, stanfit = TRUE)

} # }
```
