# Rejection Sampler for the Variety Normal Distribution

Perform rejection sampling to generate draws from a variety normal
distribution.

## Usage

``` r
rejection_sampler(
  n,
  poly,
  sd = 0.01,
  vars = sort(mpoly::vars(poly)),
  w = 1.25,
  output = "simple",
  dist = c("norm", "unif"),
  homo = TRUE,
  correct_p_coefficients = FALSE,
  correct_dp_coefficients = FALSE,
  message = FALSE
)
```

## Arguments

- n:

  The number of accepted draws to return.

- poly:

  An `mpoly` object or an `mpolyList` object.

- sd:

  The "standard deviation" component of the normal kernel.

- vars:

  A character vector of the indeterminates in the distribution.

- w:

  Proposal box constraints. If a single number, a box window `(-w, w)`
  is applied to all variables. If length 2, the same interval is used
  for all variables. A named list can be used to specify bounds for each
  variable.

- output:

  Either `"simple"` or `"tibble"` output format.

- dist:

  Either `"norm"` (normal) or `"unif"` (uniform).

- homo:

  If `TRUE`, sampling is done from a homoskedastic variety normal
  distribution.

- correct_p_coefficients:

  If `TRUE`, normalize polynomial coefficients.

- correct_dp_coefficients:

  If `TRUE`, normalize derivative coefficients.

- message:

  If `TRUE`, print progress messages showing remaining samples.

## Value

A matrix or tibble containing the accepted samples.

## Examples

``` r
if (FALSE) { # \dontrun{
library("mpoly")

# single polynomial (circle)
p1 <- mp("x^2 + y^2 - 1")
set.seed(1)
rejection_sampler(100, p1, sd = 0.05, w = 1.5)



# uniform band proposal around the variety, returning a tibble
rejection_sampler(
  100, p1, sd = 0.05, w = c(-1.5, 1.5),
  dist = "unif", output = "tibble"
)



# two-polynomial system (upper/lower acceptance geometry differs by `homo`)
p2 <- mp(c("x^2 + y^2 - 1", "y"))
rejection_sampler(50, p2, sd = 0.05, w = 1.5, homo = TRUE)
rejection_sampler(50, p2, sd = c(0.05, 0.05), w = 1.5, homo = FALSE)
} # }
```
