# Remove Files Created by [`compile_stan_code()`](https://dkahle.github.io/vnorm/reference/compile_stan_code.md)

Delete Stan source and executable files tracked by the internal compile
cache created by
[`compile_stan_code()`](https://dkahle.github.io/vnorm/reference/compile_stan_code.md).

## Usage

``` r
remove_stan_files(path = getwd())
```

## Arguments

- path:

  Directory used in
  [`compile_stan_code()`](https://dkahle.github.io/vnorm/reference/compile_stan_code.md).
  Defaults to the current working directory.
