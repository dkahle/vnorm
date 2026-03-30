#' @import mpoly
#' @importFrom ggplot2 ggplot aes ggproto Stat GeomPath layer after_stat labs
#' @importFrom cmdstanr write_stan_file cmdstan_model
#' @importFrom instantiate stan_package_model
#' @importFrom glue glue
#' @importFrom tibble as_tibble tibble
#' @importFrom stats coef deriv optimize quantile reorder rnorm runif
#' @importFrom utils modifyList head tail
#' @name vnorm
#' @aliases vnorm vnorm-package
#' @keywords internal
"_PACKAGE"

utils::globalVariables(c(
  "Polynomial",
  "group",
  "coef",
  "deriv",
  "reorder"
))
