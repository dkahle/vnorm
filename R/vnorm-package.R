#' @import mpoly
#' @importFrom ggplot2 aes ggproto Stat GeomPath layer after_stat
#' @importFrom cmdstanr write_stan_file cmdstan_model
#' @importFrom instantiate stan_package_model
#' @importFrom glue glue
#' @importFrom tibble as_tibble tibble
#' @importFrom stats deriv optimize quantile reorder runif
#' @importFrom utils modifyList tail
#' @name vnorm
#' @aliases vnorm vnorm-package
#' @keywords internal
"_PACKAGE"

utils::globalVariables(c(
  "polynomial",
  "group"
))
