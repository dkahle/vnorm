# Compilation
instantiate::stan_package_configure()
instantiate::stan_package_compile()
document()

# Installation
# Run this from the package root.
install.packages(pkgs = ".", type = "source", repos = NULL)
library(vnorm)

# Stan Files Generator

# num_of_vars_vals <- 1:3
# totaldeg_vals <- 1:3
# homo_vals <- c(TRUE, FALSE)
# w_vals <- c(TRUE, FALSE)
#
# param_grid <- expand.grid(
#   num_of_vars = num_of_vars_vals,
#   totaldeg = totaldeg_vals,
#   homo = homo_vals,
#   w = w_vals
# )
#
# results <- apply(param_grid, 1, function(params) {
#   do.call(make_stan_files, as.list(params))
# })



# Compilation code




