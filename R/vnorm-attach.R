.onLoad <- function(libname, pkgname) {
  # silence cmdStanr version-check startup chatter for package users
  options(cmdstanr_no_ver_check = TRUE)
}
