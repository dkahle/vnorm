#' Remove Files Created by [compile_stan_code()]
#'
#' Delete Stan source and executable files tracked by the internal compile
#' cache created by [compile_stan_code()].
#'
#' @param path Directory used in [compile_stan_code()].
#'   Deprecated and ignored; compiled files are removed from the internal cache.
#'
#'
#' @export
remove_stan_files <- function(path = getwd()) {
  # remove cached .Stan files and corresponding executables, then clear cache
  if (!missing(path)) {
    warning(
      "`path` is deprecated and ignored; removing files tracked in the internal cache.",
      call. = FALSE
    )
  }
  compiled_stan_info <- get_compiled_stan_info()

  if (nrow(compiled_stan_info) > 0) {
    stan_file_names <- compiled_stan_info$path
    exe_file_names <- sub("\\.stan$", "", stan_file_names)

    file.remove(stan_file_names[file.exists(stan_file_names)])
    file.remove(exe_file_names[file.exists(exe_file_names)])
    clear_compiled_stan_info()
  }
  message("Files deleted.")
}
