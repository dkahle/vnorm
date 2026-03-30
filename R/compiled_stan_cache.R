.vnorm_state <- new.env(parent = emptyenv())

# retrieve cached model metadata as a data frame
get_compiled_stan_info <- function() {
  info <- .vnorm_state$compiled_stan_info
  if (is.null(info)) {
    return(data.frame(
      name = character(),
      path = character(),
      stringsAsFactors = FALSE
    ))
  }
  info
}

# replace the full cache table
set_compiled_stan_info <- function(info) {
  .vnorm_state$compiled_stan_info <- info
  invisible(info)
}

# append a model entry, deduplicating by row
add_compiled_stan_info <- function(name, path) {
  info <- get_compiled_stan_info()
  new_row <- data.frame(name = name, path = path, stringsAsFactors = FALSE)
  info <- unique(rbind(info, new_row))
  set_compiled_stan_info(info)
}

# drop all cached model metadata
clear_compiled_stan_info <- function() {
  .vnorm_state$compiled_stan_info <- NULL
  invisible(NULL)
}
