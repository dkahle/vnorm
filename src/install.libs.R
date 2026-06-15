libs <- file.path(R_PACKAGE_DIR, "libs", R_ARCH)
dir.create(libs, recursive = TRUE, showWarnings = FALSE)
for (file in c("symbols.rds", Sys.glob(paste0("*", SHLIB_EXT)))) {
  if (file.exists(file)) {
    file.copy(file, file.path(libs, file))
  }
}
inst_stan <- file.path("..", "inst", "stan")
if (dir.exists(inst_stan)) {
  warning(
    "Stan models in inst/stan/ are deprecated in {instantiate} ",
    ">= 0.0.4.9001 (2024-01-03). Please put them in src/stan/ instead."
  )
  if (file.exists("stan")) {
    warning("src/stan/ already exists. Not copying models from inst/stan/.")
  } else {
    message("Copying inst/stan/ to src/stan/.")
    fs::dir_copy(path = inst_stan, new_path = "stan")
  }
}
bin <- file.path(R_PACKAGE_DIR, "bin")
if (!file.exists(bin)) {
  dir.create(bin, recursive = TRUE, showWarnings = FALSE)
}
bin_stan <- file.path(bin, "stan")
dir.create(bin_stan, recursive = TRUE, showWarnings = FALSE)
src_stan <- normalizePath("stan", winslash = "/", mustWork = TRUE)
bin_stan <- normalizePath(bin_stan, winslash = "/", mustWork = TRUE)
running_r_cmd_check_parent <- nzchar(Sys.getenv("_R_CHECK_PACKAGE_NAME_"))
callr::r(
  func = function(src_stan, bin_stan, running_r_cmd_check_parent) {
    sanitize_component <- function(x) {
      gsub("[^[:alnum:]_.-]", "_", x)
    }
    running_r_cmd_check <- isTRUE(running_r_cmd_check_parent) ||
      nzchar(Sys.getenv("_R_CHECK_PACKAGE_NAME_"))
    precompile_opt <- tolower(Sys.getenv("VNORM_PRECOMPILE_STAN", "auto"))
    if (!precompile_opt %in% c("auto", "true", "false")) precompile_opt <- "auto"
    enable_stan_binaries <- switch(
      precompile_opt,
      "true" = TRUE,
      "false" = FALSE,
      "auto" = !running_r_cmd_check
    )

    model_binary_exists <- function(dir, model_stem) {
      candidates <- c(
        file.path(dir, model_stem),
        file.path(dir, paste0(model_stem, ".exe"))
      )
      any(file.exists(candidates))
    }

    find_cache_root <- function() {
      user_cache <- tryCatch(
        tools::R_user_dir("vnorm", which = "cache"),
        error = function(e) NA_character_
      )
      candidates <- unique(c(user_cache, file.path(tempdir(), "vnorm-cache")))
      candidates <- candidates[!is.na(candidates)]
      for (candidate in candidates) {
        ok <- tryCatch({
          dir.create(candidate, recursive = TRUE, showWarnings = FALSE)
          if (!dir.exists(candidate)) return(FALSE)
          probe <- tempfile("write-probe-", tmpdir = candidate)
          on.exit(unlink(probe, force = TRUE), add = TRUE)
          saveRDS(TRUE, probe)
          file.exists(probe)
        }, error = function(e) FALSE)
        if (isTRUE(ok)) return(candidate)
      }
      stop("Could not create a cache directory for Stan model binaries.")
    }

    cmdstan_version <- tryCatch(
      as.character(cmdstanr::cmdstan_version()),
      error = function(e) "unknown"
    )
    cmdstan_path <- tryCatch(
      as.character(cmdstanr::cmdstan_path()),
      error = function(e) "unknown"
    )
    # Keep cache entries scoped to platform + CmdStan version + CmdStan path.
    cache_stan <- file.path(
      find_cache_root(),
      "stan",
      sanitize_component(R.version$platform),
      paste0("cmdstan-", sanitize_component(cmdstan_version)),
      paste0("path-", sanitize_component(cmdstan_path))
    )
    dir.create(cache_stan, recursive = TRUE, showWarnings = FALSE)
    cache_probe <- tempfile("write-probe-", tmpdir = cache_stan)
    tryCatch(
      {
        saveRDS(TRUE, cache_probe)
        unlink(cache_probe, force = TRUE)
      },
      error = function(e) {
        stop(
          "Could not write to the Stan model cache directory: ",
          cache_stan,
          call. = FALSE
        )
      }
    )

    source_models <- sort(
      list.files(
        path = src_stan,
        pattern = "[.]stan$",
        full.names = TRUE
      )
    )
    if (length(source_models) == 0) {
      stop("No .stan model files found in source directory: ", src_stan)
    }

    model_files <- basename(source_models)
    model_stems <- tools::file_path_sans_ext(model_files)

    model_hashes <- as.character(tools::md5sum(source_models))
    names(model_hashes) <- model_files

    cache_index_path <- file.path(cache_stan, "compile-index.rds")
    if (file.exists(cache_index_path)) {
      cache_index <- readRDS(cache_index_path)
    } else {
      cache_index <- data.frame(model = character(), hash = character(), stringsAsFactors = FALSE)
    }
    if (!is.data.frame(cache_index) || !all(c("model", "hash") %in% names(cache_index))) {
      cache_index <- data.frame(model = character(), hash = character(), stringsAsFactors = FALSE)
    }
    cache_hash_lookup <- setNames(cache_index$hash, cache_index$model)

    file.copy(source_models, file.path(cache_stan, model_files), overwrite = TRUE)

    needs_compile <- vapply(seq_along(model_files), function(i) {
      model_file <- model_files[[i]]
      model_stem <- model_stems[[i]]
      cached_hash <- unname(cache_hash_lookup[model_file])
      current_hash <- unname(model_hashes[model_file])
      hash_changed <- length(cached_hash) == 0L ||
        is.na(cached_hash) ||
        !identical(cached_hash, current_hash)
      # Recompile only when source changed or cached binary is missing.
      hash_changed ||
        !model_binary_exists(cache_stan, model_stem)
    }, logical(1))

    models_to_compile <- file.path(cache_stan, model_files[needs_compile])
    if (enable_stan_binaries && length(models_to_compile) > 0) {
      tryCatch(
        instantiate::stan_package_compile(models = models_to_compile),
        error = function(e) {
          warning(
            "Skipping Stan model precompilation during installation: ",
            conditionMessage(e),
            call. = FALSE
          )
        }
      )
    } else if (!enable_stan_binaries) {
      message(
        "Skipping Stan binary precompilation/copy during installation ",
        "(set VNORM_PRECOMPILE_STAN=true to override)."
      )
    }

    has_binary <- vapply(model_stems, function(model_stem) {
      model_binary_exists(cache_stan, model_stem)
    }, logical(1))
    cache_index <- data.frame(
      model = model_files[has_binary],
      hash = unname(model_hashes[model_files[has_binary]]),
      stringsAsFactors = FALSE
    )
    saveRDS(cache_index, cache_index_path)

    existing_bin_files <- list.files(
      path = bin_stan,
      all.files = TRUE,
      full.names = TRUE,
      no.. = TRUE
    )
    if (length(existing_bin_files) > 0) {
      unlink(existing_bin_files, recursive = TRUE, force = TRUE)
    }

    file.copy(file.path(cache_stan, model_files), bin_stan, overwrite = TRUE)
    if (enable_stan_binaries) {
      for (model_stem in model_stems) {
        for (suffix in c("", ".exe")) {
          model_binary <- file.path(cache_stan, paste0(model_stem, suffix))
          if (file.exists(model_binary)) {
            file.copy(model_binary, bin_stan, overwrite = TRUE)
          }
        }
      }
    }
  },
  args = list(
    src_stan = src_stan,
    bin_stan = bin_stan,
    running_r_cmd_check_parent = running_r_cmd_check_parent
  )
)
