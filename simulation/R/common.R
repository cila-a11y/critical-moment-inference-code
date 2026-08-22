options(warn = 1)

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L || identical(x, "")) y else x
}

stopf <- function(fmt, ...) {
  stop(sprintf(fmt, ...), call. = FALSE)
}

assert_true <- function(condition, message) {
  if (!isTRUE(condition)) {
    stop(message, call. = FALSE)
  }
  invisible(TRUE)
}

stop_numerical <- function(message, call = NULL) {
  condition <- structure(
    list(message = as.character(message), call = call),
    class = c("cmo_expected_numerical_failure", "error", "condition")
  )
  stop(condition)
}

is_expected_numerical_failure <- function(condition) {
  inherits(condition, "cmo_expected_numerical_failure")
}

cmo_root <- function() {
  root <- Sys.getenv("CMO_ROOT", unset = "")
  if (!nzchar(root)) {
    root <- normalizePath(".", winslash = "/", mustWork = TRUE)
  }
  normalizePath(root, winslash = "/", mustWork = TRUE)
}

source_cmo <- function(relative_path, local = .GlobalEnv) {
  source(file.path(cmo_root(), relative_path), local = local, chdir = FALSE)
}

read_config <- function(path) {
  full_path <- if (grepl("^/", path)) path else file.path(cmo_root(), path)
  assert_true(file.exists(full_path), sprintf("Configuration file not found: %s", full_path))
  cfg <- dget(full_path)
  assert_true(is.list(cfg), "The configuration must evaluate to a list.")
  cfg$config_path <- normalizePath(full_path, winslash = "/", mustWork = TRUE)
  cfg
}

parse_named_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  out <- list()
  for (arg in args) {
    pieces <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1L]]
    if (length(pieces) != 2L || !nzchar(pieces[1L])) {
      stopf("Arguments must have the form --name=value; received: %s", arg)
    }
    out[[pieces[1L]]] <- pieces[2L]
  }
  out
}

atomic_publish <- function(temporary, path, overwrite) {
  if (isTRUE(overwrite)) {
    if (!file.rename(temporary, path)) {
      stopf("Could not atomically publish %s", path)
    }
  } else {
    if (!file.link(temporary, path)) {
      stopf("Refusing to overwrite existing immutable artifact %s", path)
    }
    unlink(temporary)
  }
  invisible(path)
}

atomic_save_rds <- function(object, path, compress = "xz", overwrite = TRUE) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(pattern = paste0(".", basename(path), "."), tmpdir = dirname(path))
  on.exit(unlink(temporary), add = TRUE)
  saveRDS(object, temporary, compress = compress)
  atomic_publish(temporary, path, overwrite = overwrite)
}

atomic_write_csv <- function(object, path, overwrite = TRUE) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(pattern = paste0(".", basename(path), "."), tmpdir = dirname(path))
  on.exit(unlink(temporary), add = TRUE)
  write.csv(object, temporary, row.names = FALSE, na = "")
  atomic_publish(temporary, path, overwrite = overwrite)
}

atomic_write_lines <- function(text, path, overwrite = TRUE) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(pattern = paste0(".", basename(path), "."), tmpdir = dirname(path))
  on.exit(unlink(temporary), add = TRUE)
  writeLines(as.character(text), con = temporary)
  atomic_publish(temporary, path, overwrite = overwrite)
}

timestamp_utc <- function() {
  format(Sys.time(), tz = "UTC", format = "%Y-%m-%dT%H:%M:%SZ")
}

log_message <- function(...) {
  text <- paste(..., collapse = "")
  cat(sprintf("[%s] %s\n", timestamp_utc(), text))
  flush.console()
  invisible(text)
}

safe_numeric <- function(x, name) {
  value <- suppressWarnings(as.numeric(x))
  if (length(value) != 1L || !is.finite(value)) {
    stopf("%s must be one finite number; received '%s'.", name, x)
  }
  value
}

safe_integer <- function(x, name, lower = 1L) {
  numeric_value <- suppressWarnings(as.numeric(x))
  if (length(numeric_value) != 1L || !is.finite(numeric_value) ||
      numeric_value != floor(numeric_value) || numeric_value < lower ||
      numeric_value > .Machine$integer.max) {
    stopf("%s must be an integer not smaller than %d; received '%s'.", name, lower, x)
  }
  as.integer(numeric_value)
}

valid_identifier <- function(value) {
  length(value) == 1L && !is.na(value) &&
    grepl("^[A-Za-z0-9][A-Za-z0-9_.-]*$", value)
}

normalize_tokens <- function(value, name, allow_empty = FALSE) {
  if (is.null(value)) {
    value <- character(0L)
  }
  value <- as.character(unlist(value, use.names = FALSE))
  value <- value[nzchar(value)]
  if (!allow_empty) {
    assert_true(length(value) > 0L, sprintf("%s must not be empty.", name))
  }
  assert_true(!anyNA(value), sprintf("%s contains NA.", name))
  assert_true(!anyDuplicated(value), sprintf("%s contains duplicates.", name))
  value
}

collapse_tokens <- function(value) {
  paste(normalize_tokens(value, "token vector", allow_empty = TRUE), collapse = ";")
}

split_tokens <- function(value) {
  if (length(value) == 0L || is.na(value) || !nzchar(value)) {
    return(character(0L))
  }
  strsplit(value, ";", fixed = TRUE)[[1L]]
}

normalize_order_intervals <- function(config) {
  intervals <- config$order_intervals
  if (is.null(intervals)) {
    assert_true(
      !is.null(config$p_min) && !is.null(config$p_max),
      "The configuration needs order_intervals or both p_min and p_max."
    )
    intervals <- list(K = c(config$p_min, config$p_max))
  }

  if (is.data.frame(intervals)) {
    required <- c("interval_id", "p_min", "p_max")
    assert_true(all(required %in% names(intervals)),
                "The order-interval table needs interval_id, p_min, and p_max.")
    out <- lapply(seq_len(nrow(intervals)), function(index) {
      c(intervals$p_min[index], intervals$p_max[index])
    })
    names(out) <- as.character(intervals$interval_id)
  } else {
    assert_true(is.list(intervals) && length(intervals) > 0L,
                "order_intervals must be a nonempty named list or data frame.")
    out <- vector("list", length(intervals))
    interval_names <- names(intervals)
    for (index in seq_along(intervals)) {
      entry <- intervals[[index]]
      if (is.list(entry) && !is.null(entry$interval_id)) {
        interval_names[index] <- as.character(entry$interval_id)
        lower <- entry$p_min %||% entry$lower
        upper <- entry$p_max %||% entry$upper
        out[[index]] <- c(lower, upper)
      } else {
        out[[index]] <- as.numeric(entry)
      }
    }
    if (is.null(interval_names) || any(!nzchar(interval_names))) {
      stop("Every order interval must have a nonempty identifier.", call. = FALSE)
    }
    names(out) <- interval_names
  }

  assert_true(!anyDuplicated(names(out)), "Order-interval identifiers must be unique.")
  for (interval_id in names(out)) {
    bounds <- as.numeric(out[[interval_id]])
    assert_true(valid_identifier(interval_id),
                sprintf("Unsafe interval identifier: %s", interval_id))
    assert_true(length(bounds) == 2L && all(is.finite(bounds)) &&
                  bounds[1L] > 0 && bounds[1L] < bounds[2L],
                sprintf("Invalid bounds for order interval %s.", interval_id))
    out[[interval_id]] <- unname(bounds)
  }
  out
}

merge_named_lists <- function(...) {
  inputs <- list(...)
  output <- list()
  for (input in inputs) {
    if (is.null(input)) next
    assert_true(is.list(input), "Only lists can be merged.")
    for (name in names(input)) {
      output[[name]] <- input[[name]]
    }
  }
  output
}

set_single_thread_math <- function() {
  variables <- c(
    "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
    "BLIS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS",
    "RCPP_PARALLEL_NUM_THREADS"
  )
  for (variable in variables) {
    do.call(Sys.setenv, setNames(list("1"), variable))
  }
  invisible(TRUE)
}

capture_session <- function() {
  list(
    timestamp_utc = timestamp_utc(),
    host = Sys.info()[["nodename"]],
    pid = Sys.getpid(),
    r_version = R.version.string,
    platform = R.version$platform,
    slurm_job_id = Sys.getenv("SLURM_JOB_ID", unset = NA_character_),
    slurm_array_job_id = Sys.getenv("SLURM_ARRAY_JOB_ID", unset = NA_character_),
    slurm_array_task_id = Sys.getenv("SLURM_ARRAY_TASK_ID", unset = NA_character_),
    slurm_cpus_per_task = Sys.getenv("SLURM_CPUS_PER_TASK", unset = NA_character_),
    session_info = utils::capture.output(sessionInfo())
  )
}

project_fingerprint <- function(config_path) {
  root <- cmo_root()
  code_files <- c(
    list.files(file.path(root, "R"), pattern = "\\.R$", full.names = TRUE),
    list.files(file.path(root, "bin"), pattern = "\\.sh$", full.names = TRUE),
    list.files(file.path(root, "slurm"), pattern = "\\.(sh|sbatch)$", full.names = TRUE),
    normalizePath(config_path, winslash = "/", mustWork = TRUE)
  )
  code_files <- sort(unique(code_files[file.exists(code_files)]))
  hashes <- unname(tools::md5sum(code_files))
  names(hashes) <- sub(paste0("^", root, "/?"), "", code_files)
  list(files = names(hashes), md5 = hashes)
}

fingerprints_identical <- function(first, second) {
  is.list(first) && is.list(second) &&
    identical(first$files, second$files) && identical(first$md5, second$md5)
}

design_artifact_paths <- function(design_directory) {
  file.path(
    design_directory,
    c(
      "units.csv", "tasks.csv", "queue.rds", "cells.rds", "truth.csv",
      "seed_registry.csv", "seeds.rds", "config_snapshot.rds",
      "fingerprint.rds"
    )
  )
}

compute_design_fingerprint <- function(design_directory) {
  paths <- design_artifact_paths(design_directory)
  assert_true(all(file.exists(paths)), "The frozen design is incomplete.")
  hashes <- unname(tools::md5sum(paths))
  names(hashes) <- basename(paths)
  list(files = names(hashes), md5 = hashes)
}
