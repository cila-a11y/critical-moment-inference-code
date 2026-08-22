#!/usr/bin/env Rscript

source(file.path(Sys.getenv("CMO_ROOT", unset = "."), "R", "common.R"))

arguments <- parse_named_args()
config_argument <- arguments$config %||% Sys.getenv("CMO_CONFIG", unset = "")
assert_true(nzchar(config_argument), "Supply --config=... or set CMO_CONFIG.")
config <- read_config(config_argument)
design_directory <- file.path(cmo_root(), "design", config$study_id)
required_paths <- c(
  design_artifact_paths(design_directory),
  file.path(design_directory, "design_fingerprint.rds")
)
assert_true(all(file.exists(required_paths)), "The frozen design is incomplete.")

queue <- readRDS(file.path(design_directory, "queue.rds"))
units <- queue$units
tasks <- queue$tasks
cells_object <- readRDS(file.path(design_directory, "cells.rds"))
seed_registry <- read.csv(
  file.path(design_directory, "seed_registry.csv"),
  stringsAsFactors = FALSE
)
seeds <- readRDS(file.path(design_directory, "seeds.rds"))
stored_config <- readRDS(file.path(design_directory, "config_snapshot.rds"))
stored_fingerprint <- readRDS(file.path(design_directory, "fingerprint.rds"))
stored_design_fingerprint <- readRDS(file.path(design_directory, "design_fingerprint.rds"))
current_fingerprint <- project_fingerprint(config$config_path)
assert_true(identical(stored_config, config),
            "The supplied configuration differs from the frozen snapshot.")
assert_true(
  fingerprints_identical(stored_fingerprint, current_fingerprint),
  "The current code/configuration differs from the frozen study fingerprint."
)
assert_true(
  fingerprints_identical(
    stored_design_fingerprint,
    compute_design_fingerprint(design_directory)
  ),
  "A frozen design artifact differs from its recorded fingerprint."
)

required_unit_columns <- c(
  "unit_id", "seed_id", "bootstrap_seed_id", "task_id", "queue_position",
  "cell_index", "cell_id", "cell_class", "scenario_id", "n_x", "n_y", "rep_id",
  "pairing_key", "multiplier_pairing_key", "variant_id", "interval_ids",
  "moment_types"
)
required_task_columns <- c(
  "task_id", "first_queue_position", "last_queue_position", "first_unit",
  "last_unit", "number_of_units", "number_of_cells", "number_of_seed_streams"
)
assert_true(all(required_unit_columns %in% names(units)),
            "units.csv has an incompatible schema.")
assert_true(all(required_task_columns %in% names(tasks)),
            "tasks.csv has an incompatible schema.")
assert_true(nrow(units) > 0L && nrow(tasks) > 0L, "The design is empty.")
assert_true(identical(as.integer(units$unit_id), seq_len(nrow(units))),
            "unit_id must be contiguous and ordered.")
assert_true(identical(as.integer(units$queue_position), seq_len(nrow(units))),
            "queue_position must be contiguous and ordered.")
assert_true(identical(as.integer(tasks$task_id), seq_len(nrow(tasks))),
            "task_id must be contiguous and ordered.")
assert_true(!anyDuplicated(units$unit_id), "The design has duplicate unit IDs.")
assert_true(all(units$task_id %in% tasks$task_id),
            "A unit refers to an unknown task.")

task_schema_valid <- all(vapply(seq_len(nrow(tasks)), function(index) {
  task <- tasks[index, , drop = FALSE]
  task_units <- units[units$task_id == task$task_id, , drop = FALSE]
  nrow(task_units) == task$number_of_units &&
    min(task_units$queue_position) == task$first_queue_position &&
    max(task_units$queue_position) == task$last_queue_position &&
    min(task_units$unit_id) == task$first_unit &&
    max(task_units$unit_id) == task$last_unit &&
    length(unique(task_units$cell_id)) == task$number_of_cells &&
    length(unique(task_units$seed_id)) +
      length(unique(task_units$bootstrap_seed_id)) == task$number_of_seed_streams
}, logical(1L)))
assert_true(task_schema_valid, "tasks.csv does not match the global unit queue.")

assert_true(
  all(c("stream_role", "seed_id", "stream_key", "rep_id") %in% names(seed_registry)) &&
    setequal(unique(seed_registry$stream_role), c("sample", "bootstrap")),
  "The seed registry has an incompatible role schema."
)
sample_registry <- seed_registry[seed_registry$stream_role == "sample", , drop = FALSE]
bootstrap_registry <- seed_registry[seed_registry$stream_role == "bootstrap", , drop = FALSE]
assert_true(
  identical(as.integer(sample_registry$seed_id), seq_len(nrow(sample_registry))) &&
    identical(as.integer(bootstrap_registry$seed_id), seq_len(nrow(bootstrap_registry))),
  "The sample and bootstrap seed registries must each be contiguous."
)
assert_true(
  is.list(seeds) && is.matrix(seeds$sample) && is.matrix(seeds$bootstrap) &&
    ncol(seeds$sample) == nrow(sample_registry) &&
    ncol(seeds$bootstrap) == nrow(bootstrap_registry) &&
    identical(colnames(seeds$sample), as.character(sample_registry$seed_id)) &&
    identical(colnames(seeds$bootstrap), as.character(bootstrap_registry$seed_id)),
  "seeds.rds does not match the role-specific seed registries."
)
assert_true(all(units$seed_id %in% sample_registry$seed_id) &&
              all(units$bootstrap_seed_id %in% bootstrap_registry$seed_id),
            "A unit refers to an unknown sample or bootstrap stream.")
sample_map <- unique(units[, c("pairing_key", "rep_id", "seed_id")])
names(sample_map)[names(sample_map) == "pairing_key"] <- "stream_key"
bootstrap_map <- unique(
  units[, c("multiplier_pairing_key", "rep_id", "bootstrap_seed_id")]
)
names(bootstrap_map) <- c("stream_key", "rep_id", "seed_id")
assert_true(!anyDuplicated(paste(sample_map$stream_key, sample_map$rep_id)) &&
              !anyDuplicated(paste(bootstrap_map$stream_key, bootstrap_map$rep_id)),
            "A pairing key and replication map to multiple RNG streams.")
sample_merge <- merge(
  sample_map,
  sample_registry[, c("seed_id", "stream_key", "rep_id")],
  by = c("seed_id", "stream_key", "rep_id"),
  all = TRUE
)
bootstrap_merge <- merge(
  bootstrap_map,
  bootstrap_registry[, c("seed_id", "stream_key", "rep_id")],
  by = c("seed_id", "stream_key", "rep_id"),
  all = TRUE
)
assert_true(nrow(sample_merge) == nrow(sample_registry) &&
              nrow(bootstrap_merge) == nrow(bootstrap_registry),
            "The frozen role-specific seed mapping is inconsistent.")

cell_ids <- vapply(cells_object$cells, `[[`, character(1L), "cell_id")
assert_true(!anyDuplicated(cell_ids) && all(units$cell_id %in% cell_ids),
            "The unit queue and cells.rds disagree.")

raw_directory <- file.path(cmo_root(), "results", "raw", config$study_id)
shard_paths <- if (dir.exists(raw_directory)) {
  sort(list.files(
    raw_directory,
    pattern = "^wave_[0-9]{5}\\.rds$",
    recursive = TRUE,
    full.names = TRUE
  ))
} else {
  character(0L)
}

objects <- vector("list", length(shard_paths))
unreadable_paths <- character(0L)
for (index in seq_along(shard_paths)) {
  objects[[index]] <- tryCatch(readRDS(shard_paths[index]), error = function(error) NULL)
  if (is.null(objects[[index]])) {
    unreadable_paths <- c(unreadable_paths, shard_paths[index])
  }
}
readable <- !vapply(objects, is.null, logical(1L))
objects <- objects[readable]
readable_paths <- shard_paths[readable]

compatible <- if (length(objects) == 0L) logical(0L) else vapply(objects, function(object) {
  is.list(object) &&
    identical(object$shard_schema_version, "global-queue-superwave-v2") &&
    identical(object$schema_version, config$schema_version) &&
    identical(object$study_id, config$study_id) &&
    fingerprints_identical(object$fingerprint, stored_fingerprint) &&
    fingerprints_identical(object$design_fingerprint, stored_design_fingerprint)
}, logical(1L))
incompatible_count <- sum(!compatible)
objects <- objects[compatible]
compatible_paths <- readable_paths[compatible]

expected_result_keys <- function(unit) {
  moment_types <- split_tokens(unit$moment_types)
  interval_ids <- split_tokens(unit$interval_ids)
  grid <- expand.grid(
    moment_type = moment_types,
    interval_id = interval_ids,
    stringsAsFactors = FALSE
  )
  paste(grid$moment_type, grid$interval_id, sep = "::")
}

metadata_matches <- function(rows, unit) {
  if (nrow(rows) == 0L) return(TRUE)
  comparisons <- list(
    rows$unit_id == unit$unit_id,
    rows$seed_id == unit$seed_id,
    rows$bootstrap_seed_id == unit$bootstrap_seed_id,
    rows$task_id == unit$task_id,
    as.character(rows$cell_id) == as.character(unit$cell_id),
    as.character(rows$cell_class) == as.character(unit$cell_class),
    as.character(rows$scenario_id) == as.character(unit$scenario_id),
    rows$rep_id == unit$rep_id,
    as.character(rows$pairing_key) == as.character(unit$pairing_key),
    as.character(rows$multiplier_pairing_key) ==
      as.character(unit$multiplier_pairing_key),
    as.character(rows$variant_id) == as.character(unit$variant_id)
  )
  all(vapply(comparisons, function(value) all(!is.na(value) & value), logical(1L)))
}

semantic_valid <- if (length(objects) == 0L) logical(0L) else
  vapply(seq_along(objects), function(index) {
    object <- objects[[index]]
    path <- compatible_paths[index]
    task_id <- as.integer(object$task_id)
    if (length(task_id) != 1L || is.na(task_id) || !task_id %in% tasks$task_id) {
      return(FALSE)
    }
    expected_task_directory <- sprintf("task_%05d", task_id)
    if (!identical(basename(dirname(path)), expected_task_directory)) return(FALSE)
    file_wave_id <- suppressWarnings(as.integer(sub(
      "^wave_([0-9]{5})\\.rds$", "\\1", basename(path)
    )))
    if (is.na(file_wave_id) || !identical(as.integer(object$wave_id), file_wave_id)) {
      return(FALSE)
    }
    if (!is.data.frame(object$results) || !is.data.frame(object$errors)) return(FALSE)
    shard_ids <- as.integer(object$unit_ids)
    if (length(shard_ids) == 0L || anyNA(shard_ids) || anyDuplicated(shard_ids) ||
        !all(shard_ids %in% units$unit_id[units$task_id == task_id])) return(FALSE)

    result_ids <- if (nrow(object$results) > 0L) {
      unique(as.integer(object$results$unit_id))
    } else integer(0L)
    error_ids <- if (nrow(object$errors) > 0L) {
      as.integer(object$errors$unit_id)
    } else integer(0L)
    if (anyDuplicated(error_ids) || length(intersect(result_ids, error_ids)) > 0L ||
        !setequal(c(result_ids, error_ids), shard_ids)) return(FALSE)

    required_metadata <- c(
      "unit_id", "seed_id", "bootstrap_seed_id", "task_id", "cell_id",
      "cell_class", "scenario_id", "rep_id", "pairing_key", "multiplier_pairing_key",
      "variant_id"
    )
    if (nrow(object$results) > 0L &&
        (!all(c(required_metadata, "moment_type") %in% names(object$results)))) {
      return(FALSE)
    }
    if (nrow(object$errors) > 0L &&
        !all(c(required_metadata, "error_class", "error_message") %in% names(object$errors))) {
      return(FALSE)
    }

    result_row_indices <- if (nrow(object$results) > 0L) {
      split(seq_len(nrow(object$results)), as.character(object$results$unit_id))
    } else list()
    error_row_indices <- if (nrow(object$errors) > 0L) {
      split(seq_len(nrow(object$errors)), as.character(object$errors$unit_id))
    } else list()

    for (unit_id in shard_ids) {
      unit <- units[unit_id, , drop = FALSE]
      if (unit$unit_id != unit_id) return(FALSE)
      result_indices <- result_row_indices[[as.character(unit_id)]] %||% integer(0L)
      error_indices <- error_row_indices[[as.character(unit_id)]] %||% integer(0L)
      result_rows <- object$results[result_indices, , drop = FALSE]
      error_rows <- object$errors[error_indices, , drop = FALSE]
      if (!metadata_matches(result_rows, unit) || !metadata_matches(error_rows, unit)) {
        return(FALSE)
      }
      if (nrow(error_rows) == 1L) next
      if (nrow(result_rows) == 0L) return(FALSE)

      interval_ids <- split_tokens(unit$interval_ids)
      if ("interval_id" %in% names(result_rows)) {
        observed_keys <- paste(
          as.character(result_rows$moment_type),
          as.character(result_rows$interval_id),
          sep = "::"
        )
      } else {
        if (length(interval_ids) != 1L) return(FALSE)
        observed_keys <- paste(
          as.character(result_rows$moment_type), interval_ids, sep = "::"
        )
      }
      expected_keys <- expected_result_keys(unit)
      if (anyDuplicated(observed_keys) || !setequal(observed_keys, expected_keys)) {
        return(FALSE)
      }
    }
    TRUE
  }, logical(1L))
malformed_count <- sum(!semantic_valid)
objects <- objects[semantic_valid]
valid_paths <- compatible_paths[semantic_valid]

if (length(objects) > 0L) {
  shard_task_wave_keys <- vapply(objects, function(object) {
    paste(object$task_id, object$wave_id, sep = "::")
  }, character(1L))
  duplicated_task_wave <- unique(shard_task_wave_keys[duplicated(shard_task_wave_keys)])
} else {
  duplicated_task_wave <- character(0L)
}

observed_ids <- if (length(objects) == 0L) integer(0L) else
  as.integer(unlist(lapply(objects, `[[`, "unit_ids")))
duplicate_ids <- sort(unique(observed_ids[duplicated(observed_ids)]))
unexpected_ids <- sort(setdiff(unique(observed_ids), units$unit_id))
missing_ids <- sort(setdiff(units$unit_id, unique(observed_ids)))
missing_task_ids <- sort(unique(units$task_id[units$unit_id %in% missing_ids]))

error_tables <- lapply(objects, `[[`, "errors")
error_tables <- error_tables[vapply(error_tables, nrow, integer(1L)) > 0L]
numerical_errors <- if (length(error_tables) > 0L) do.call(rbind, error_tables) else data.frame()

status_directory <- file.path(cmo_root(), "status", config$study_id)
dir.create(status_directory, recursive = TRUE, showWarnings = FALSE)
fatal_diagnostics <- list.files(
  status_directory,
  pattern = "^fatal_task_[0-9]{5}_job_.*\\.rds$",
  full.names = TRUE
)
fatal_objects <- lapply(fatal_diagnostics, function(path) {
  tryCatch(readRDS(path), error = function(error) NULL)
})
fatal_unreadable <- if (length(fatal_diagnostics) == 0L) {
  character(0L)
} else {
  fatal_diagnostics[vapply(fatal_objects, is.null, logical(1L))]
}
fatal_unresolved <- if (length(fatal_diagnostics) == 0L) {
  character(0L)
} else {
  unresolved <- vapply(seq_along(fatal_diagnostics), function(index) {
    object <- fatal_objects[[index]]
    if (is.null(object) || is.null(object$affected_unit_ids)) return(TRUE)
    affected_ids <- suppressWarnings(as.integer(object$affected_unit_ids))
    length(affected_ids) == 0L || anyNA(affected_ids) ||
      length(setdiff(affected_ids, unique(observed_ids))) > 0L
  }, logical(1L))
  fatal_diagnostics[unresolved]
}
fatal_recovered <- setdiff(fatal_diagnostics, fatal_unresolved)
done_markers <- list.files(
  status_directory,
  pattern = "^task_[0-9]{5}\\.done$",
  full.names = TRUE
)
done_task_ids <- suppressWarnings(as.integer(sub(
  "^task_([0-9]{5})\\.done$", "\\1", basename(done_markers)
)))
premature_done_tasks <- sort(intersect(done_task_ids, missing_task_ids))
complete_task_ids <- setdiff(tasks$task_id, missing_task_ids)
missing_done_markers <- sort(setdiff(complete_task_ids, done_task_ids))

missing_table <- units[units$unit_id %in% missing_ids, , drop = FALSE]
atomic_write_csv(missing_table, file.path(status_directory, "missing_units.csv"))
atomic_write_lines(as.character(missing_task_ids),
                   file.path(status_directory, "missing_task_ids.txt"))
atomic_write_csv(
  data.frame(path = unreadable_paths, stringsAsFactors = FALSE),
  file.path(status_directory, "unreadable_shards.csv")
)
if (nrow(numerical_errors) > 0L) {
  atomic_write_csv(numerical_errors, file.path(status_directory, "numerical_failures.csv"))
} else {
  atomic_write_csv(data.frame(), file.path(status_directory, "numerical_failures.csv"))
}

validation <- list(
  schema_version = config$schema_version,
  validation_schema_version = "global-queue-validation-v2",
  study_id = config$study_id,
  timestamp_utc = timestamp_utc(),
  expected_units = nrow(units),
  expected_tasks = nrow(tasks),
  independent_sample_streams = nrow(sample_registry),
  independent_bootstrap_streams = nrow(bootstrap_registry),
  observed_unique_units = length(unique(observed_ids)),
  missing_units = length(missing_ids),
  missing_task_ids = missing_task_ids,
  duplicate_units = duplicate_ids,
  unexpected_units = unexpected_ids,
  shard_count = length(shard_paths),
  valid_shards = length(valid_paths),
  unreadable_shards = unreadable_paths,
  incompatible_shards = incompatible_count,
  malformed_shards = malformed_count,
  duplicated_task_wave_keys = duplicated_task_wave,
  numerical_failures = nrow(numerical_errors),
  numerically_clean = nrow(numerical_errors) == 0L,
  fatal_diagnostics = fatal_diagnostics,
  fatal_diagnostics_total = length(fatal_diagnostics),
  fatal_diagnostics_unreadable = fatal_unreadable,
  fatal_diagnostics_unresolved = fatal_unresolved,
  fatal_diagnostics_unresolved_count = length(fatal_unresolved),
  fatal_diagnostics_recovered = fatal_recovered,
  fatal_diagnostics_recovered_count = length(fatal_recovered),
  premature_completion_markers = premature_done_tasks,
  missing_completion_markers = missing_done_markers,
  complete = length(missing_ids) == 0L &&
    length(duplicate_ids) == 0L &&
    length(unexpected_ids) == 0L &&
    length(unreadable_paths) == 0L &&
    incompatible_count == 0L &&
    malformed_count == 0L &&
    length(duplicated_task_wave) == 0L &&
    length(fatal_unresolved) == 0L &&
    length(premature_done_tasks) == 0L
)
atomic_save_rds(validation, file.path(status_directory, "validation.rds"), compress = "xz")

log_message(
  "Validation study=", config$study_id,
  ": expected=", validation$expected_units,
  ", observed=", validation$observed_unique_units,
  ", missing=", validation$missing_units,
  ", duplicates=", length(validation$duplicate_units),
  ", numerical_failures=", validation$numerical_failures,
  ", numerically_clean=", validation$numerically_clean,
  ", malformed_shards=", validation$malformed_shards,
  ", fatal_diagnostics=", validation$fatal_diagnostics_total,
  ", unresolved_fatal_diagnostics=",
  length(validation$fatal_diagnostics_unresolved)
)
if (length(missing_done_markers) > 0L && validation$complete) {
  log_message(
    "Note: ", length(missing_done_markers),
    " completion marker(s) are absent, but all corresponding atomic checkpoints are complete."
  )
}

if (!validation$complete) {
  quit(save = "no", status = 2L)
}
