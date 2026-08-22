#!/usr/bin/env Rscript

source(file.path(Sys.getenv("CMO_ROOT", unset = "."), "R", "common.R"))
source_cmo("R/dgp.R")
source_cmo("R/method.R")
set_single_thread_math()

arguments <- parse_named_args()
config_argument <- arguments$config %||% Sys.getenv("CMO_CONFIG", unset = "")
assert_true(nzchar(config_argument), "Supply --config=... or set CMO_CONFIG.")
config <- read_config(config_argument)
task_value <- arguments$task %||% Sys.getenv("SLURM_ARRAY_TASK_ID", unset = "")
task_id <- safe_integer(task_value, "task", lower = 1L)

design_directory <- file.path(cmo_root(), "design", config$study_id)
design_paths <- c(
  design_artifact_paths(design_directory),
  file.path(design_directory, "design_fingerprint.rds")
)
for (path in design_paths) {
  assert_true(file.exists(path), sprintf("Missing frozen design file: %s", path))
}

stored_config <- readRDS(file.path(design_directory, "config_snapshot.rds"))
stored_fingerprint <- readRDS(file.path(design_directory, "fingerprint.rds"))
stored_design_fingerprint <- readRDS(file.path(design_directory, "design_fingerprint.rds"))
current_fingerprint <- project_fingerprint(config$config_path)
assert_true(identical(stored_config, config),
            "The supplied configuration differs from the frozen configuration.")
assert_true(
  fingerprints_identical(stored_fingerprint, current_fingerprint),
  "The code or configuration changed after design freezing. Create a new study_id."
)
assert_true(
  fingerprints_identical(
    stored_design_fingerprint,
    compute_design_fingerprint(design_directory)
  ),
  "A frozen design artifact changed after preparation."
)

queue <- readRDS(file.path(design_directory, "queue.rds"))
all_units <- queue$units
all_tasks <- queue$tasks
cells_object <- readRDS(file.path(design_directory, "cells.rds"))
seeds <- readRDS(file.path(design_directory, "seeds.rds"))
assert_true(task_id %in% all_tasks$task_id,
            sprintf("Task %d is absent from the frozen design.", task_id))
units <- all_units[all_units$task_id == task_id, , drop = FALSE]
assert_true(nrow(units) == all_tasks$number_of_units[match(task_id, all_tasks$task_id)],
            "Task metadata and unit allocation disagree.")
assert_true(!anyDuplicated(units$unit_id), "A task contains duplicate unit IDs.")
assert_true(
  is.list(seeds) && is.matrix(seeds$sample) && is.matrix(seeds$bootstrap) &&
    all(as.character(unique(units$seed_id)) %in% colnames(seeds$sample)) &&
    all(as.character(unique(units$bootstrap_seed_id)) %in% colnames(seeds$bootstrap)),
  "At least one sample or bootstrap seed is absent from the immutable registry."
)

cell_lookup <- cells_object$cells
names(cell_lookup) <- vapply(cell_lookup, `[[`, character(1L), "cell_id")
assert_true(all(units$cell_id %in% names(cell_lookup)),
            "A task refers to a cell absent from cells.rds.")

result_directory <- file.path(
  cmo_root(), "results", "raw", config$study_id, sprintf("task_%05d", task_id)
)
status_directory <- file.path(cmo_root(), "status", config$study_id)
dir.create(result_directory, recursive = TRUE, showWarnings = FALSE)
dir.create(status_directory, recursive = TRUE, showWarnings = FALSE)
completion_marker <- file.path(status_directory, sprintf("task_%05d.done", task_id))

empty_error_table <- function() {
  data.frame(
    unit_id = integer(0L),
    seed_id = integer(0L),
    bootstrap_seed_id = integer(0L),
    task_id = integer(0L),
    cell_id = character(0L),
    cell_class = character(0L),
    scenario_id = character(0L),
    rep_id = integer(0L),
    pairing_key = character(0L),
    multiplier_pairing_key = character(0L),
    variant_id = character(0L),
    error_class = character(0L),
    error_message = character(0L),
    error_call = character(0L),
    stringsAsFactors = FALSE
  )
}

existing_paths <- sort(list.files(
  result_directory,
  pattern = "^wave_[0-9]{5}\\.rds$",
  full.names = TRUE
))
completed_units <- integer(0L)
next_wave_id <- 1L
if (length(existing_paths) > 0L) {
  existing_objects <- lapply(existing_paths, function(path) {
    tryCatch(
      readRDS(path),
      error = function(error) stopf("Unreadable existing checkpoint %s: %s", path, conditionMessage(error))
    )
  })
  valid_existing <- vapply(existing_objects, function(object) {
    is.list(object) &&
      identical(object$schema_version, config$schema_version) &&
      identical(object$study_id, config$study_id) &&
      identical(as.integer(object$task_id), task_id) &&
      fingerprints_identical(object$fingerprint, stored_fingerprint) &&
      fingerprints_identical(object$design_fingerprint, stored_design_fingerprint)
  }, logical(1L))
  assert_true(all(valid_existing), "At least one existing checkpoint is incompatible.")
  existing_id_lists <- lapply(existing_objects, function(object) as.integer(object$unit_ids))
  all_existing_ids <- as.integer(unlist(existing_id_lists))
  assert_true(!anyDuplicated(all_existing_ids),
              "Existing checkpoints contain duplicate unit IDs.")
  assert_true(all(all_existing_ids %in% units$unit_id),
              "An existing checkpoint contains a unit outside this task.")
  partition_valid <- vapply(existing_objects, function(object) {
    result_ids <- if (is.data.frame(object$results) && nrow(object$results) > 0L) {
      unique(as.integer(object$results$unit_id))
    } else integer(0L)
    error_ids <- if (is.data.frame(object$errors) && nrow(object$errors) > 0L) {
      unique(as.integer(object$errors$unit_id))
    } else integer(0L)
    length(intersect(result_ids, error_ids)) == 0L &&
      setequal(c(result_ids, error_ids), as.integer(object$unit_ids))
  }, logical(1L))
  assert_true(all(partition_valid),
              "An existing checkpoint does not partition its units into results and errors.")
  completed_units <- all_existing_ids
  existing_wave_ids <- as.integer(sub(
    "^wave_([0-9]{5})\\.rds$", "\\1", basename(existing_paths)
  ))
  assert_true(!anyNA(existing_wave_ids) && !anyDuplicated(existing_wave_ids),
              "Existing checkpoint wave identifiers are malformed or duplicated.")
  next_wave_id <- max(existing_wave_ids) + 1L
}

pending <- units[!units$unit_id %in% completed_units, , drop = FALSE]
if (file.exists(completion_marker) && nrow(pending) > 0L) {
  stop("A completion marker exists although units are missing; inspect the task manually.",
       call. = FALSE)
}
if (nrow(pending) == 0L) {
  log_message("Task ", task_id, " was already complete; nothing to do.")
  if (!file.exists(completion_marker)) {
    atomic_write_lines(
      c(
        paste0("study_id=", config$study_id),
        paste0("task_id=", task_id),
        paste0("completed_utc=", timestamp_utc()),
        "recovered_from_complete_checkpoints=yes"
      ),
      completion_marker,
      overwrite = FALSE
    )
  }
  quit(save = "no", status = 0L)
}

allocated_cpus <- safe_integer(
  Sys.getenv("SLURM_CPUS_PER_TASK", unset = as.character(parallel::detectCores())),
  "SLURM_CPUS_PER_TASK",
  lower = 2L
)
worker_count <- min(
  as.integer(config$maximum_workers),
  allocated_cpus - 1L,
  nrow(pending)
)
chunk_setting <- Sys.getenv("CMO_UNITS_PER_WORKER_CHUNK", unset = "")
if (!nzchar(chunk_setting)) {
  chunk_setting <- as.character(config$units_per_worker_chunk %||% 32L)
}
units_per_worker_chunk <- safe_integer(
  chunk_setting,
  "CMO_UNITS_PER_WORKER_CHUNK",
  lower = 1L
)
assert_true(units_per_worker_chunk <= 1024L,
            "The per-worker checkpoint chunk cannot exceed 1024 units.")
superwave_capacity <- worker_count * units_per_worker_chunk
wall_budget <- as.numeric(config$internal_wall_budget_seconds)
checkpoint_reserve <- as.numeric(config$checkpoint_reserve_seconds %||% 60)
assert_true(is.finite(checkpoint_reserve) && checkpoint_reserve >= 15 &&
              checkpoint_reserve < wall_budget,
            "checkpoint_reserve_seconds must be at least 15 and below the internal wall budget.")
job_started <- proc.time()[["elapsed"]]
wave_times <- numeric(0L)

analysis_formals <- names(formals(analyse_one_replication))
assert_true(length(analysis_formals) > 0L,
            "analyse_one_replication has no declared interface.")
log_message(
  "Method interface: analyse_one_replication(",
  paste(analysis_formals, collapse = ", "), ")"
)

controlled_order_grid <- function(
    bounds, interval_endpoints, spacing = NULL, size = NULL,
    label = "analysis", minimum_size = 2L) {
  use_spacing <- !is.null(spacing) && length(spacing) == 1L &&
    is.finite(spacing) && spacing > 0
  if (!is.null(spacing)) {
    assert_true(
      length(spacing) == 1L && is.finite(spacing) && spacing >= 0,
      sprintf("%s_grid_spacing must be one finite nonnegative value.", label)
    )
  }
  use_size <- !is.null(size)
  if (use_size) {
    size <- safe_integer(
      size,
      sprintf("%s_grid_size", label),
      lower = as.integer(minimum_size)
    )
  }
  assert_true(use_spacing || use_size,
              sprintf("The %s grid needs a spacing or a size.", label))
  size_matches_spacing <- use_spacing && use_size &&
    abs(
      (bounds[2L] - bounds[1L]) / (size - 1L) - as.numeric(spacing)
    ) <= 64 * .Machine$double.eps * max(
      1, abs(bounds), abs(as.numeric(spacing))
    )
  grid <- if (use_size && (!use_spacing || size_matches_spacing)) {
    regular_interval_grid(bounds, size = size)
  } else {
    regular_interval_grid(bounds, spacing = as.numeric(spacing))
  }
  grid <- sort(unique(c(as.numeric(grid), as.numeric(interval_endpoints))))
  assert_true(
    length(grid) >= minimum_size && all(is.finite(grid)) &&
      all(grid > 0) && all(diff(grid) > 0),
    sprintf("The %s grid is not finite, positive, and strictly increasing.", label)
  )
  grid
}

arguments_for_unit <- function(unit) {
  cell <- cell_lookup[[as.character(unit$cell_id)]]
  interval_ids <- split_tokens(unit$interval_ids)
  moment_types <- split_tokens(unit$moment_types)
  selected_intervals <- cells_object$order_intervals[interval_ids]
  assert_true(length(selected_intervals) == length(interval_ids) &&
                !any(vapply(selected_intervals, is.null, logical(1L))),
              sprintf("Unit %d has an invalid interval selection.", unit$unit_id))

  parameters <- merge_named_lists(config, cell$analysis_overrides)
  scenario <- get_scenario(as.character(unit$scenario_id))
  parameters$scenario <- scenario
  parameters$n_x <- as.integer(unit$n_x)
  parameters$n_y <- as.integer(unit$n_y)
  parameters$order_intervals <- selected_intervals
  parameters$interval_ids <- interval_ids
  parameters$moment_types <- moment_types
  parameters$variant_id <- as.character(unit$variant_id)
  parameters$analysis_spec <- cell$analysis_overrides
  alias_values <- list(
    guard_max_levels = parameters$maximum_enclosure_bisections,
    enclosure_max_levels =
      parameters$maximum_continuum_enclosure_bisections %||%
        parameters$maximum_enclosure_bisections,
    enclosure_max_nodes = parameters$maximum_enclosure_nodes,
    enclosure_safety_margin = parameters$enclosure_absolute_tolerance,
    root_max_levels = parameters$maximum_root_bisections,
    root_max_evaluations = parameters$maximum_root_evaluations,
    root_tolerance = parameters$root_absolute_tolerance
  )
  for (alias_name in names(alias_values)) {
    if (is.null(parameters[[alias_name]]) && !is.null(alias_values[[alias_name]])) {
      parameters[[alias_name]] <- alias_values[[alias_name]]
    }
  }
  if (is.null(parameters$root_intervals) &&
      as.character(unit$scenario_id) %in% c(
        "TWO_ROOT", "TWO_ROOT_STRONG", "TWO_ROOT_SV_POWER"
      ) &&
      !is.null(parameters$two_root_intervals)) {
    parameters$root_intervals <- parameters$two_root_intervals
  }
  if ("orders" %in% analysis_formals) {
    interval_endpoints <- unlist(selected_intervals, use.names = FALSE)
    full_bounds <- range(interval_endpoints)
    parameters$orders <- controlled_order_grid(
      bounds = full_bounds,
      interval_endpoints = interval_endpoints,
      spacing = parameters$bootstrap_grid_spacing,
      size = parameters$bootstrap_grid_size %||% parameters$grid_size,
      label = sprintf("bootstrap for unit %d", unit$unit_id),
      minimum_size = 3L
    )
    if ("enclosure_orders" %in% analysis_formals) {
      enclosure_has_spacing <- !is.null(parameters$enclosure_grid_spacing) &&
        length(parameters$enclosure_grid_spacing) == 1L &&
        is.finite(parameters$enclosure_grid_spacing) &&
        parameters$enclosure_grid_spacing > 0
      enclosure_has_size <- !is.null(parameters$enclosure_grid_size)
      parameters$enclosure_orders <- if (
        enclosure_has_spacing || enclosure_has_size
      ) {
        controlled_order_grid(
          bounds = full_bounds,
          interval_endpoints = interval_endpoints,
          spacing = parameters$enclosure_grid_spacing,
          size = parameters$enclosure_grid_size,
          label = sprintf("enclosure for unit %d", unit$unit_id),
          minimum_size = 2L
        )
      } else {
        NULL
      }
    }
  }
  parameters
}

run_unit <- function(unit) {
  unit_id <- as.integer(unit$unit_id)
  seed_id <- as.integer(unit$seed_id)
  bootstrap_seed_id <- as.integer(unit$bootstrap_seed_id)
  stream <- seeds$sample[, as.character(seed_id)]
  bootstrap_stream <- seeds$bootstrap[, as.character(bootstrap_seed_id)]
  assert_true(length(stream) == nrow(seeds$sample),
              sprintf("Stored RNG stream %d has the wrong length.", seed_id))
  assert_true(length(bootstrap_stream) == nrow(seeds$bootstrap),
              sprintf("Stored bootstrap RNG stream %d has the wrong length.", bootstrap_seed_id))
  tryCatch({
    assign(".Random.seed", stream, envir = .GlobalEnv)
    parameters <- arguments_for_unit(unit)
    parameters$bootstrap_seed <- bootstrap_stream
    call_arguments <- parameters[names(parameters) %in% setdiff(analysis_formals, "...")]
    result <- do.call(analyse_one_replication, call_arguments)
    assert_true(is.data.frame(result) && nrow(result) > 0L,
                "analyse_one_replication must return a nonempty data frame.")
    assert_true(
      all(c("interval_id", "moment_type") %in% names(result)),
      "The method result must contain interval_id and moment_type."
    )
    actual_keys <- paste(
      as.character(result$interval_id),
      as.character(result$moment_type),
      sep = "::"
    )
    expected_keys <- as.vector(outer(
      names(parameters$order_intervals),
      parameters$moment_types,
      FUN = function(interval_id, moment_type) {
        paste(interval_id, moment_type, sep = "::")
      }
    ))
    assert_true(
      !anyNA(result$interval_id) && !anyNA(result$moment_type) &&
        !anyDuplicated(actual_keys) && setequal(actual_keys, expected_keys),
      paste0(
        "The method must return exactly one row for each selected ",
        "interval_id by moment_type combination."
      )
    )
    reserved <- c(
      "schema_version", "study_id", "unit_id", "seed_id", "task_id",
      "bootstrap_seed_id", "cell_id", "scenario_id", "rep_id", "pairing_key",
      "multiplier_pairing_key", "variant_id", "cell_class"
    )
    assert_true(length(intersect(names(result), reserved)) == 0L,
                "The method returned a reserved orchestration column.")
    metadata <- data.frame(
      schema_version = config$schema_version,
      study_id = config$study_id,
      unit_id = unit_id,
      seed_id = seed_id,
      bootstrap_seed_id = bootstrap_seed_id,
      task_id = task_id,
      cell_id = as.character(unit$cell_id),
      cell_class = as.character(unit$cell_class),
      scenario_id = as.character(unit$scenario_id),
      rep_id = as.integer(unit$rep_id),
      pairing_key = as.character(unit$pairing_key),
      multiplier_pairing_key = as.character(unit$multiplier_pairing_key),
      variant_id = as.character(unit$variant_id),
      stringsAsFactors = FALSE
    )
    metadata <- metadata[rep(1L, nrow(result)), , drop = FALSE]
    list(
      status = "success",
      unit_id = unit_id,
      result = cbind(metadata, result),
      error = NULL
    )
  }, error = function(error) {
    if (!is_expected_numerical_failure(error)) {
      stop(error)
    }
    error_call <- conditionCall(error)
    list(
      status = "numerical_error",
      unit_id = unit_id,
      result = NULL,
      error = data.frame(
        unit_id = unit_id,
        seed_id = seed_id,
        bootstrap_seed_id = bootstrap_seed_id,
        task_id = task_id,
        cell_id = as.character(unit$cell_id),
        cell_class = as.character(unit$cell_class),
        scenario_id = as.character(unit$scenario_id),
        rep_id = as.integer(unit$rep_id),
        pairing_key = as.character(unit$pairing_key),
        multiplier_pairing_key = as.character(unit$multiplier_pairing_key),
        variant_id = as.character(unit$variant_id),
        error_class = class(error)[1L],
        error_message = conditionMessage(error),
        error_call = if (is.null(error_call)) "" else paste(deparse(error_call), collapse = " "),
        stringsAsFactors = FALSE
      )
    )
  })
}

run_chunk <- function(indices, superwave) {
  lapply(indices, function(index) run_unit(superwave[index, , drop = FALSE]))
}

log_message(
  "Starting study=", config$study_id,
  ", task=", task_id,
  ", pending_units=", nrow(pending),
  ", workers=", worker_count,
  ", units_per_worker_chunk=", units_per_worker_chunk,
  ", superwave_capacity=", superwave_capacity
)

while (nrow(pending) > 0L) {
  elapsed <- proc.time()[["elapsed"]] - job_started
  predicted_next <- if (length(wave_times) == 0L) {
    0
  } else {
    1.5 * max(tail(wave_times, 3L))
  }
  if (elapsed + predicted_next + checkpoint_reserve > wall_budget) {
    log_message("Stopping before a superwave to preserve atomic checkpoints.")
    break
  }

  superwave_size <- min(superwave_capacity, nrow(pending))
  superwave <- pending[seq_len(superwave_size), , drop = FALSE]
  chunk_count <- min(worker_count, superwave_size)
  # Round-robin assignment mixes cells and sample sizes within persistent forks.
  chunk_indices <- split(
    seq_len(superwave_size),
    rep(seq_len(chunk_count), length.out = superwave_size)
  )
  wave_started_utc <- timestamp_utc()
  wave_started <- proc.time()[["elapsed"]]
  chunk_outputs <- parallel::mclapply(
    chunk_indices,
    FUN = run_chunk,
    superwave = superwave,
    mc.cores = chunk_count,
    mc.preschedule = TRUE,
    mc.set.seed = FALSE
  )
  wave_elapsed <- proc.time()[["elapsed"]] - wave_started
  wave_finished_utc <- timestamp_utc()
  wave_times <- c(wave_times, wave_elapsed)

  valid_chunk <- vapply(chunk_outputs, function(output) {
    is.list(output) && all(vapply(output, function(item) {
      is.list(item) && length(item$unit_id) == 1L &&
        (identical(item$status, "success") ||
           identical(item$status, "numerical_error"))
    }, logical(1L)))
  }, logical(1L))
  valid_outputs <- if (any(valid_chunk)) {
    unlist(chunk_outputs[valid_chunk], recursive = FALSE, use.names = FALSE)
  } else {
    list()
  }

  result_rows <- lapply(valid_outputs, `[[`, "result")
  result_rows <- result_rows[!vapply(result_rows, is.null, logical(1L))]
  result_table <- if (length(result_rows) > 0L) do.call(rbind, result_rows) else data.frame()
  error_rows <- lapply(valid_outputs, `[[`, "error")
  error_rows <- error_rows[!vapply(error_rows, is.null, logical(1L))]
  error_table <- if (length(error_rows) > 0L) do.call(rbind, error_rows) else empty_error_table()
  unit_ids <- if (length(valid_outputs) > 0L) {
    vapply(valid_outputs, `[[`, integer(1L), "unit_id")
  } else integer(0L)
  assert_true(!anyDuplicated(unit_ids), "A superwave returned duplicate unit IDs.")

  if (length(unit_ids) > 0L) {
    shard <- list(
      shard_schema_version = "global-queue-superwave-v2",
      schema_version = config$schema_version,
      study_id = config$study_id,
      task_id = task_id,
      wave_id = next_wave_id,
      unit_ids = unit_ids,
      fingerprint = stored_fingerprint,
      design_fingerprint = stored_design_fingerprint,
      started_utc = wave_started_utc,
      finished_utc = wave_finished_utc,
      wave_elapsed_seconds = wave_elapsed,
      worker_count = worker_count,
      worker_chunk_count = chunk_count,
      configured_units_per_worker_chunk = units_per_worker_chunk,
      results = result_table,
      errors = error_table,
      session = capture_session()
    )
    shard_path <- file.path(result_directory, sprintf("wave_%05d.rds", next_wave_id))
    atomic_save_rds(
      shard,
      shard_path,
      compress = FALSE,
      overwrite = FALSE
    )
    log_message(
      "Published ", basename(shard_path),
      ": units=", length(unit_ids),
      ", numerical_failures=", nrow(error_table),
      ", seconds=", sprintf("%.2f", wave_elapsed)
    )
    pending <- pending[!pending$unit_id %in% unit_ids, , drop = FALSE]
    next_wave_id <- next_wave_id + 1L
  }

  if (!all(valid_chunk)) {
    bad_indices <- unlist(chunk_indices[!valid_chunk], use.names = FALSE)
    failure_path <- file.path(
      status_directory,
      sprintf("fatal_task_%05d_job_%s.rds", task_id,
              Sys.getenv("SLURM_JOB_ID", unset = "local"))
    )
    atomic_save_rds(
      list(
        schema_version = config$schema_version,
        study_id = config$study_id,
        task_id = task_id,
        timestamp_utc = timestamp_utc(),
        affected_unit_ids = as.integer(superwave$unit_id[bad_indices]),
        raw_chunk_outputs = chunk_outputs[!valid_chunk],
        session = capture_session()
      ),
      failure_path,
      compress = "xz",
      overwrite = FALSE
    )
    stopf(
      "%d worker chunk(s) failed structurally; good chunks were checkpointed. Diagnostics: %s",
      sum(!valid_chunk), failure_path
    )
  }
}

if (nrow(pending) > 0L) {
  log_message("Task ", task_id, " remains incomplete: ", nrow(pending), " units pending.")
  quit(save = "no", status = 75L)
}

atomic_write_lines(
  c(
    paste0("study_id=", config$study_id),
    paste0("task_id=", task_id),
    paste0("completed_utc=", timestamp_utc()),
    paste0("completed_units=", nrow(units))
  ),
  completion_marker,
  overwrite = FALSE
)
log_message("Task ", task_id, " completed.")
