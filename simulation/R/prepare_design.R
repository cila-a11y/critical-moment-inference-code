#!/usr/bin/env Rscript

source(file.path(Sys.getenv("CMO_ROOT", unset = "."), "R", "common.R"))
source_cmo("R/dgp.R")

arguments <- parse_named_args()
config_argument <- arguments$config %||% ""
assert_true(nzchar(config_argument), "Supply --config=... explicitly.")
config <- read_config(config_argument)

required_names <- c(
  "schema_version", "study_id", "release_status", "master_seed", "alpha",
  "bootstrap_reps", "bootstrap_batch_size", "multiplier_distribution",
  "variance_tolerance", "contrast_variance_tolerance",
  "root_variance_tolerance", "augmented_variance_tolerance",
  "slope_tolerance", "numerical_tolerance",
  "roundoff_inflation_required", "run_derivative_assisted_ablation",
  "run_grid_only_ablation", "run_pointwise_ablation",
  "run_wald_refinement", "run_root_multiplier_refinement",
  "maximum_enclosure_bisections", "maximum_enclosure_nodes",
  "maximum_truth_bisection_levels", "maximum_truth_subintervals",
  "mc_reps_per_cell",
  "units_per_task", "units_per_worker_chunk", "maximum_workers",
  "checkpoint_reserve_seconds", "internal_wall_budget_seconds"
)
missing_names <- setdiff(required_names, names(config))
assert_true(
  length(missing_names) == 0L,
  sprintf("The configuration is missing: %s", paste(missing_names, collapse = ", "))
)
if (!identical(config$release_status, "ready")) {
  stopf(
    "Configuration '%s' has release_status='%s' and is not released for design generation.",
    config$config_path,
    config$release_status
  )
}
assert_true(valid_identifier(config$study_id), "study_id contains an unsafe character.")
assert_true(
  length(config$master_seed) == 1L && is.finite(config$master_seed) &&
    config$master_seed == as.integer(config$master_seed) &&
    config$master_seed >= 1 && config$master_seed <= .Machine$integer.max,
  "master_seed must be a positive 32-bit integer."
)

integer_fields <- intersect(
  c(
    "grid_size", "bootstrap_grid_size", "enclosure_grid_size",
    "audit_grid_size", "root_scan_size", "bootstrap_reps",
    "bootstrap_batch_size", "mc_reps_per_cell", "units_per_task",
    "maximum_workers", "maximum_enclosure_bisections",
    "maximum_continuum_enclosure_bisections",
    "maximum_enclosure_nodes", "maximum_enclosure_subintervals",
    "maximum_truth_bisection_levels", "maximum_truth_subintervals",
    "maximum_root_bisections",
    "maximum_root_evaluations", "total_units", "units_per_worker_chunk"
  ),
  names(config)
)
for (field in integer_fields) {
  value <- config[[field]]
  assert_true(
    length(value) == 1L && is.finite(value) && value == as.integer(value),
    sprintf("%s must be a finite integer.", field)
  )
}
assert_true(config$mc_reps_per_cell >= 1L, "Each cell needs at least one replication.")
assert_true(config$units_per_task >= 1L, "units_per_task must be positive.")
assert_true(config$units_per_worker_chunk >= 1L &&
              config$units_per_worker_chunk <= 1024L,
            "units_per_worker_chunk must lie between one and 1024.")
assert_true(config$maximum_workers >= 1L && config$maximum_workers <= 127L,
            "maximum_workers must lie between one and 127.")
assert_true(
  is.finite(config$internal_wall_budget_seconds) &&
    config$internal_wall_budget_seconds > 0 &&
    config$internal_wall_budget_seconds <= 42000,
  "internal_wall_budget_seconds must lie in (0,42000]."
)
assert_true(
  is.finite(config$checkpoint_reserve_seconds) &&
    config$checkpoint_reserve_seconds >= 15 &&
    config$checkpoint_reserve_seconds < config$internal_wall_budget_seconds,
  "checkpoint_reserve_seconds must be at least 15 and below the internal budget."
)
assert_true(config$alpha > 0 && config$alpha < 1, "alpha must lie in (0,1).")
assert_true(config$bootstrap_reps >= 2L, "bootstrap_reps must be at least two.")
assert_true(config$multiplier_distribution %in% c("rademacher", "normal"),
            "The multiplier distribution must be rademacher or normal.")
if (!is.null(config$bootstrap_batch_size)) {
  assert_true(
    config$bootstrap_batch_size >= 1L &&
      config$bootstrap_batch_size <= config$bootstrap_reps,
    "bootstrap_batch_size must lie between one and bootstrap_reps."
  )
}
for (field in intersect(
  c(
    "variance_tolerance", "contrast_variance_tolerance",
    "root_variance_tolerance", "augmented_variance_tolerance",
    "slope_tolerance", "numerical_tolerance", "bootstrap_grid_spacing",
    "enclosure_grid_spacing", "audit_grid_spacing"
  ),
  names(config)
)) {
  assert_true(length(config[[field]]) == 1L && is.finite(config[[field]]) &&
                config[[field]] >= 0,
              sprintf("%s must be finite and nonnegative.", field))
}
logical_fields <- c(
  "roundoff_inflation_required", "run_derivative_assisted_ablation",
  "run_grid_only_ablation", "run_pointwise_ablation",
  "run_wald_refinement", "run_root_multiplier_refinement"
)
for (field in logical_fields) {
  assert_true(
    is.logical(config[[field]]) && length(config[[field]]) == 1L &&
      !is.na(config[[field]]),
    sprintf("%s must be one nonmissing logical value.", field)
  )
}

order_intervals <- normalize_order_intervals(config)
default_interval_ids <- names(order_intervals)
default_moment_types <- normalize_tokens(
  config$moment_types %||% character(0L),
  "config$moment_types",
  allow_empty = TRUE
)
assert_true(
  length(default_moment_types) > 0L || !is.null(config$cells),
  "Specify global moment_types or moment_types in every explicit cell."
)
if (length(default_moment_types) > 0L) {
  assert_true(all(default_moment_types %in% c("absolute", "signed")),
              "Only absolute and signed moment types are supported.")
}

as_cell_list <- function(config) {
  if (!is.null(config$cells)) {
    if (is.data.frame(config$cells)) {
      cells <- lapply(seq_len(nrow(config$cells)), function(index) {
        lapply(config$cells[index, , drop = FALSE], `[[`, 1L)
      })
      names(cells) <- row.names(config$cells)
      return(cells)
    }
    assert_true(is.list(config$cells) && length(config$cells) > 0L,
                "cells must be a nonempty list or data frame.")
    return(config$cells)
  }

  assert_true(length(config$scenario_ids) > 0L,
              "A legacy configuration needs scenario_ids.")
  assert_true(length(config$sample_sizes) > 0L,
              "A legacy configuration needs sample_sizes.")
  cells <- list()
  counter <- 0L
  for (scenario_id in config$scenario_ids) {
    for (sample_size in config$sample_sizes) {
      counter <- counter + 1L
      cells[[counter]] <- list(
        scenario_id = scenario_id,
        n_x = sample_size[["n_x"]],
        n_y = sample_size[["n_y"]],
        interval_ids = default_interval_ids,
        moment_types = default_moment_types
      )
    }
  }
  cells
}

raw_cells <- as_cell_list(config)
raw_cell_names <- names(raw_cells)
if (is.null(raw_cell_names)) {
  raw_cell_names <- rep("", length(raw_cells))
}
variants <- config$variants %||% config$analysis_variants %||% list()
assert_true(is.list(variants), "variants must be a named list when supplied.")

canonical_fields <- c(
  "cell_id", "scenario_id", "n_x", "n_y", "interval_id", "interval_ids",
  "moment_type", "moment_types", "pairing_key", "multiplier_pairing_key",
  "variant_id", "cell_class"
)
normalized_cells <- vector("list", length(raw_cells))
for (cell_index in seq_along(raw_cells)) {
  raw_cell <- raw_cells[[cell_index]]
  assert_true(is.list(raw_cell), sprintf("Cell %d must be a list.", cell_index))

  variant_id <- as.character(raw_cell$variant_id %||% "primary")
  assert_true(valid_identifier(variant_id),
              sprintf("Cell %d has an unsafe variant_id.", cell_index))
  variant_spec <- list()
  if (length(variants) > 0L && variant_id %in% names(variants)) {
    variant_spec <- variants[[variant_id]]
    assert_true(is.list(variant_spec),
                sprintf("Variant %s must be a list.", variant_id))
  }
  effective <- merge_named_lists(variant_spec, raw_cell)

  for (field in intersect(
    c("bootstrap_grid_spacing", "enclosure_grid_spacing"),
    names(effective)
  )) {
    value <- effective[[field]]
    assert_true(
      length(value) == 1L && is.finite(value) && value >= 0,
      sprintf("cells[[%d]]$%s must be finite and nonnegative.",
              cell_index, field)
    )
  }
  for (field in intersect(
    c(
      "grid_size", "bootstrap_grid_size", "enclosure_grid_size",
      "maximum_enclosure_bisections",
      "maximum_continuum_enclosure_bisections"
    ),
    names(effective)
  )) {
    lower <- if (field %in% c(
      "enclosure_grid_size", "maximum_enclosure_bisections",
      "maximum_continuum_enclosure_bisections"
    )) 0L else 3L
    safe_integer(
      effective[[field]],
      sprintf("cells[[%d]]$%s", cell_index, field),
      lower = lower
    )
  }

  scenario_id <- as.character(effective$scenario_id %||% "")
  assert_true(valid_identifier(scenario_id),
              sprintf("Cell %d has an invalid scenario_id.", cell_index))
  get_scenario(scenario_id)
  n_x <- safe_integer(effective$n_x, sprintf("cells[[%d]]$n_x", cell_index), lower = 2L)
  n_y <- safe_integer(effective$n_y, sprintf("cells[[%d]]$n_y", cell_index), lower = 2L)

  interval_ids <- effective$interval_ids %||% effective$interval_id %||%
    default_interval_ids
  interval_ids <- normalize_tokens(
    interval_ids,
    sprintf("cells[[%d]]$interval_ids", cell_index)
  )
  assert_true(all(interval_ids %in% default_interval_ids),
              sprintf("Cell %d refers to an unknown order interval.", cell_index))

  moment_types <- effective$moment_types %||% effective$moment_type %||%
    default_moment_types
  moment_types <- normalize_tokens(
    moment_types,
    sprintf("cells[[%d]]$moment_types", cell_index)
  )
  assert_true(all(moment_types %in% c("absolute", "signed")),
              sprintf("Cell %d has an unsupported moment type.", cell_index))

  suggested_name <- raw_cell_names[cell_index]
  if (!nzchar(suggested_name)) {
    suggested_name <- sprintf(
      "%03d_%s_nx%d_ny%d_%s",
      cell_index, scenario_id, n_x, n_y, variant_id
    )
  }
  cell_id <- as.character(effective$cell_id %||% suggested_name)
  assert_true(valid_identifier(cell_id),
              sprintf("Cell %d has an unsafe cell_id: %s", cell_index, cell_id))
  pairing_key <- as.character(effective$pairing_key %||% cell_id)
  assert_true(valid_identifier(pairing_key),
              sprintf("Cell %d has an unsafe pairing_key.", cell_index))
  multiplier_pairing_key <- as.character(
    effective$multiplier_pairing_key %||% pairing_key
  )
  assert_true(valid_identifier(multiplier_pairing_key),
              sprintf("Cell %d has an unsafe multiplier_pairing_key.", cell_index))
  cell_class <- as.character(effective$cell_class %||% "unclassified")
  allowed_cell_classes <- c(
    "benchmark", "stress", "power", "diagnostic", "operational",
    "unclassified"
  )
  assert_true(
    length(cell_class) == 1L && !is.na(cell_class) &&
      cell_class %in% allowed_cell_classes,
    sprintf(
      "Cell %d has cell_class='%s'; expected benchmark, stress, power, diagnostic, operational, or unclassified.",
      cell_index, paste(cell_class, collapse = ";")
    )
  )

  analysis_overrides <- effective[setdiff(names(effective), canonical_fields)]
  normalized_cells[[cell_index]] <- list(
    cell_index = as.integer(cell_index),
    cell_id = cell_id,
    scenario_id = scenario_id,
    n_x = n_x,
    n_y = n_y,
    interval_ids = interval_ids,
    moment_types = moment_types,
    pairing_key = pairing_key,
    multiplier_pairing_key = multiplier_pairing_key,
    variant_id = variant_id,
    cell_class = cell_class,
    analysis_overrides = analysis_overrides
  )
}

cell_ids <- vapply(normalized_cells, `[[`, character(1L), "cell_id")
assert_true(!anyDuplicated(cell_ids), "cell_id values must be unique.")
if (!is.null(config$paired_contrasts)) {
  pair_specs <- config$paired_contrasts
  pair_columns <- c(
    "contrast_id", "contrast_definition", "reference_variant_id",
    "comparison_variant_id"
  )
  assert_true(is.data.frame(pair_specs) &&
                all(pair_columns %in% names(pair_specs)) &&
                nrow(pair_specs) > 0L,
              "paired_contrasts must be a nonempty data frame with the required columns.")
  assert_true(!anyDuplicated(as.character(pair_specs$contrast_id)) &&
                all(nzchar(as.character(pair_specs$contrast_id))),
              "paired_contrasts needs unique nonempty contrast_id values.")
  declared_variants <- unique(c(
    as.character(pair_specs$reference_variant_id),
    as.character(pair_specs$comparison_variant_id)
  ))
  available_variants <- vapply(
    normalized_cells, `[[`, character(1L), "variant_id"
  )
  assert_true(all(declared_variants %in% available_variants),
              "paired_contrasts refers to an absent variant_id.")
  assert_true(all(as.character(pair_specs$reference_variant_id) !=
                    as.character(pair_specs$comparison_variant_id)),
              "A paired contrast cannot compare a variant with itself.")
}
pairing_keys <- vapply(normalized_cells, `[[`, character(1L), "pairing_key")
pairing_groups <- split(seq_along(normalized_cells), pairing_keys)
for (pairing_key in names(pairing_groups)) {
  indices <- pairing_groups[[pairing_key]]
  signatures <- vapply(indices, function(index) {
    cell <- normalized_cells[[index]]
    paste(cell$scenario_id, cell$n_x, cell$n_y, sep = "::")
  }, character(1L))
  assert_true(
    length(unique(signatures)) == 1L,
    sprintf(
      "Cells sharing pairing_key '%s' must have the same scenario_id, n_x, and n_y.",
      pairing_key
    )
  )
}
multiplier_pairing_keys <- vapply(
  normalized_cells, `[[`, character(1L), "multiplier_pairing_key"
)
multiplier_pairing_groups <- split(seq_along(normalized_cells), multiplier_pairing_keys)
for (pairing_key in names(multiplier_pairing_groups)) {
  indices <- multiplier_pairing_groups[[pairing_key]]
  signatures <- vapply(indices, function(index) {
    cell <- normalized_cells[[index]]
    paste(cell$n_x, cell$n_y, sep = "::")
  }, character(1L))
  assert_true(
    length(unique(signatures)) == 1L,
    sprintf(
      "Cells sharing multiplier_pairing_key '%s' must have the same n_x and n_y.",
      pairing_key
    )
  )
}

design_directory <- file.path(cmo_root(), "design", config$study_id)
dir.create(design_directory, recursive = TRUE, showWarnings = FALSE)
artifact_paths <- design_artifact_paths(design_directory)
design_fingerprint_path <- file.path(design_directory, "design_fingerprint.rds")
all_frozen_paths <- c(artifact_paths, design_fingerprint_path)
if (any(file.exists(all_frozen_paths)) && !all(file.exists(all_frozen_paths))) {
  stop(
    "A partial frozen design already exists. Never repair it in place; use a new study_id.",
    call. = FALSE
  )
}
if (all(file.exists(all_frozen_paths))) {
  frozen_config <- readRDS(file.path(design_directory, "config_snapshot.rds"))
  frozen_fingerprint <- readRDS(file.path(design_directory, "fingerprint.rds"))
  frozen_design_fingerprint <- readRDS(design_fingerprint_path)
  current_fingerprint <- project_fingerprint(config$config_path)
  if (identical(frozen_config, config) &&
      fingerprints_identical(frozen_fingerprint, current_fingerprint) &&
      fingerprints_identical(
        frozen_design_fingerprint,
        compute_design_fingerprint(design_directory)
      )) {
    log_message("The frozen design for ", config$study_id, " already exists and is unchanged.")
    quit(save = "no", status = 0L)
  }
  stop(
    "This study_id already has a different frozen design or code fingerprint. Use a new study_id.",
    call. = FALSE
  )
}

# Replication-major ordering interleaves cells before globally packing units
# into Slurm tasks.  Task boundaries therefore never define the RNG streams.
number_of_cells <- length(normalized_cells)
number_of_replications <- as.integer(config$mc_reps_per_cell)
cell_sequence <- rep(seq_len(number_of_cells), times = number_of_replications)
rep_sequence <- rep(seq_len(number_of_replications), each = number_of_cells)
cell_column <- function(name, FUN, FUN.VALUE) {
  values <- vapply(normalized_cells, `[[`, FUN.VALUE, name)
  FUN(values[cell_sequence])
}
cell_interval_values <- vapply(
  normalized_cells, function(cell) collapse_tokens(cell$interval_ids), character(1L)
)
cell_moment_values <- vapply(
  normalized_cells, function(cell) collapse_tokens(cell$moment_types), character(1L)
)
units <- data.frame(
  unit_id = seq_along(cell_sequence),
  queue_position = seq_along(cell_sequence),
  cell_index = cell_sequence,
  cell_id = cell_column("cell_id", as.character, character(1L)),
  cell_class = cell_column("cell_class", as.character, character(1L)),
  scenario_id = cell_column("scenario_id", as.character, character(1L)),
  n_x = cell_column("n_x", as.integer, integer(1L)),
  n_y = cell_column("n_y", as.integer, integer(1L)),
  rep_id = rep_sequence,
  pairing_key = cell_column("pairing_key", as.character, character(1L)),
  multiplier_pairing_key = cell_column(
    "multiplier_pairing_key", as.character, character(1L)
  ),
  variant_id = cell_column("variant_id", as.character, character(1L)),
  interval_ids = cell_interval_values[cell_sequence],
  moment_types = cell_moment_values[cell_sequence],
  stringsAsFactors = FALSE
)
if (!is.null(config$total_units)) {
  assert_true(
    safe_integer(config$total_units, "total_units", lower = 1L) == nrow(units),
    sprintf("Configured total_units does not equal the generated total %d.", nrow(units))
  )
}
if (identical(config$task_packing %||% "global", "all_cells_single_task")) {
  assert_true(
    config$units_per_task >= nrow(units),
    "task_packing='all_cells_single_task' requires units_per_task >= total_units."
  )
}

sample_seed_keys <- paste(units$pairing_key, units$rep_id, sep = "::rep=")
unique_sample_seed_keys <- unique(sample_seed_keys)
units$seed_id <- match(sample_seed_keys, unique_sample_seed_keys)
bootstrap_seed_keys <- paste(
  units$multiplier_pairing_key, units$rep_id, sep = "::rep="
)
unique_bootstrap_seed_keys <- unique(bootstrap_seed_keys)
units$bootstrap_seed_id <- match(bootstrap_seed_keys, unique_bootstrap_seed_keys)

task_count <- ceiling(nrow(units) / as.integer(config$units_per_task))
base_task_size <- nrow(units) %/% task_count
extra_units <- nrow(units) %% task_count
task_sizes <- rep(base_task_size, task_count)
if (extra_units > 0L) {
  task_sizes[seq_len(extra_units)] <- task_sizes[seq_len(extra_units)] + 1L
}
assert_true(max(task_sizes) <= config$units_per_task,
            "Internal error: a packed task exceeds units_per_task.")
units$task_id <- rep(seq_len(task_count), times = task_sizes)
units <- units[, c(
  "unit_id", "seed_id", "task_id", "queue_position", "cell_index",
  "bootstrap_seed_id", "cell_id", "cell_class", "scenario_id", "n_x", "n_y", "rep_id",
  "pairing_key", "multiplier_pairing_key", "variant_id", "interval_ids",
  "moment_types"
)]

task_rows <- lapply(seq_len(task_count), function(task_id) {
  subset <- units[units$task_id == task_id, , drop = FALSE]
  data.frame(
    task_id = task_id,
    first_queue_position = min(subset$queue_position),
    last_queue_position = max(subset$queue_position),
    first_unit = min(subset$unit_id),
    last_unit = max(subset$unit_id),
    number_of_units = nrow(subset),
    number_of_cells = length(unique(subset$cell_id)),
    number_of_seed_streams = length(unique(subset$seed_id)) +
      length(unique(subset$bootstrap_seed_id)),
    stringsAsFactors = FALSE
  )
})
tasks <- do.call(rbind, task_rows)
row.names(tasks) <- NULL

sample_seed_registry <- unique(units[, c("seed_id", "pairing_key", "rep_id")])
sample_seed_registry <- sample_seed_registry[order(sample_seed_registry$seed_id), , drop = FALSE]
names(sample_seed_registry)[names(sample_seed_registry) == "pairing_key"] <- "stream_key"
sample_seed_registry$stream_role <- "sample"
bootstrap_seed_registry <- unique(
  units[, c("bootstrap_seed_id", "multiplier_pairing_key", "rep_id")]
)
bootstrap_seed_registry <- bootstrap_seed_registry[
  order(bootstrap_seed_registry$bootstrap_seed_id), , drop = FALSE
]
names(bootstrap_seed_registry) <- c("seed_id", "stream_key", "rep_id")
bootstrap_seed_registry$stream_role <- "bootstrap"
seed_registry <- rbind(sample_seed_registry, bootstrap_seed_registry)
seed_registry <- seed_registry[, c("stream_role", "seed_id", "stream_key", "rep_id")]
assert_true(
  identical(sample_seed_registry$seed_id, seq_len(nrow(sample_seed_registry))) &&
    identical(bootstrap_seed_registry$seed_id, seq_len(nrow(bootstrap_seed_registry))),
  "The sample and bootstrap seed registries must each be contiguous."
)

RNGkind(
  kind = "L'Ecuyer-CMRG",
  normal.kind = "Inversion",
  sample.kind = "Rejection"
)
set.seed(as.integer(config$master_seed))
sample_seeds <- matrix(
  NA_integer_, nrow = length(.Random.seed), ncol = nrow(sample_seed_registry)
)
current_seed <- .Random.seed
for (seed_id in seq_len(ncol(sample_seeds))) {
  sample_seeds[, seed_id] <- current_seed
  current_seed <- parallel::nextRNGStream(current_seed)
}
bootstrap_seeds <- matrix(
  NA_integer_, nrow = length(.Random.seed), ncol = nrow(bootstrap_seed_registry)
)
for (seed_id in seq_len(ncol(bootstrap_seeds))) {
  bootstrap_seeds[, seed_id] <- current_seed
  current_seed <- parallel::nextRNGStream(current_seed)
}
colnames(sample_seeds) <- as.character(sample_seed_registry$seed_id)
colnames(bootstrap_seeds) <- as.character(bootstrap_seed_registry$seed_id)
seeds <- list(sample = sample_seeds, bootstrap = bootstrap_seeds)

collapse_truth_value <- function(value, digits = NULL) {
  if (is.null(value) || length(value) == 0L) {
    return("")
  }
  if (is.numeric(value) && !is.null(digits)) {
    value <- format(value, digits = digits, trim = TRUE, scientific = FALSE)
  }
  paste(as.character(value), collapse = ";")
}

truth_rows <- list()
truth_counter <- 0L
for (cell in normalized_cells) {
  scenario <- get_scenario(cell$scenario_id)
  for (moment_type in cell$moment_types) {
    for (interval_id in cell$interval_ids) {
      bounds <- order_intervals[[interval_id]]
      truth_formals <- names(formals(truth_on_interval))
      truth_arguments <- list(
        scenario = scenario,
        p_min = bounds[1L],
        p_max = bounds[2L],
        interval = bounds,
        moment_type = moment_type
      )
      truth <- do.call(
        truth_on_interval,
        truth_arguments[names(truth_arguments) %in% truth_formals]
      )
      type_root_truth <- scenario$root_truth[[moment_type]]
      type_root_truth <- type_root_truth[
        type_root_truth$order >= bounds[1L] &
          type_root_truth$order <= bounds[2L],
        ,
        drop = FALSE
      ]
      truth_roots <- truth$roots %||% type_root_truth$order
      truth_multiplicities <- truth$root_multiplicities %||%
        type_root_truth$multiplicity
      truth_directions <- truth$root_directions %||% type_root_truth$direction
      distinct_root_count <- truth$root_count %||% nrow(type_root_truth)
      total_root_multiplicity <- truth$total_multiplicity %||%
        sum(truth_multiplicities)
      type_structural_budget <- scenario$structural_budgets[[moment_type]] %||%
        truth$structural_budget %||% NA_integer_
      type_tail_crossing_scales <- scenario$tail_crossing_scales[[moment_type]]
      truth_counter <- truth_counter + 1L
      truth_rows[[truth_counter]] <- data.frame(
        cell_id = cell$cell_id,
        cell_class = cell$cell_class,
        scenario_id = cell$scenario_id,
        variant_id = cell$variant_id,
        moment_type = moment_type,
        interval_id = interval_id,
        p_min = bounds[1L],
        p_max = bounds[2L],
        label = as.character(scenario$label %||% cell$scenario_id),
        structural_budget = collapse_truth_value(type_structural_budget),
        root_direction = collapse_truth_value(truth_directions),
        roots_in_interval = collapse_truth_value(truth_roots, digits = 16L),
        root_count_in_interval = collapse_truth_value(distinct_root_count),
        distinct_root_count_in_interval = collapse_truth_value(
          distinct_root_count
        ),
        total_multiplicity_in_interval = collapse_truth_value(
          total_root_multiplicity
        ),
        root_multiplicities = collapse_truth_value(truth_multiplicities),
        tail_crossing_scale = collapse_truth_value(
          type_tail_crossing_scales,
          digits = 16L
        ),
        stringsAsFactors = FALSE
      )
    }
  }
}
truth_table <- do.call(rbind, truth_rows)

cells_object <- list(
  schema_version = config$schema_version,
  study_id = config$study_id,
  order_intervals = order_intervals,
  cells = normalized_cells
)
fingerprint <- project_fingerprint(config$config_path)
atomic_write_csv(units, file.path(design_directory, "units.csv"), overwrite = FALSE)
atomic_write_csv(tasks, file.path(design_directory, "tasks.csv"), overwrite = FALSE)
atomic_save_rds(list(units = units, tasks = tasks),
                file.path(design_directory, "queue.rds"),
                compress = FALSE, overwrite = FALSE)
atomic_save_rds(cells_object, file.path(design_directory, "cells.rds"),
                compress = "xz", overwrite = FALSE)
atomic_write_csv(truth_table, file.path(design_directory, "truth.csv"), overwrite = FALSE)
atomic_write_csv(seed_registry, file.path(design_directory, "seed_registry.csv"), overwrite = FALSE)
atomic_save_rds(seeds, file.path(design_directory, "seeds.rds"),
                compress = "xz", overwrite = FALSE)
atomic_save_rds(config, file.path(design_directory, "config_snapshot.rds"),
                compress = "xz", overwrite = FALSE)
atomic_save_rds(fingerprint, file.path(design_directory, "fingerprint.rds"),
                compress = "xz", overwrite = FALSE)
design_fingerprint <- compute_design_fingerprint(design_directory)
atomic_save_rds(
  design_fingerprint,
  design_fingerprint_path,
  compress = "xz",
  overwrite = FALSE
)

log_message("Prepared study ", config$study_id)
log_message(
  "Cells: ", length(normalized_cells),
  "; tasks: ", nrow(tasks),
  "; units: ", nrow(units),
  "; sample RNG streams: ", nrow(sample_seed_registry),
  "; bootstrap RNG streams: ", nrow(bootstrap_seed_registry)
)
log_message(
  "Packed task sizes: ", min(tasks$number_of_units),
  "--", max(tasks$number_of_units),
  " units (configured maximum ", config$units_per_task, ")"
)
log_message("Design directory: ", design_directory)
