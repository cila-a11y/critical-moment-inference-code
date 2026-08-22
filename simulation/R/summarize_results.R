#!/usr/bin/env Rscript

# Study-wide summaries for the v2 schema.  K1/K2/K3 and sensitivity variants
# are never pooled.  Binomial outputs include event counts and Wilson intervals,
# so a zero observed rate is not reported as if it had zero uncertainty.

source(file.path(Sys.getenv("CMO_ROOT", unset = "."), "R", "common.R"))
source_cmo("R/reporting.R")

arguments <- parse_named_args()
config_argument <- arguments$config %||% Sys.getenv("CMO_CONFIG", unset = "")
assert_true(nzchar(config_argument), "Supply --config=... or set CMO_CONFIG.")
config <- read_config(config_argument)

status_directory <- file.path(cmo_root(), "status", config$study_id)
validation_path <- file.path(status_directory, "validation.rds")
assert_true(file.exists(validation_path), "Run validate_results.R first.")
validation <- readRDS(validation_path)
assert_true(isTRUE(validation$complete),
            "The result set is incomplete or internally inconsistent.")

raw_directory <- file.path(cmo_root(), "results", "raw", config$study_id)
shard_paths <- sort(list.files(
  raw_directory, pattern = "^wave_[0-9]{5}\\.rds$",
  recursive = TRUE, full.names = TRUE
))
assert_true(length(shard_paths) > 0L, "No result checkpoint was found.")
objects <- lapply(shard_paths, readRDS)
result_tables <- lapply(objects, `[[`, "results")
result_tables <- result_tables[vapply(result_tables, nrow, integer(1L)) > 0L]
results <- if (length(result_tables) > 0L) {
  do.call(rbind, result_tables)
} else {
  data.frame()
}
row.names(results) <- NULL

# Root errors exist only after a prespecified isolating interval succeeds.
# Promote self-describing names before any table or persisted replication file
# is produced.  The engine-level legacy columns remain available solely for
# backward compatibility with completed bridge.3 evidence.
if (nrow(results) > 0L && "root_error" %in% names(results) &&
    !"root_error_conditional_on_isolation" %in% names(results)) {
  results$root_error_conditional_on_isolation <- results$root_error
  if ("root_isolation_success" %in% names(results)) {
    results$root_error_conditional_on_isolation[
      is.na(results$root_isolation_success) |
        !as.logical(results$root_isolation_success)
    ] <- NA_real_
  }
}
if (nrow(results) > 0L && "wald_collection_bias" %in% names(results) &&
    !"wald_collection_bias_conditional_on_isolation" %in% names(results)) {
  results$wald_collection_bias_conditional_on_isolation <-
    results$wald_collection_bias
  if ("root_isolation_success" %in% names(results)) {
    results$wald_collection_bias_conditional_on_isolation[
      is.na(results$root_isolation_success) |
        !as.logical(results$root_isolation_success)
    ] <- NA_real_
  }
}
if (nrow(results) > 0L && "wald_collection_rmse" %in% names(results) &&
    !"wald_collection_rmse_conditional_on_isolation" %in% names(results)) {
  results$wald_collection_rmse_conditional_on_isolation <-
    results$wald_collection_rmse
  if ("root_isolation_success" %in% names(results)) {
    results$wald_collection_rmse_conditional_on_isolation[
      is.na(results$root_isolation_success) |
        !as.logical(results$root_isolation_success)
    ] <- NA_real_
  }
}

error_tables <- lapply(objects, `[[`, "errors")
error_tables <- error_tables[vapply(error_tables, nrow, integer(1L)) > 0L]
numerical_failures <- if (length(error_tables)) {
  do.call(rbind, error_tables)
} else {
  data.frame()
}

group_columns <- c(
  "cell_id", "cell_class", "variant_id", "scenario_id", "n_x", "n_y",
  "moment_type", "interval_id"
)
if (nrow(results) > 0L) {
  assert_true(all(group_columns %in% names(results)),
              "The result table does not have the v2 grouping schema.")
  assert_true(!anyDuplicated(results[, c("unit_id", "moment_type", "interval_id")]),
              "A unit has duplicate moment-type/order-interval rows.")
}

# Build the aggregation universe from the immutable design rather than from
# successful rows.  Consequently, a cell with no successful replication is
# retained with a zero success count instead of silently disappearing.
design_directory <- file.path(cmo_root(), "design", config$study_id)
cells_object <- readRDS(file.path(design_directory, "cells.rds"))
truth_table <- read.csv(
  file.path(design_directory, "truth.csv"), stringsAsFactors = FALSE
)
cell_metadata <- do.call(rbind, lapply(cells_object$cells, function(cell) {
  tangency_neighborhood <-
    cell$analysis_overrides$tangency_neighborhood %||%
    config$tangency_neighborhood %||% c(0.50, 1.45)
  tangency_neighborhood <- as.numeric(tangency_neighborhood)
  assert_true(
    length(tangency_neighborhood) == 2L &&
      all(is.finite(tangency_neighborhood)) &&
      tangency_neighborhood[2L] > tangency_neighborhood[1L],
    sprintf("Cell %s has an invalid tangency_neighborhood.", cell$cell_id)
  )
  data.frame(
    cell_id = cell$cell_id,
    cell_class = cell$cell_class %||% "unclassified",
    variant_id = cell$variant_id,
    scenario_id = cell$scenario_id,
    pairing_key = cell$pairing_key,
    n_x = cell$n_x,
    n_y = cell$n_y,
    configured_bootstrap_reps = as.integer(
      cell$analysis_overrides$bootstrap_reps %||% config$bootstrap_reps
    ),
    configured_tangency_neighborhood_width =
      tangency_neighborhood[2L] - tangency_neighborhood[1L],
    stringsAsFactors = FALSE
  )
}))
expected_groups <- merge(
  truth_table[, c(
    "cell_id", "variant_id", "scenario_id", "moment_type", "interval_id",
    "p_min", "p_max"
  )],
  cell_metadata,
  by = c("cell_id", "variant_id", "scenario_id"),
  all.x = TRUE,
  sort = FALSE
)
assert_true(
  nrow(expected_groups) > 0L &&
    !anyNA(expected_groups[, c("n_x", "n_y", "configured_bootstrap_reps")]),
  "The frozen truth table and cell metadata do not define every expected group."
)
expected_groups <- unique(expected_groups)

make_group_key <- function(data) {
  if (nrow(data) == 0L) return(character(0L))
  do.call(
    paste,
    c(lapply(group_columns, function(name) as.character(data[[name]])),
      list(sep = "::"))
  )
}
expected_group_keys <- make_group_key(expected_groups)
assert_true(!anyDuplicated(expected_group_keys),
            "The frozen design contains duplicate scientific groups.")
observed_group_keys <- if (nrow(results) > 0L) make_group_key(results) else character(0L)
assert_true(all(observed_group_keys %in% expected_group_keys),
            "At least one successful result belongs to an unexpected scientific group.")
groups <- setNames(lapply(expected_group_keys, function(key) {
  if (nrow(results) == 0L) return(results)
  results[observed_group_keys == key, , drop = FALSE]
}), expected_group_keys)
group_metadata <- setNames(lapply(seq_len(nrow(expected_groups)), function(index) {
  expected_groups[index, , drop = FALSE]
}), expected_group_keys)

wilson_interval <- function(successes, total, confidence = 0.95) {
  if (!is.finite(total) || total <= 0L) {
    return(c(lower = NA_real_, upper = NA_real_))
  }
  estimate <- successes / total
  z <- stats::qnorm(1 - (1 - confidence) / 2)
  denominator <- 1 + z^2 / total
  centre <- (estimate + z^2 / (2 * total)) / denominator
  radius <- z / denominator * sqrt(
    estimate * (1 - estimate) / total + z^2 / (4 * total^2)
  )
  c(lower = max(0, centre - radius), upper = min(1, centre + radius))
}

binary_summary <- function(values) {
  values <- as.logical(values)
  observed <- values[!is.na(values)]
  total <- length(observed)
  successes <- if (total) sum(observed) else 0L
  estimate <- if (total) successes / total else NA_real_
  interval <- wilson_interval(successes, total)
  c(
    successes = successes, denominator = total, estimate = estimate,
    mcse = if (total) sqrt(estimate * (1 - estimate) / total) else NA_real_,
    wilson_lower = interval[["lower"]],
    wilson_upper = interval[["upper"]]
  )
}

finite_summary <- function(values) {
  values <- suppressWarnings(as.numeric(values))
  finite <- values[is.finite(values)]
  c(
    rows = length(values), finite = length(finite),
    missing = sum(is.na(values)),
    positive_infinite = sum(is.infinite(values) & values > 0, na.rm = TRUE),
    negative_infinite = sum(is.infinite(values) & values < 0, na.rm = TRUE),
    mean = if (length(finite)) mean(finite) else NA_real_,
    sd = if (length(finite) >= 2L) stats::sd(finite) else NA_real_,
    mean_mcse = if (length(finite) >= 2L) {
      stats::sd(finite) / sqrt(length(finite))
    } else NA_real_,
    q05 = if (length(finite)) unname(stats::quantile(finite, 0.05)) else NA_real_,
    q50 = if (length(finite)) unname(stats::quantile(finite, 0.50)) else NA_real_,
    q95 = if (length(finite)) unname(stats::quantile(finite, 0.95)) else NA_real_
  )
}

column_or <- function(data, name, default = NA) {
  if (name %in% names(data)) data[[name]] else rep(default, nrow(data))
}

mean_finite_or_na <- function(values) {
  values <- suppressWarnings(as.numeric(values))
  values <- values[is.finite(values)]
  if (length(values)) mean(values) else NA_real_
}

rmse_finite_or_na <- function(values) {
  values <- suppressWarnings(as.numeric(values))
  values <- values[is.finite(values)]
  if (length(values)) sqrt(mean(values^2)) else NA_real_
}

binary_candidates <- c(
  "band_covers_truth_on_grid", "band_covers_truth_on_audit_grid",
  "continuum_band_covers_truth", "guard_failed",
  "roundoff_inflation_requested", "roundoff_inflation_margin_applied",
  "continuum_truth_limit_hit", "continuum_truth_endpoint_violation",
  "enclosure_limit_hit", "enclosure_depth_limit_hit",
  "enclosure_node_limit_hit", "enclosure_variance_limit_hit",
  "enclosure_statistical_limit_hit", "continuum_no_root_certified",
  "grid_only_no_root", "grid_continuum_no_root_disagreement",
  "reversal_certified", "grid_only_reversal_certified",
  "direction_correct", "false_direction_report", "false_sign_report",
  "false_reversal_report", "false_no_root_report", "budget_rejected",
  "false_budget_rejection", "zero_budget_diagnostic_applicable",
  "zero_budget_rejected", "exact_count_certified", "simplicity_certified",
  "exact_count_statement_true", "false_exact_count_report",
  "bracket_covers_root", "brackets_cover_all_distinct_roots",
  "false_bracket_report", "outer_band_contains_root",
  "outer_set_covers_truth", "numerical_outer_contains_root",
  "numerical_outer_contains_all_roots", "false_outer_report",
  "tangency_localization_applicable", "tangency_local_outer_nonempty",
  "tangency_truncated_local_outer_nonempty",
  "tangency_truncated_local_touches_left_boundary",
  "tangency_truncated_local_touches_right_boundary",
  "tangency_truncated_local_censored_by_window",
  "tangency_global_geometry_applicable",
  "tangency_global_geometry_conditioned_on_continuum_band_coverage",
  "tangency_global_outer_nonempty_conditional_on_continuum_band_coverage",
  "tangency_global_outer_contains_root_conditional_on_continuum_band_coverage",
  "tangency_global_outer_touches_K_left_conditional_on_continuum_band_coverage",
  "tangency_global_outer_touches_K_right_conditional_on_continuum_band_coverage",
  "tangency_global_outer_touches_K_boundary_conditional_on_continuum_band_coverage",
  "tangency_global_root_component_present_conditional_on_continuum_band_coverage",
  "tangency_global_root_component_touches_K_left_conditional_on_continuum_band_coverage",
  "tangency_global_root_component_touches_K_right_conditional_on_continuum_band_coverage",
  "tangency_global_root_component_touches_K_boundary_conditional_on_continuum_band_coverage",
  "root_isolation_applicable", "root_isolation_success",
  "root_isolation_limit_hit", "root_isolation_all_tolerances_met",
  "wald_reported", "wald_conditional_cover",
  "wald_formal_unconditional_cover", "wald_report_and_cover",
  "root_multiplier_applicable", "root_multiplier_reported",
  "root_multiplier_conditional_cover",
  "root_multiplier_formal_unconditional_cover",
  "root_multiplier_report_and_cover", "derivative_assisted_requested",
  "derivative_assisted_implemented", "derivative_assisted_applicable",
  "derivative_assisted_reported", "derivative_assisted_failure",
  "derivative_assisted_same_multiplier_seed",
  "derivative_assisted_guard_failed",
  "derivative_assisted_enclosure_limit_hit",
  "derivative_assisted_derivative_positive_throughout_K",
  "derivative_assisted_derivative_negative_throughout_K",
  "derivative_assisted_exact_root_certified",
  "derivative_assisted_local_exact_root_certified",
  "derivative_assisted_exact_root_on_K_certified",
  "derivative_assisted_exact_root_on_K_statement_true",
  "derivative_assisted_false_exact_root_on_K_report",
  "derivative_assisted_certificate_reported",
  "derivative_assisted_certificate_statement_true",
  "derivative_assisted_false_exact_root_report",
  "derivative_assisted_value_grid_coverage",
  "derivative_assisted_derivative_grid_coverage",
  "derivative_assisted_joint_grid_coverage",
  "pointwise_grid_curve_coverage", "pointwise_root_coverage",
  "pointwise_false_sign_report", "pointwise_reversal_report",
  "pointwise_no_root_report", "any_primary_false_report",
  "joint_band_covers_truth_on_grid",
  "joint_band_covers_truth_on_audit_grid",
  "joint_continuum_band_covers_truth",
  "joint_any_primary_false_report",
  "joint_pointwise_grid_curve_coverage",
  "joint_derivative_assisted_grid_coverage",
  "joint_derivative_assisted_false_exact_root_report",
  "across_interval_joint_band_covers_truth_on_grid",
  "across_interval_joint_band_covers_truth_on_audit_grid",
  "across_interval_joint_continuum_band_covers_truth",
  "across_interval_joint_any_primary_false_report"
)
binary_metrics <- if (nrow(results) > 0L) {
  intersect(binary_candidates, names(results))
} else {
  binary_candidates
}

metric_rows <- list()
metric_counter <- 0L
for (group_key in names(groups)) {
  group <- groups[[group_key]]
  metadata <- group_metadata[[group_key]][, group_columns, drop = FALSE]
  for (metric in binary_metrics) {
    score <- binary_summary(group[[metric]])
    metric_counter <- metric_counter + 1L
    metric_rows[[metric_counter]] <- cbind(
      metadata,
      data.frame(
        metric = metric,
        successes = as.integer(score[["successes"]]),
        denominator = as.integer(score[["denominator"]]),
        estimate = score[["estimate"]], mcse = score[["mcse"]],
        wilson_lower = score[["wilson_lower"]],
        wilson_upper = score[["wilson_upper"]],
        stringsAsFactors = FALSE
      )
    )
  }
}
metric_table <- do.call(rbind, metric_rows)
row.names(metric_table) <- NULL

across_interval_metrics <- c(
  "across_interval_joint_band_covers_truth_on_grid",
  "across_interval_joint_band_covers_truth_on_audit_grid",
  "across_interval_joint_continuum_band_covers_truth",
  "across_interval_joint_any_primary_false_report"
)
across_group_columns <- setdiff(group_columns, "interval_id")
across_expected_groups <- unique(
  expected_groups[, across_group_columns, drop = FALSE]
)
across_interval_rows <- list()
across_metric_rows <- list()
across_interval_counter <- 0L
across_metric_counter <- 0L
for (index in seq_len(nrow(across_expected_groups))) {
  metadata <- across_expected_groups[index, , drop = FALSE]
  expected_match <- rep(TRUE, nrow(expected_groups))
  observed_match <- rep(TRUE, nrow(results))
  for (name in across_group_columns) {
    expected_match <- expected_match &
      as.character(expected_groups[[name]]) == as.character(metadata[[name]])
    if (nrow(results)) {
      observed_match <- observed_match &
        as.character(results[[name]]) == as.character(metadata[[name]])
    }
  }
  expected_intervals <- sort(unique(as.character(
    expected_groups$interval_id[expected_match]
  )))
  observed <- if (nrow(results)) {
    results[observed_match, , drop = FALSE]
  } else results
  across_deduplication <- deduplicate_across_interval_events(
    observed, expected_intervals
  )
  deduplicated <- across_deduplication$events
  incomplete_units <- across_deduplication$incomplete_units
  expected_replications <- as.integer(config$mc_reps_per_cell)
  successful_replications <- across_deduplication$successful_units
  complete_replications <- nrow(deduplicated)
  assert_true(
    successful_replications <= expected_replications &&
      complete_replications + incomplete_units == successful_replications,
    "The ALL_K replication accounting does not partition successful units."
  )
  scores <- lapply(across_interval_metrics, function(metric) {
    binary_summary(if (complete_replications) {
      deduplicated[[metric]]
    } else logical(0L))
  })
  names(scores) <- across_interval_metrics
  across_interval_counter <- across_interval_counter + 1L
  across_interval_rows[[across_interval_counter]] <- data.frame(
    study_id = config$study_id,
    cell_id = metadata$cell_id,
    cell_class = metadata$cell_class,
    variant_id = metadata$variant_id,
    scenario_id = metadata$scenario_id,
    n_x = metadata$n_x, n_y = metadata$n_y,
    moment_type = metadata$moment_type,
    interval_scope = "ALL_CONFIGURED_INTERVALS",
    configured_interval_ids = paste(expected_intervals, collapse = ";"),
    configured_interval_count = length(expected_intervals),
    expected_replications = expected_replications,
    successful_replications = successful_replications,
    complete_interval_replications = complete_replications,
    missing_component_replications = incomplete_units,
    failed_or_missing_unit_replications =
      expected_replications - complete_replications,
    across_interval_joint_construction_grid_band_coverage =
      scores[[1L]][["estimate"]],
    across_interval_joint_audit_grid_band_coverage =
      scores[[2L]][["estimate"]],
    across_interval_joint_continuum_band_coverage =
      scores[[3L]][["estimate"]],
    across_interval_joint_any_primary_false_report_rate =
      scores[[4L]][["estimate"]],
    stringsAsFactors = FALSE
  )
  for (metric in across_interval_metrics) {
    score <- scores[[metric]]
    across_metric_counter <- across_metric_counter + 1L
    across_metric_rows[[across_metric_counter]] <- data.frame(
      study_id = config$study_id,
      cell_id = metadata$cell_id,
      cell_class = metadata$cell_class,
      variant_id = metadata$variant_id,
      scenario_id = metadata$scenario_id,
      n_x = metadata$n_x, n_y = metadata$n_y,
      moment_type = metadata$moment_type,
      interval_scope = "ALL_CONFIGURED_INTERVALS",
      configured_interval_ids = paste(expected_intervals, collapse = ";"),
      metric = metric,
      successes = as.integer(score[["successes"]]),
      denominator = as.integer(score[["denominator"]]),
      estimate = score[["estimate"]], mcse = score[["mcse"]],
      wilson_lower = score[["wilson_lower"]],
      wilson_upper = score[["wilson_upper"]],
      stringsAsFactors = FALSE
    )
  }
}
across_interval_summary <- do.call(rbind, across_interval_rows)
row.names(across_interval_summary) <- NULL
across_interval_binary_table <- do.call(rbind, across_metric_rows)
row.names(across_interval_binary_table) <- NULL

continuous_candidates <- c(
  "bootstrap_grid_size", "bootstrap_grid_spacing",
  "enclosure_grid_size", "enclosure_grid_spacing",
  "configured_guard_maximum_levels",
  "configured_continuum_enclosure_maximum_levels",
  "guard_minimum_ratio", "guard_nodes_used", "guard_maximum_level",
  "critical_value", "maximum_half_width", "integrated_half_width",
  "continuum_truth_nodes_used", "continuum_truth_maximum_level",
  "continuum_truth_unresolved_cells", "enclosure_nodes_used",
  "enclosure_maximum_level", "enclosure_variance_unresolved_cells",
  "enclosure_unresolved_cells", "enclosure_limit_hit_cells",
  "enclosure_limit_hit_total_width",
  "enclosure_limit_hit_width_proportion",
  "enclosure_limit_hit_maximum_cell_width",
  "enclosure_depth_limit_cells", "enclosure_depth_limit_total_width",
  "enclosure_depth_limit_width_proportion",
  "enclosure_depth_limit_maximum_cell_width",
  "enclosure_node_limit_cells", "enclosure_node_limit_total_width",
  "enclosure_node_limit_width_proportion",
  "enclosure_node_limit_maximum_cell_width",
  "enclosure_variance_limit_cells", "enclosure_variance_limit_total_width",
  "enclosure_variance_limit_width_proportion",
  "enclosure_variance_limit_maximum_cell_width",
  "enclosure_statistical_limit_cells",
  "enclosure_statistical_limit_total_width",
  "enclosure_statistical_limit_width_proportion",
  "enclosure_statistical_limit_maximum_cell_width",
  "certified_alternation_count",
  "bracket_count", "bracket_length", "numerical_outer_component_count",
  "numerical_outer_total_length", "numerical_outer_retained_proportion",
  "numerical_outer_hausdorff", "tangency_local_outer_total_length",
  "tangency_scaled_outer_radius", "tangency_scaled_outer_hausdorff",
  "tangency_truncated_local_outer_total_length",
  "tangency_truncated_local_scaled_outer_radius",
  "tangency_truncated_local_scaled_outer_hausdorff",
  "tangency_global_outer_component_count_conditional_on_continuum_band_coverage",
  "tangency_global_outer_total_length_conditional_on_continuum_band_coverage",
  "tangency_global_outer_hausdorff_conditional_on_continuum_band_coverage",
  "tangency_global_outer_left_radius_conditional_on_continuum_band_coverage",
  "tangency_global_outer_right_radius_conditional_on_continuum_band_coverage",
  "tangency_global_scaled_outer_total_length_conditional_on_continuum_band_coverage",
  "tangency_global_scaled_outer_hausdorff_conditional_on_continuum_band_coverage",
  "tangency_global_scaled_outer_left_radius_conditional_on_continuum_band_coverage",
  "tangency_global_scaled_outer_right_radius_conditional_on_continuum_band_coverage",
  "tangency_global_root_component_left_conditional_on_continuum_band_coverage",
  "tangency_global_root_component_right_conditional_on_continuum_band_coverage",
  "tangency_global_root_component_length_conditional_on_continuum_band_coverage",
  "tangency_global_root_component_left_radius_conditional_on_continuum_band_coverage",
  "tangency_global_root_component_right_radius_conditional_on_continuum_band_coverage",
  "tangency_global_scaled_root_component_length_conditional_on_continuum_band_coverage",
  "tangency_global_scaled_root_component_left_radius_conditional_on_continuum_band_coverage",
  "tangency_global_scaled_root_component_right_radius_conditional_on_continuum_band_coverage",
  "empirical_root_count", "root_isolation_maximum_bisections",
  "root_error_conditional_on_isolation",
  "derivative_hat", "variance_at_root_hat", "root_standard_error",
  "wald_length", "wald_collection_bias_conditional_on_isolation",
  "wald_collection_rmse_conditional_on_isolation",
  "wald_collection_total_length", "root_multiplier_critical_value",
  "root_multiplier_total_length",
  "derivative_assisted_critical_value",
  "derivative_assisted_critical_value_ratio",
  "derivative_assisted_value_maximum_half_width",
  "derivative_assisted_derivative_maximum_half_width",
  "derivative_assisted_value_width_inflation",
  "derivative_assisted_guard_minimum_ratio",
  "derivative_assisted_guard_minimum_ratio_d0",
  "derivative_assisted_guard_minimum_ratio_d1",
  "derivative_assisted_guard_nodes_used",
  "derivative_assisted_guard_maximum_level",
  "derivative_assisted_value_enclosure_nodes_used",
  "derivative_assisted_derivative_enclosure_nodes_used",
  "derivative_assisted_value_unresolved_cells",
  "derivative_assisted_derivative_unresolved_cells",
  "derivative_assisted_value_variance_unresolved_cells",
  "derivative_assisted_derivative_variance_unresolved_cells",
  "derivative_assisted_certificate_count",
  "derivative_assisted_covariance_01_min",
  "derivative_assisted_covariance_01_max",
  "derivative_assisted_correlation_01_min",
  "derivative_assisted_correlation_01_max",
  "x_ess_ratio_p1", "y_ess_ratio_p1",
  "x_max_share_p1", "y_max_share_p1", "x_top1pct_share_p1",
  "y_top1pct_share_p1", "x_ess_ratio_upper", "y_ess_ratio_upper",
  "x_max_share_upper", "y_max_share_upper", "x_top1pct_share_upper",
  "y_top1pct_share_upper", "x_count_above_q95", "y_count_above_q95",
  "x_count_above_q99", "y_count_above_q99", "x_latent_high_count",
  "y_latent_high_count", "elapsed_seconds", "peak_r_heap_mb"
)
continuous_metrics <- if (nrow(results) > 0L) {
  intersect(continuous_candidates, names(results))
} else {
  continuous_candidates
}
continuous_rows <- list()
continuous_counter <- 0L
for (group_key in names(groups)) {
  group <- groups[[group_key]]
  metadata <- group_metadata[[group_key]][, group_columns, drop = FALSE]
  for (metric in continuous_metrics) {
    score <- finite_summary(group[[metric]])
    continuous_counter <- continuous_counter + 1L
    continuous_rows[[continuous_counter]] <- cbind(
      metadata,
      data.frame(
        metric = metric,
        rows = as.integer(score[["rows"]]),
        finite = as.integer(score[["finite"]]),
        missing = as.integer(score[["missing"]]),
        positive_infinite = as.integer(score[["positive_infinite"]]),
        negative_infinite = as.integer(score[["negative_infinite"]]),
        mean = score[["mean"]], sd = score[["sd"]],
        mean_mcse = score[["mean_mcse"]],
        q05 = score[["q05"]], q50 = score[["q50"]], q95 = score[["q95"]],
        stringsAsFactors = FALSE
      )
    )
  }
}
continuous_table <- do.call(rbind, continuous_rows)
row.names(continuous_table) <- NULL

metric_estimate <- function(group, metric, subset = rep(TRUE, nrow(group))) {
  if (!metric %in% names(group)) return(NA_real_)
  binary_summary(group[[metric]][subset])[["estimate"]]
}

metric_bound <- function(group, metric, side) {
  if (!metric %in% names(group)) return(NA_real_)
  binary_summary(group[[metric]])[[paste0("wilson_", side)]]
}

wide_rows <- lapply(names(groups), function(group_key) {
  group <- groups[[group_key]]
  metadata <- group_metadata[[group_key]]
  expected <- as.integer(config$mc_reps_per_cell)
  successful <- length(unique(group$unit_id))
  root_applicable <-
    !is.na(column_or(group, "root_isolation_applicable", NA)) &
    column_or(group, "root_isolation_applicable", FALSE)
  zero_budget_applicable <-
    !is.na(column_or(group, "zero_budget_diagnostic_applicable", NA)) &
    column_or(group, "zero_budget_diagnostic_applicable", FALSE)
  wald_reported <- !is.na(column_or(group, "wald_reported", NA)) &
    column_or(group, "wald_reported", FALSE)
  root_multiplier_reported <-
    !is.na(column_or(group, "root_multiplier_reported", NA)) &
    column_or(group, "root_multiplier_reported", FALSE)
  derivative_assisted_applicable <-
    !is.na(column_or(group, "derivative_assisted_applicable", NA)) &
    column_or(group, "derivative_assisted_applicable", FALSE)
  tangency_applicable <-
    !is.na(column_or(group, "tangency_localization_applicable", NA)) &
    column_or(group, "tangency_localization_applicable", FALSE)
  tangency_local_outer_length <- column_or(
    group, "tangency_local_outer_total_length", NA_real_
  )
  tangency_local_outer_length[!tangency_applicable] <- NA_real_
  tangency_neighborhood_width <- as.numeric(
    metadata$configured_tangency_neighborhood_width[1L]
  )
  tangency_local_outer_length_normalized <-
    tangency_local_outer_length / tangency_neighborhood_width
  tangency_local_outer_length_summary <- finite_summary(
    tangency_local_outer_length
  )
  tangency_local_outer_length_normalized_summary <- finite_summary(
    tangency_local_outer_length_normalized
  )
  tangency_scaled_outer_radius_summary <- finite_summary(column_or(
    group, "tangency_scaled_outer_radius", NA_real_
  ))
  tangency_scaled_outer_hausdorff_summary <- finite_summary(column_or(
    group, "tangency_scaled_outer_hausdorff", NA_real_
  ))
  tangency_global_masks <- tangency_geometry_condition_masks(group)
  tangency_global_summary <- function(metric) {
    values <- column_or(group, metric, NA_real_)
    finite_summary(values[tangency_global_masks$condition])
  }
  tangency_global_component_count_summary <- tangency_global_summary(
    "tangency_global_outer_component_count_conditional_on_continuum_band_coverage"
  )
  tangency_global_total_length_summary <- tangency_global_summary(
    "tangency_global_outer_total_length_conditional_on_continuum_band_coverage"
  )
  tangency_global_hausdorff_summary <- tangency_global_summary(
    "tangency_global_outer_hausdorff_conditional_on_continuum_band_coverage"
  )
  tangency_global_root_component_present <- as.logical(column_or(
    group,
    "tangency_global_root_component_present_conditional_on_continuum_band_coverage",
    NA
  ))
  assert_true(
    !anyNA(tangency_global_root_component_present[
      tangency_global_masks$condition
    ]),
    paste(
      "Root-component presence is missing inside the tangency",
      "continuum-coverage denominator."
    )
  )
  tangency_global_root_component_condition <-
    tangency_global_masks$condition &
      !is.na(tangency_global_root_component_present) &
      tangency_global_root_component_present
  tangency_global_root_component_present_score <- binary_summary(
    tangency_global_root_component_present[tangency_global_masks$condition]
  )
  tangency_global_root_component_summary <- function(metric) {
    values <- column_or(group, metric, NA_real_)
    finite_summary(values[tangency_global_root_component_condition])
  }
  tangency_global_left_radius_summary <- tangency_global_summary(
    "tangency_global_outer_left_radius_conditional_on_continuum_band_coverage"
  )
  tangency_global_right_radius_summary <- tangency_global_summary(
    "tangency_global_outer_right_radius_conditional_on_continuum_band_coverage"
  )
  tangency_global_scaled_total_length_summary <- tangency_global_summary(
    "tangency_global_scaled_outer_total_length_conditional_on_continuum_band_coverage"
  )
  tangency_global_scaled_hausdorff_summary <- tangency_global_summary(
    "tangency_global_scaled_outer_hausdorff_conditional_on_continuum_band_coverage"
  )
  tangency_global_scaled_left_radius_summary <- tangency_global_summary(
    "tangency_global_scaled_outer_left_radius_conditional_on_continuum_band_coverage"
  )
  tangency_global_scaled_right_radius_summary <- tangency_global_summary(
    "tangency_global_scaled_outer_right_radius_conditional_on_continuum_band_coverage"
  )
  tangency_root_component_length_summary <- tangency_global_root_component_summary(
    "tangency_global_root_component_length_conditional_on_continuum_band_coverage"
  )
  tangency_root_component_left_summary <- tangency_global_root_component_summary(
    "tangency_global_root_component_left_conditional_on_continuum_band_coverage"
  )
  tangency_root_component_right_summary <- tangency_global_root_component_summary(
    "tangency_global_root_component_right_conditional_on_continuum_band_coverage"
  )
  tangency_root_component_left_radius_summary <- tangency_global_root_component_summary(
    "tangency_global_root_component_left_radius_conditional_on_continuum_band_coverage"
  )
  tangency_root_component_right_radius_summary <- tangency_global_root_component_summary(
    "tangency_global_root_component_right_radius_conditional_on_continuum_band_coverage"
  )
  tangency_scaled_root_component_length_summary <- tangency_global_root_component_summary(
    "tangency_global_scaled_root_component_length_conditional_on_continuum_band_coverage"
  )
  tangency_scaled_root_component_left_radius_summary <- tangency_global_root_component_summary(
    "tangency_global_scaled_root_component_left_radius_conditional_on_continuum_band_coverage"
  )
  tangency_scaled_root_component_right_radius_summary <- tangency_global_root_component_summary(
    "tangency_global_scaled_root_component_right_radius_conditional_on_continuum_band_coverage"
  )
  root_errors <- column_or(
    group, "root_error_conditional_on_isolation", NA_real_
  )
  root_collection <- root_collection_errors_for_group(group)
  wald_report_score <- binary_summary(column_or(group, "wald_reported", NA))
  wald_conditional_score <- binary_summary(column_or(
    group[wald_reported, , drop = FALSE], "wald_conditional_cover", NA
  ))
  wald_report_and_cover_score <- binary_summary(column_or(
    group, "wald_report_and_cover", NA
  ))
  root_multiplier_report_score <- binary_summary(column_or(
    group, "root_multiplier_reported", NA
  ))
  root_multiplier_conditional_score <- binary_summary(column_or(
    group[root_multiplier_reported, , drop = FALSE],
    "root_multiplier_conditional_cover", NA
  ))
  root_multiplier_report_and_cover_score <- binary_summary(column_or(
    group, "root_multiplier_report_and_cover", NA
  ))
  data.frame(
    study_id = config$study_id, cell_id = metadata$cell_id[1L],
    cell_class = metadata$cell_class[1L],
    variant_id = metadata$variant_id[1L],
    scenario_id = metadata$scenario_id[1L],
    n_x = metadata$n_x[1L], n_y = metadata$n_y[1L],
    moment_type = metadata$moment_type[1L],
    interval_id = metadata$interval_id[1L],
    p_min = metadata$p_min[1L], p_max = metadata$p_max[1L],
    bootstrap_reps = metadata$configured_bootstrap_reps[1L],
    bootstrap_grid_size = mean_finite_or_na(
      column_or(group, "bootstrap_grid_size", NA_real_)
    ),
    bootstrap_grid_spacing = mean_finite_or_na(
      column_or(group, "bootstrap_grid_spacing", NA_real_)
    ),
    enclosure_grid_size = mean_finite_or_na(
      column_or(group, "enclosure_grid_size", NA_real_)
    ),
    enclosure_grid_spacing = mean_finite_or_na(
      column_or(group, "enclosure_grid_spacing", NA_real_)
    ),
    expected_replications = expected, successful_replications = successful,
    replication_failure_rate = max(0, expected - successful) / expected,
    continuum_band_coverage = metric_estimate(group, "continuum_band_covers_truth"),
    continuum_band_coverage_wilson_lower = metric_bound(
      group, "continuum_band_covers_truth", "lower"
    ),
    continuum_band_coverage_wilson_upper = metric_bound(
      group, "continuum_band_covers_truth", "upper"
    ),
    construction_grid_band_coverage = metric_estimate(group, "band_covers_truth_on_grid"),
    audit_grid_band_coverage = metric_estimate(group, "band_covers_truth_on_audit_grid"),
    joint_continuum_band_coverage = metric_estimate(
      group, "joint_continuum_band_covers_truth"
    ),
    across_interval_joint_construction_grid_band_coverage = metric_estimate(
      group, "across_interval_joint_band_covers_truth_on_grid"
    ),
    across_interval_joint_audit_grid_band_coverage = metric_estimate(
      group, "across_interval_joint_band_covers_truth_on_audit_grid"
    ),
    across_interval_joint_continuum_band_coverage = metric_estimate(
      group, "across_interval_joint_continuum_band_covers_truth"
    ),
    any_primary_false_report_rate = metric_estimate(group, "any_primary_false_report"),
    familywise_primary_error_rate = metric_estimate(
      group, "any_primary_false_report"
    ),
    joint_any_primary_false_report_rate = metric_estimate(
      group, "joint_any_primary_false_report"
    ),
    across_interval_joint_any_primary_false_report_rate = metric_estimate(
      group, "across_interval_joint_any_primary_false_report"
    ),
    any_primary_false_wilson_lower = metric_bound(
      group, "any_primary_false_report", "lower"
    ),
    any_primary_false_wilson_upper = metric_bound(
      group, "any_primary_false_report", "upper"
    ),
    familywise_primary_error_wilson_lower = metric_bound(
      group, "any_primary_false_report", "lower"
    ),
    familywise_primary_error_wilson_upper = metric_bound(
      group, "any_primary_false_report", "upper"
    ),
    reversal_certification_rate = metric_estimate(group, "reversal_certified"),
    exact_count_certification_rate = metric_estimate(group, "exact_count_certified"),
    continuum_no_root_certification_rate = metric_estimate(
      group, "continuum_no_root_certified"
    ),
    tangency_localization_applicable_rate = metric_estimate(
      group, "tangency_localization_applicable"
    ),
    tangency_local_outer_nonempty_rate_conditional_on_applicability =
      metric_estimate(
        group, "tangency_local_outer_nonempty", tangency_applicable
      ),
    tangency_local_geometry_scope = "TRUNCATED_TO_FIXED_NEIGHBORHOOD",
    tangency_truncated_local_censoring_rate_conditional_on_applicability =
      metric_estimate(
        group, "tangency_truncated_local_censored_by_window",
        tangency_applicable
      ),
    tangency_truncated_local_left_boundary_touch_rate_conditional_on_applicability =
      metric_estimate(
        group, "tangency_truncated_local_touches_left_boundary",
        tangency_applicable
      ),
    tangency_truncated_local_right_boundary_touch_rate_conditional_on_applicability =
      metric_estimate(
        group, "tangency_truncated_local_touches_right_boundary",
        tangency_applicable
      ),
    tangency_neighborhood_width = if (any(tangency_applicable)) {
      tangency_neighborhood_width
    } else NA_real_,
    mean_tangency_local_outer_total_length =
      tangency_local_outer_length_summary[["mean"]],
    q05_tangency_local_outer_total_length =
      tangency_local_outer_length_summary[["q05"]],
    q50_tangency_local_outer_total_length =
      tangency_local_outer_length_summary[["q50"]],
    q95_tangency_local_outer_total_length =
      tangency_local_outer_length_summary[["q95"]],
    mean_tangency_local_outer_total_length_normalized_by_neighborhood_width =
      tangency_local_outer_length_normalized_summary[["mean"]],
    q05_tangency_local_outer_total_length_normalized_by_neighborhood_width =
      tangency_local_outer_length_normalized_summary[["q05"]],
    q50_tangency_local_outer_total_length_normalized_by_neighborhood_width =
      tangency_local_outer_length_normalized_summary[["q50"]],
    q95_tangency_local_outer_total_length_normalized_by_neighborhood_width =
      tangency_local_outer_length_normalized_summary[["q95"]],
    mean_tangency_scaled_outer_radius_conditional_on_continuum_band_coverage =
      tangency_scaled_outer_radius_summary[["mean"]],
    q50_tangency_scaled_outer_radius_conditional_on_continuum_band_coverage =
      tangency_scaled_outer_radius_summary[["q50"]],
    mean_tangency_scaled_outer_hausdorff_conditional_on_continuum_band_coverage =
      tangency_scaled_outer_hausdorff_summary[["mean"]],
    q50_tangency_scaled_outer_hausdorff_conditional_on_continuum_band_coverage =
      tangency_scaled_outer_hausdorff_summary[["q50"]],
    tangency_global_geometry_applicable_replications =
      tangency_global_masks$applicable_denominator,
    tangency_global_geometry_continuum_coverage_conditioned_replications =
      tangency_global_masks$coverage_conditioned_denominator,
    tangency_global_outer_total_length_conditioned_rows =
      tangency_global_total_length_summary[["rows"]],
    tangency_global_outer_total_length_finite_count =
      tangency_global_total_length_summary[["finite"]],
    tangency_global_outer_total_length_missing_count =
      tangency_global_total_length_summary[["missing"]],
    tangency_global_outer_total_length_positive_infinite_count =
      tangency_global_total_length_summary[["positive_infinite"]],
    tangency_global_outer_hausdorff_conditioned_rows =
      tangency_global_hausdorff_summary[["rows"]],
    tangency_global_outer_hausdorff_finite_count =
      tangency_global_hausdorff_summary[["finite"]],
    tangency_global_outer_hausdorff_missing_count =
      tangency_global_hausdorff_summary[["missing"]],
    tangency_global_outer_hausdorff_positive_infinite_count =
      tangency_global_hausdorff_summary[["positive_infinite"]],
    tangency_global_root_component_geometry_scope =
      "CONDITIONAL_ON_CONTINUUM_COVERAGE_AND_ROOT_COMPONENT_PRESENT",
    tangency_global_root_component_present_replications =
      tangency_global_root_component_present_score[["successes"]],
    tangency_global_root_component_presence_denominator =
      tangency_global_root_component_present_score[["denominator"]],
    tangency_global_root_component_absent_replications_under_continuum_coverage =
      tangency_global_root_component_present_score[["denominator"]] -
        tangency_global_root_component_present_score[["successes"]],
    tangency_global_root_component_length_continuum_coverage_conditioned_rows =
      tangency_global_masks$coverage_conditioned_denominator,
    tangency_global_root_component_geometry_conditioned_replications =
      sum(tangency_global_root_component_condition),
    tangency_global_root_component_length_root_component_present_conditioned_rows =
      tangency_root_component_length_summary[["rows"]],
    tangency_global_root_component_length_finite_count =
      tangency_root_component_length_summary[["finite"]],
    tangency_global_root_component_length_missing_count =
      tangency_root_component_length_summary[["missing"]],
    tangency_global_root_component_length_positive_infinite_count =
      tangency_root_component_length_summary[["positive_infinite"]],
    tangency_global_geometry_conditioning_rate = if (
      tangency_global_masks$applicable_denominator > 0L
    ) {
      tangency_global_masks$coverage_conditioned_denominator /
        tangency_global_masks$applicable_denominator
    } else NA_real_,
    tangency_global_outer_nonempty_rate_conditional_on_continuum_band_coverage =
      metric_estimate(
        group,
        "tangency_global_outer_nonempty_conditional_on_continuum_band_coverage",
        tangency_global_masks$condition
      ),
    tangency_global_outer_contains_root_rate_conditional_on_continuum_band_coverage =
      metric_estimate(
        group,
        "tangency_global_outer_contains_root_conditional_on_continuum_band_coverage",
        tangency_global_masks$condition
      ),
    tangency_global_root_component_present_rate_conditional_on_continuum_band_coverage =
      metric_estimate(
        group,
        "tangency_global_root_component_present_conditional_on_continuum_band_coverage",
        tangency_global_masks$condition
      ),
    mean_tangency_global_outer_component_count_conditional_on_continuum_band_coverage =
      tangency_global_component_count_summary[["mean"]],
    q50_tangency_global_outer_component_count_conditional_on_continuum_band_coverage =
      tangency_global_component_count_summary[["q50"]],
    mean_tangency_global_outer_total_length_conditional_on_continuum_band_coverage =
      tangency_global_total_length_summary[["mean"]],
    q05_tangency_global_outer_total_length_conditional_on_continuum_band_coverage =
      tangency_global_total_length_summary[["q05"]],
    q50_tangency_global_outer_total_length_conditional_on_continuum_band_coverage =
      tangency_global_total_length_summary[["q50"]],
    q95_tangency_global_outer_total_length_conditional_on_continuum_band_coverage =
      tangency_global_total_length_summary[["q95"]],
    mean_tangency_global_outer_hausdorff_conditional_on_continuum_band_coverage =
      tangency_global_hausdorff_summary[["mean"]],
    q05_tangency_global_outer_hausdorff_conditional_on_continuum_band_coverage =
      tangency_global_hausdorff_summary[["q05"]],
    q50_tangency_global_outer_hausdorff_conditional_on_continuum_band_coverage =
      tangency_global_hausdorff_summary[["q50"]],
    q95_tangency_global_outer_hausdorff_conditional_on_continuum_band_coverage =
      tangency_global_hausdorff_summary[["q95"]],
    mean_tangency_global_outer_left_radius_conditional_on_continuum_band_coverage =
      tangency_global_left_radius_summary[["mean"]],
    q50_tangency_global_outer_left_radius_conditional_on_continuum_band_coverage =
      tangency_global_left_radius_summary[["q50"]],
    mean_tangency_global_outer_right_radius_conditional_on_continuum_band_coverage =
      tangency_global_right_radius_summary[["mean"]],
    q50_tangency_global_outer_right_radius_conditional_on_continuum_band_coverage =
      tangency_global_right_radius_summary[["q50"]],
    mean_tangency_global_scaled_outer_total_length_conditional_on_continuum_band_coverage =
      tangency_global_scaled_total_length_summary[["mean"]],
    q50_tangency_global_scaled_outer_total_length_conditional_on_continuum_band_coverage =
      tangency_global_scaled_total_length_summary[["q50"]],
    mean_tangency_global_scaled_outer_hausdorff_conditional_on_continuum_band_coverage =
      tangency_global_scaled_hausdorff_summary[["mean"]],
    q50_tangency_global_scaled_outer_hausdorff_conditional_on_continuum_band_coverage =
      tangency_global_scaled_hausdorff_summary[["q50"]],
    mean_tangency_global_scaled_outer_left_radius_conditional_on_continuum_band_coverage =
      tangency_global_scaled_left_radius_summary[["mean"]],
    q50_tangency_global_scaled_outer_left_radius_conditional_on_continuum_band_coverage =
      tangency_global_scaled_left_radius_summary[["q50"]],
    mean_tangency_global_scaled_outer_right_radius_conditional_on_continuum_band_coverage =
      tangency_global_scaled_right_radius_summary[["mean"]],
    q50_tangency_global_scaled_outer_right_radius_conditional_on_continuum_band_coverage =
      tangency_global_scaled_right_radius_summary[["q50"]],
    mean_tangency_global_root_component_length_conditional_on_continuum_band_coverage =
      tangency_root_component_length_summary[["mean"]],
    q50_tangency_global_root_component_length_conditional_on_continuum_band_coverage =
      tangency_root_component_length_summary[["q50"]],
    q50_tangency_global_root_component_left_conditional_on_continuum_band_coverage =
      tangency_root_component_left_summary[["q50"]],
    q50_tangency_global_root_component_right_conditional_on_continuum_band_coverage =
      tangency_root_component_right_summary[["q50"]],
    mean_tangency_global_root_component_left_radius_conditional_on_continuum_band_coverage =
      tangency_root_component_left_radius_summary[["mean"]],
    q50_tangency_global_root_component_left_radius_conditional_on_continuum_band_coverage =
      tangency_root_component_left_radius_summary[["q50"]],
    mean_tangency_global_root_component_right_radius_conditional_on_continuum_band_coverage =
      tangency_root_component_right_radius_summary[["mean"]],
    q50_tangency_global_root_component_right_radius_conditional_on_continuum_band_coverage =
      tangency_root_component_right_radius_summary[["q50"]],
    mean_tangency_global_scaled_root_component_length_conditional_on_continuum_band_coverage =
      tangency_scaled_root_component_length_summary[["mean"]],
    q50_tangency_global_scaled_root_component_length_conditional_on_continuum_band_coverage =
      tangency_scaled_root_component_length_summary[["q50"]],
    mean_tangency_global_scaled_root_component_left_radius_conditional_on_continuum_band_coverage =
      tangency_scaled_root_component_left_radius_summary[["mean"]],
    q50_tangency_global_scaled_root_component_left_radius_conditional_on_continuum_band_coverage =
      tangency_scaled_root_component_left_radius_summary[["q50"]],
    mean_tangency_global_scaled_root_component_right_radius_conditional_on_continuum_band_coverage =
      tangency_scaled_root_component_right_radius_summary[["mean"]],
    q50_tangency_global_scaled_root_component_right_radius_conditional_on_continuum_band_coverage =
      tangency_scaled_root_component_right_radius_summary[["q50"]],
    mean_tangency_global_root_component_length_conditional_on_continuum_band_coverage_and_root_component_present =
      tangency_root_component_length_summary[["mean"]],
    q50_tangency_global_root_component_length_conditional_on_continuum_band_coverage_and_root_component_present =
      tangency_root_component_length_summary[["q50"]],
    mean_tangency_global_root_component_left_radius_conditional_on_continuum_band_coverage_and_root_component_present =
      tangency_root_component_left_radius_summary[["mean"]],
    q50_tangency_global_root_component_left_radius_conditional_on_continuum_band_coverage_and_root_component_present =
      tangency_root_component_left_radius_summary[["q50"]],
    mean_tangency_global_root_component_right_radius_conditional_on_continuum_band_coverage_and_root_component_present =
      tangency_root_component_right_radius_summary[["mean"]],
    q50_tangency_global_root_component_right_radius_conditional_on_continuum_band_coverage_and_root_component_present =
      tangency_root_component_right_radius_summary[["q50"]],
    mean_tangency_global_scaled_root_component_length_conditional_on_continuum_band_coverage_and_root_component_present =
      tangency_scaled_root_component_length_summary[["mean"]],
    q50_tangency_global_scaled_root_component_length_conditional_on_continuum_band_coverage_and_root_component_present =
      tangency_scaled_root_component_length_summary[["q50"]],
    q50_tangency_global_root_component_left_conditional_on_continuum_band_coverage_and_root_component_present =
      tangency_root_component_left_summary[["q50"]],
    q50_tangency_global_root_component_right_conditional_on_continuum_band_coverage_and_root_component_present =
      tangency_root_component_right_summary[["q50"]],
    mean_tangency_global_scaled_root_component_left_radius_conditional_on_continuum_band_coverage_and_root_component_present =
      tangency_scaled_root_component_left_radius_summary[["mean"]],
    q50_tangency_global_scaled_root_component_left_radius_conditional_on_continuum_band_coverage_and_root_component_present =
      tangency_scaled_root_component_left_radius_summary[["q50"]],
    mean_tangency_global_scaled_root_component_right_radius_conditional_on_continuum_band_coverage_and_root_component_present =
      tangency_scaled_root_component_right_radius_summary[["mean"]],
    q50_tangency_global_scaled_root_component_right_radius_conditional_on_continuum_band_coverage_and_root_component_present =
      tangency_scaled_root_component_right_radius_summary[["q50"]],
    tangency_global_outer_K_boundary_touch_rate_conditional_on_continuum_band_coverage =
      metric_estimate(
        group,
        "tangency_global_outer_touches_K_boundary_conditional_on_continuum_band_coverage",
        tangency_global_masks$condition
      ),
    tangency_global_outer_K_left_touch_rate_conditional_on_continuum_band_coverage =
      metric_estimate(
        group,
        "tangency_global_outer_touches_K_left_conditional_on_continuum_band_coverage",
        tangency_global_masks$condition
      ),
    tangency_global_outer_K_right_touch_rate_conditional_on_continuum_band_coverage =
      metric_estimate(
        group,
        "tangency_global_outer_touches_K_right_conditional_on_continuum_band_coverage",
        tangency_global_masks$condition
      ),
    tangency_global_root_component_K_boundary_touch_rate_conditional_on_continuum_band_coverage =
      metric_estimate(
        group,
        "tangency_global_root_component_touches_K_boundary_conditional_on_continuum_band_coverage",
        tangency_global_masks$condition
      ),
    zero_budget_rejection_rate = metric_estimate(
      group, "zero_budget_rejected", zero_budget_applicable
    ),
    outer_set_truth_coverage = metric_estimate(group, "outer_set_covers_truth"),
    root_isolation_report_rate = metric_estimate(
      group, "root_isolation_success", root_applicable
    ),
    root_bias_conditional_on_isolation = mean_finite_or_na(root_errors),
    root_rmse_conditional_on_isolation = rmse_finite_or_na(root_errors),
    root_collection_applicable_replications = root_collection$applicable,
    root_collection_joint_isolation_reports = root_collection$joint,
    root_collection_joint_isolation_report_rate = if (
      root_collection$applicable > 0L
    ) root_collection$joint / root_collection$applicable else NA_real_,
    root_collection_conditional_replication_denominator =
      root_collection$joint,
    root_collection_error_denominator = length(root_collection$errors),
    root_collection_bias_conditional_on_joint_isolation =
      mean_finite_or_na(root_collection$errors),
    root_collection_rmse_conditional_on_joint_isolation =
      rmse_finite_or_na(root_collection$errors),
    wald_report_rate = metric_estimate(group, "wald_reported"),
    wald_applicable_replications = wald_report_score[["denominator"]],
    wald_report_count = wald_report_score[["successes"]],
    wald_conditional_coverage = metric_estimate(
      group, "wald_conditional_cover", wald_reported
    ),
    wald_conditional_coverage_denominator =
      wald_conditional_score[["denominator"]],
    wald_conditional_coverage_count =
      wald_conditional_score[["successes"]],
    wald_vacuous_fallback_set_coverage = metric_estimate(
      group, "wald_formal_unconditional_cover"
    ),
    wald_report_and_cover_rate = metric_estimate(group, "wald_report_and_cover"),
    wald_report_and_cover_denominator =
      wald_report_and_cover_score[["denominator"]],
    wald_report_and_cover_count =
      wald_report_and_cover_score[["successes"]],
    root_multiplier_report_rate = metric_estimate(group, "root_multiplier_reported"),
    root_multiplier_applicable_replications =
      root_multiplier_report_score[["denominator"]],
    root_multiplier_report_count =
      root_multiplier_report_score[["successes"]],
    root_multiplier_conditional_coverage = metric_estimate(
      group, "root_multiplier_conditional_cover", root_multiplier_reported
    ),
    root_multiplier_conditional_coverage_denominator =
      root_multiplier_conditional_score[["denominator"]],
    root_multiplier_conditional_coverage_count =
      root_multiplier_conditional_score[["successes"]],
    root_multiplier_vacuous_fallback_set_coverage = metric_estimate(
      group, "root_multiplier_formal_unconditional_cover"
    ),
    root_multiplier_report_and_cover_rate = metric_estimate(
      group, "root_multiplier_report_and_cover"
    ),
    root_multiplier_report_and_cover_denominator =
      root_multiplier_report_and_cover_score[["denominator"]],
    root_multiplier_report_and_cover_count =
      root_multiplier_report_and_cover_score[["successes"]],
    derivative_assisted_band_report_rate = metric_estimate(
      group, "derivative_assisted_reported", derivative_assisted_applicable
    ),
    derivative_assisted_guard_failure_rate = metric_estimate(
      group, "derivative_assisted_guard_failed", derivative_assisted_applicable
    ),
    derivative_assisted_operational_failure_rate = metric_estimate(
      group, "derivative_assisted_failure", derivative_assisted_applicable
    ),
    derivative_assisted_local_exact_root_certification_rate = metric_estimate(
      group, "derivative_assisted_local_exact_root_certified",
      derivative_assisted_applicable
    ),
    derivative_assisted_exact_root_on_K_certification_rate = metric_estimate(
      group, "derivative_assisted_exact_root_on_K_certified",
      derivative_assisted_applicable
    ),
    derivative_assisted_false_exact_root_report_rate = metric_estimate(
      group, "derivative_assisted_false_exact_root_report",
      derivative_assisted_applicable
    ),
    derivative_assisted_false_exact_root_on_K_report_rate = metric_estimate(
      group, "derivative_assisted_false_exact_root_on_K_report",
      derivative_assisted_applicable
    ),
    derivative_assisted_joint_grid_coverage = metric_estimate(
      group, "derivative_assisted_joint_grid_coverage",
      derivative_assisted_applicable
    ),
    mean_derivative_assisted_critical_value = mean_finite_or_na(column_or(
      group, "derivative_assisted_critical_value", NA_real_
    )),
    mean_derivative_assisted_value_width_inflation = mean_finite_or_na(column_or(
      group, "derivative_assisted_value_width_inflation", NA_real_
    )),
    pointwise_grid_curve_coverage = metric_estimate(group, "pointwise_grid_curve_coverage"),
    pointwise_false_sign_report_rate = metric_estimate(group, "pointwise_false_sign_report"),
    guard_failure_rate = metric_estimate(group, "guard_failed"),
    enclosure_limit_rate = metric_estimate(group, "enclosure_limit_hit"),
    enclosure_depth_limit_rate = metric_estimate(
      group, "enclosure_depth_limit_hit"
    ),
    enclosure_node_limit_rate = metric_estimate(
      group, "enclosure_node_limit_hit"
    ),
    enclosure_variance_limit_rate = metric_estimate(
      group, "enclosure_variance_limit_hit"
    ),
    enclosure_statistical_limit_rate = metric_estimate(
      group, "enclosure_statistical_limit_hit"
    ),
    mean_enclosure_limit_hit_cells = mean_finite_or_na(column_or(
      group, "enclosure_limit_hit_cells", NA_real_
    )),
    mean_enclosure_limit_hit_total_width = mean_finite_or_na(column_or(
      group, "enclosure_limit_hit_total_width", NA_real_
    )),
    mean_enclosure_limit_hit_width_proportion = mean_finite_or_na(column_or(
      group, "enclosure_limit_hit_width_proportion", NA_real_
    )),
    maximum_enclosure_limit_hit_cell_width = {
      values <- column_or(
        group, "enclosure_limit_hit_maximum_cell_width", NA_real_
      )
      values <- values[is.finite(values)]
      if (length(values)) max(values) else NA_real_
    },
    mean_enclosure_depth_limit_cells = mean_finite_or_na(column_or(
      group, "enclosure_depth_limit_cells", NA_real_
    )),
    mean_enclosure_depth_limit_total_width = mean_finite_or_na(column_or(
      group, "enclosure_depth_limit_total_width", NA_real_
    )),
    mean_enclosure_node_limit_cells = mean_finite_or_na(column_or(
      group, "enclosure_node_limit_cells", NA_real_
    )),
    mean_enclosure_node_limit_total_width = mean_finite_or_na(column_or(
      group, "enclosure_node_limit_total_width", NA_real_
    )),
    mean_enclosure_variance_limit_cells = mean_finite_or_na(column_or(
      group, "enclosure_variance_limit_cells", NA_real_
    )),
    mean_enclosure_variance_limit_total_width = mean_finite_or_na(column_or(
      group, "enclosure_variance_limit_total_width", NA_real_
    )),
    mean_enclosure_statistical_limit_cells = mean_finite_or_na(column_or(
      group, "enclosure_statistical_limit_cells", NA_real_
    )),
    mean_enclosure_statistical_limit_total_width = mean_finite_or_na(column_or(
      group, "enclosure_statistical_limit_total_width", NA_real_
    )),
    truth_score_limit_rate = metric_estimate(group, "continuum_truth_limit_hit"),
    mean_critical_value = mean_finite_or_na(column_or(group, "critical_value", NA_real_)),
    mean_maximum_half_width = mean_finite_or_na(
      column_or(group, "maximum_half_width", NA_real_)
    ),
    mean_outer_set_length = mean_finite_or_na(
      column_or(group, "numerical_outer_total_length", NA_real_)
    ),
    mean_outer_set_retained_proportion = mean_finite_or_na(
      column_or(group, "numerical_outer_retained_proportion", NA_real_)
    ),
    mean_y_ess_ratio_upper = mean_finite_or_na(
      column_or(group, "y_ess_ratio_upper", NA_real_)
    ),
    mean_y_max_share_upper = mean_finite_or_na(
      column_or(group, "y_max_share_upper", NA_real_)
    ),
    mean_elapsed_seconds_per_replication = mean_finite_or_na(
      column_or(group, "elapsed_seconds", NA_real_)
    ),
    mean_peak_r_heap_mb_per_replication = mean_finite_or_na(
      column_or(group, "peak_r_heap_mb", NA_real_)
    ),
    stringsAsFactors = FALSE
  )
})
summary_table <- do.call(rbind, wide_rows)
row.names(summary_table) <- NULL
summary_table <- summary_table[order(
  summary_table$scenario_id, summary_table$n_x, summary_table$n_y,
  summary_table$variant_id, summary_table$moment_type, summary_table$p_max
), , drop = FALSE]

root_metric_template <- data.frame(
  study_id = character(0L), cell_id = character(0L),
  cell_class = character(0L), variant_id = character(0L),
  scenario_id = character(0L), n_x = integer(0L), n_y = integer(0L),
  moment_type = character(0L), interval_id = character(0L),
  root_id = character(0L), true_root = numeric(0L),
  applicable_successful_replications = integer(0L),
  joint_isolation_replications = integer(0L),
  finite_error_replications = integer(0L),
  joint_isolation_report_rate = numeric(0L),
  bias_conditional_on_joint_isolation = numeric(0L),
  rmse_conditional_on_joint_isolation = numeric(0L),
  rmse_mcse_conditional_on_joint_isolation = numeric(0L),
  error_sd_conditional_on_joint_isolation = numeric(0L),
  bias_mcse_conditional_on_joint_isolation = numeric(0L),
  error_q05_conditional_on_joint_isolation = numeric(0L),
  error_q50_conditional_on_joint_isolation = numeric(0L),
  error_q95_conditional_on_joint_isolation = numeric(0L),
  stringsAsFactors = FALSE
)
root_metric_rows <- list()
root_metric_counter <- 0L
for (group_key in names(groups)) {
  group <- groups[[group_key]]
  metadata <- group_metadata[[group_key]]
  root_details <- root_collection_details_for_group(group)
  if (root_details$applicable == 0L) next
  reference_ids <- root_details$root_ids
  reference_truth <- root_details$true_roots
  error_matrix <- root_details$error_matrix
  for (root_index in seq_along(reference_ids)) {
    errors <- if (nrow(error_matrix)) {
      error_matrix[, root_index]
    } else numeric(0L)
    assert_true(
      length(errors) == root_details$joint && all(is.finite(errors)),
      "A root-specific conditional denominator was silently reduced."
    )
    score <- finite_summary(errors)
    root_metric_counter <- root_metric_counter + 1L
    root_metric_rows[[root_metric_counter]] <- data.frame(
      study_id = config$study_id,
      cell_id = metadata$cell_id[1L],
      cell_class = metadata$cell_class[1L],
      variant_id = metadata$variant_id[1L],
      scenario_id = metadata$scenario_id[1L],
      n_x = metadata$n_x[1L], n_y = metadata$n_y[1L],
      moment_type = metadata$moment_type[1L],
      interval_id = metadata$interval_id[1L],
      root_id = reference_ids[root_index],
      true_root = reference_truth[root_index],
      applicable_successful_replications = root_details$applicable,
      joint_isolation_replications = root_details$joint,
      finite_error_replications = length(errors),
      joint_isolation_report_rate = root_details$joint / root_details$applicable,
      bias_conditional_on_joint_isolation = mean_finite_or_na(errors),
      rmse_conditional_on_joint_isolation = rmse_finite_or_na(errors),
      rmse_mcse_conditional_on_joint_isolation = if (
          length(errors) >= 2L && rmse_finite_or_na(errors) > 0) {
        stats::sd(errors^2) /
          (2 * rmse_finite_or_na(errors) * sqrt(length(errors)))
      } else if (length(errors) >= 2L && all(errors == 0)) {
        0
      } else NA_real_,
      error_sd_conditional_on_joint_isolation = score[["sd"]],
      bias_mcse_conditional_on_joint_isolation = if (length(errors) >= 2L) {
        stats::sd(errors) / sqrt(length(errors))
      } else NA_real_,
      error_q05_conditional_on_joint_isolation = score[["q05"]],
      error_q50_conditional_on_joint_isolation = score[["q50"]],
      error_q95_conditional_on_joint_isolation = score[["q95"]],
      stringsAsFactors = FALSE
    )
  }
}
root_metric_table <- if (length(root_metric_rows)) {
  do.call(rbind, root_metric_rows)
} else root_metric_template
root_metric_table <- root_metric_table[, names(root_metric_template), drop = FALSE]

# Root inference is computed on a prespecified isolating interval and can be
# repeated verbatim on several K rows that contain that interval.  Publication
# output must not count those copies as independent results.  Collapse them
# here, after asserting equality of every inferential quantity that is meant to
# be invariant to K.  K-specific band and outer-set quantities remain in the
# ordinary summary tables and are deliberately excluded from this collapse.
values_identical_for_deduplication <- function(values, tolerance = 1e-12) {
  if (is.factor(values)) values <- as.character(values)
  if (is.numeric(values) || is.integer(values)) {
    values <- as.numeric(values)
    if (length(values) <= 1L) return(TRUE)
    if (!all(is.na(values) == is.na(values[1L]))) return(FALSE)
    observed <- values[!is.na(values)]
    if (length(observed) <= 1L) return(TRUE)
    if (!all(is.infinite(observed) == is.infinite(observed[1L]))) return(FALSE)
    finite <- observed[is.finite(observed)]
    if (!length(finite)) return(length(unique(observed)) == 1L)
    max(abs(finite - finite[1L])) <= tolerance * max(1, abs(finite))
  } else {
    values <- as.character(values)
    if (length(values) <= 1L) return(TRUE)
    all(is.na(values) == is.na(values[1L])) &&
      length(unique(values[!is.na(values)])) <= 1L
  }
}

root_collection_identity_fields <- c(
  "study_id", "cell_id", "cell_class", "variant_id", "scenario_id",
  "n_x", "n_y", "moment_type"
)
root_collection_inference_fields <- c(
  "root_collection_applicable_replications",
  "root_collection_joint_isolation_reports",
  "root_collection_joint_isolation_report_rate",
  "root_collection_conditional_replication_denominator",
  "root_collection_error_denominator",
  "root_collection_bias_conditional_on_joint_isolation",
  "root_collection_rmse_conditional_on_joint_isolation",
  "wald_report_rate", "wald_applicable_replications", "wald_report_count",
  "wald_conditional_coverage", "wald_conditional_coverage_denominator",
  "wald_conditional_coverage_count", "wald_report_and_cover_rate",
  "wald_report_and_cover_denominator", "wald_report_and_cover_count",
  "root_multiplier_report_rate", "root_multiplier_applicable_replications",
  "root_multiplier_report_count", "root_multiplier_conditional_coverage",
  "root_multiplier_conditional_coverage_denominator",
  "root_multiplier_conditional_coverage_count",
  "root_multiplier_report_and_cover_rate",
  "root_multiplier_report_and_cover_denominator",
  "root_multiplier_report_and_cover_count"
)
assert_true(
  all(c(root_collection_identity_fields,
        root_collection_inference_fields) %in% names(summary_table)),
  "The root-collection publication fields are incomplete."
)
root_collection_source <- summary_table[
  is.finite(summary_table$root_collection_applicable_replications) &
    summary_table$root_collection_applicable_replications > 0L,
  , drop = FALSE
]
root_collection_groups <- if (nrow(root_collection_source)) split(
  root_collection_source,
  interaction(
    root_collection_source$cell_id,
    root_collection_source$variant_id,
    root_collection_source$moment_type,
    drop = TRUE, sep = "::"
  )
) else list()
root_collection_deduplicated_rows <- lapply(
  root_collection_groups,
  function(group) {
    for (field in c(root_collection_identity_fields,
                    root_collection_inference_fields)) {
      assert_true(
        values_identical_for_deduplication(group[[field]]),
        sprintf(
          "Root-collection field %s disagrees across K for cell %s.",
          field, group$cell_id[1L]
        )
      )
    }
    data.frame(
      group[1L, c(root_collection_identity_fields,
                  root_collection_inference_fields), drop = FALSE],
      source_interval_ids = paste(
        sort(unique(as.character(group$interval_id))), collapse = ";"
      ),
      source_interval_count = length(unique(as.character(group$interval_id))),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }
)
root_collection_summary_deduplicated <- if (
    length(root_collection_deduplicated_rows)) {
  do.call(rbind, root_collection_deduplicated_rows)
} else {
  template <- summary_table[FALSE, c(
    root_collection_identity_fields, root_collection_inference_fields
  ), drop = FALSE]
  template$source_interval_ids <- character(0L)
  template$source_interval_count <- integer(0L)
  template
}
row.names(root_collection_summary_deduplicated) <- NULL

root_metric_identity_fields <- setdiff(
  names(root_metric_template), "interval_id"
)
root_metric_groups <- if (nrow(root_metric_table)) split(
  root_metric_table,
  interaction(
    root_metric_table$cell_id, root_metric_table$variant_id,
    root_metric_table$moment_type, root_metric_table$root_id,
    drop = TRUE, sep = "::"
  )
) else list()
root_metric_deduplicated_rows <- lapply(root_metric_groups, function(group) {
  for (field in root_metric_identity_fields) {
    assert_true(
      values_identical_for_deduplication(group[[field]]),
      sprintf(
        "Root-specific field %s disagrees across K for cell %s/root %s.",
        field, group$cell_id[1L], group$root_id[1L]
      )
    )
  }
  data.frame(
    group[1L, root_metric_identity_fields, drop = FALSE],
    source_interval_ids = paste(
      sort(unique(as.character(group$interval_id))), collapse = ";"
    ),
    source_interval_count = length(unique(as.character(group$interval_id))),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
})
root_metrics_deduplicated_by_root <- if (length(root_metric_deduplicated_rows)) {
  do.call(rbind, root_metric_deduplicated_rows)
} else {
  template <- root_metric_table[FALSE, root_metric_identity_fields, drop = FALSE]
  template$source_interval_ids <- character(0L)
  template$source_interval_count <- integer(0L)
  template
}
row.names(root_metrics_deduplicated_by_root) <- NULL

tail_rows <- list()
tail_counter <- 0L
tail_binary_metric_specs <- c(
  continuum_band_coverage = "continuum_band_covers_truth",
  any_primary_false_report_rate = "any_primary_false_report",
  reversal_certification_rate = "reversal_certified",
  exact_count_certification_rate = "exact_count_certified",
  derivative_assisted_local_exact_root_certification_rate =
    "derivative_assisted_local_exact_root_certified",
  derivative_assisted_false_exact_root_report_rate =
    "derivative_assisted_false_exact_root_report"
)
if ("y_latent_high_count" %in% names(results)) {
  latent <- results[is.finite(results$y_latent_high_count), , drop = FALSE]
  if (nrow(latent)) {
    latent$tail_count_stratum <- cut(
      latent$y_latent_high_count, breaks = c(-Inf, 0, 1, 2, 5, Inf),
      labels = c("0", "1", "2", "3-5", "6+"), right = TRUE
    )
    strata <- split(
      latent,
      interaction(
        latent$cell_id, latent$variant_id, latent$moment_type,
        latent$interval_id, latent$tail_count_stratum, drop = TRUE, sep = "::"
      )
    )
    for (group in strata) {
      tail_counter <- tail_counter + 1L
      tail_binary_scores <- lapply(tail_binary_metric_specs, function(metric) {
        if (metric %in% names(group)) {
          binary_summary(group[[metric]])
        } else {
          binary_summary(rep(NA, nrow(group)))
        }
      })
      tail_row <- data.frame(
        study_id = config$study_id, cell_id = group$cell_id[1L],
        cell_class = group$cell_class[1L],
        variant_id = group$variant_id[1L], scenario_id = group$scenario_id[1L],
        n_x = group$n_x[1L], n_y = group$n_y[1L],
        moment_type = group$moment_type[1L], interval_id = group$interval_id[1L],
        tail_count_stratum = as.character(group$tail_count_stratum[1L]),
        replications = nrow(group),
        continuum_band_coverage =
          tail_binary_scores$continuum_band_coverage[["estimate"]],
        any_primary_false_report_rate =
          tail_binary_scores$any_primary_false_report_rate[["estimate"]],
        reversal_certification_rate =
          tail_binary_scores$reversal_certification_rate[["estimate"]],
        exact_count_certification_rate =
          tail_binary_scores$exact_count_certification_rate[["estimate"]],
        derivative_assisted_local_exact_root_certification_rate =
          tail_binary_scores[[
            "derivative_assisted_local_exact_root_certification_rate"
          ]][["estimate"]],
        derivative_assisted_false_exact_root_report_rate =
          tail_binary_scores[[
            "derivative_assisted_false_exact_root_report_rate"
          ]][["estimate"]],
        mean_y_ess_ratio_upper = mean_finite_or_na(group$y_ess_ratio_upper),
        mean_y_max_share_upper = mean_finite_or_na(group$y_max_share_upper),
        stringsAsFactors = FALSE
      )
      for (output_name in names(tail_binary_metric_specs)) {
        score <- tail_binary_scores[[output_name]]
        tail_row[[paste0(output_name, "_denominator")]] <-
          as.integer(score[["denominator"]])
        tail_row[[paste0(output_name, "_mcse")]] <- score[["mcse"]]
        tail_row[[paste0(output_name, "_wilson_lower")]] <-
          score[["wilson_lower"]]
        tail_row[[paste0(output_name, "_wilson_upper")]] <-
          score[["wilson_upper"]]
      }
      tail_rows[[tail_counter]] <- tail_row
    }
  }
}
tail_table <- if (length(tail_rows)) do.call(rbind, tail_rows) else data.frame()

# Studies may use common sample and multiplier streams within a pairing_key.
# Report declared changes as within-replication contrasts; marginal summaries
# alone would discard that deliberate coupling.  The original sensitivity
# contrasts remain the backward-compatible default, while a targeted study
# can declare its own `paired_contrasts` table in the frozen configuration.
paired_table_template <- data.frame(
  study_id = character(0L), contrast_id = character(0L),
  contrast_definition = character(0L), pairing_key = character(0L),
  cell_class = character(0L), scenario_id = character(0L),
  n_x = integer(0L), n_y = integer(0L),
  moment_type = character(0L), interval_id = character(0L),
  p_min = numeric(0L), p_max = numeric(0L),
  reference_variant_id = character(0L),
  comparison_variant_id = character(0L),
  reference_cell_id = character(0L), comparison_cell_id = character(0L),
  metric_type = character(0L), metric = character(0L),
  expected_replications = integer(0L),
  reference_successful_replications = integer(0L),
  comparison_successful_replications = integer(0L),
  paired_successful_replications = integer(0L),
  metric_pair_denominator = integer(0L),
  reference_event_count = integer(0L), comparison_event_count = integer(0L),
  reference_only_event_count = integer(0L),
  comparison_only_event_count = integer(0L),
  disagreement_count = integer(0L), disagreement_rate = numeric(0L),
  disagreement_wilson_lower = numeric(0L),
  disagreement_wilson_upper = numeric(0L),
  reference_mean = numeric(0L), comparison_mean = numeric(0L),
  mean_difference = numeric(0L), difference_sd = numeric(0L),
  paired_mcse = numeric(0L), difference_q05 = numeric(0L),
  difference_q50 = numeric(0L), difference_q95 = numeric(0L),
  stringsAsFactors = FALSE
)

enclosure_depth_audit_template <- data.frame(
  study_id = character(0L), contrast_id = character(0L),
  pairing_key = character(0L), cell_class = character(0L),
  scenario_id = character(0L), n_x = integer(0L), n_y = integer(0L),
  moment_type = character(0L), interval_id = character(0L),
  reference_variant_id = character(0L),
  comparison_variant_id = character(0L),
  reference_cell_id = character(0L), comparison_cell_id = character(0L),
  reference_level = integer(0L), comparison_level = integer(0L),
  shallow_side = character(0L), deep_side = character(0L),
  shallow_level = integer(0L), deep_level = integer(0L),
  paired_replications = integer(0L),
  subset_checks = integer(0L), subset_failures = integer(0L),
  root_retention_checks = integer(0L),
  root_losses_under_continuum_coverage = integer(0L),
  length_bound_checks = integer(0L), length_bound_failures = integer(0L),
  maximum_length_bound_excess = numeric(0L),
  minimum_length_reduction = numeric(0L),
  decision_comparisons = integer(0L),
  decision_disagreements = integer(0L),
  decision_disagreement_replications = integer(0L),
  decision_disagreement_detail = character(0L),
  stringsAsFactors = FALSE
)

default_sensitivity_pair_specs <- data.frame(
  contrast_id = c("bootstrap_replication_effect", "grid_refinement_effect"),
  contrast_definition = c(
    "interaction minus fine_grid: B=3999 minus B=999 at the fine grid",
    "interaction minus finite_B: fine minus coarse grid at B=3999"
  ),
  reference_variant_id = c("fine_grid", "finite_B"),
  comparison_variant_id = c("interaction", "interaction"),
  stringsAsFactors = FALSE
)
configured_pair_specs <- config$paired_contrasts
if (is.null(configured_pair_specs)) {
  paired_specs <- default_sensitivity_pair_specs
  sensitivity_variant_ids <- names(config$sensitivity_variants %||% list())
  run_paired_comparisons <- all(
    unique(c(
      paired_specs$reference_variant_id,
      paired_specs$comparison_variant_id
    )) %in% sensitivity_variant_ids
  )
} else {
  assert_true(is.data.frame(configured_pair_specs),
              "paired_contrasts must be a data frame.")
  required_pair_columns <- c(
    "contrast_id", "contrast_definition", "reference_variant_id",
    "comparison_variant_id"
  )
  assert_true(all(required_pair_columns %in% names(configured_pair_specs)),
              "paired_contrasts is missing a required column.")
  paired_specs <- configured_pair_specs[, required_pair_columns, drop = FALSE]
  assert_true(nrow(paired_specs) > 0L &&
                !anyDuplicated(as.character(paired_specs$contrast_id)),
              "paired_contrasts must contain unique nonempty contrasts.")
  declared_variants <- unique(c(
    as.character(paired_specs$reference_variant_id),
    as.character(paired_specs$comparison_variant_id)
  ))
  available_variants <- unique(as.character(expected_groups$variant_id))
  assert_true(all(declared_variants %in% available_variants),
              "A declared paired contrast refers to an absent variant.")
  run_paired_comparisons <- TRUE
}

paired_sensitivity_table <- paired_table_template
enclosure_depth_audit_table <- enclosure_depth_audit_template
if (run_paired_comparisons) {
  variants <- config$sensitivity_variants
  legacy_variant_ids <- c("finite_B", "fine_grid", "interaction")
  if (!is.null(variants) && all(legacy_variant_ids %in% names(variants))) {
    assert_true(
      as.integer(variants$interaction$bootstrap_reps) ==
        as.integer(variants$finite_B$bootstrap_reps) &&
        as.integer(variants$interaction$grid_size) ==
        as.integer(variants$fine_grid$grid_size),
      paste(
        "The sensitivity variants do not isolate the advertised finite-B and",
        "grid-resolution contrasts."
      )
    )
  }

  pair_group_columns <- c(
    "pairing_key", "cell_class", "scenario_id", "n_x", "n_y",
    "moment_type", "interval_id"
  )
  expected_pair_groups <- list()
  expected_pair_counter <- 0L
  for (spec_index in seq_len(nrow(paired_specs))) {
    spec <- paired_specs[spec_index, , drop = FALSE]
    reference <- expected_groups[
      expected_groups$variant_id == spec$reference_variant_id,
      c(pair_group_columns, "cell_id", "p_min", "p_max"), drop = FALSE
    ]
    comparison <- expected_groups[
      expected_groups$variant_id == spec$comparison_variant_id,
      c(pair_group_columns, "cell_id", "p_min", "p_max"), drop = FALSE
    ]
    names(reference)[names(reference) == "cell_id"] <- "reference_cell_id"
    names(reference)[names(reference) == "p_min"] <- "reference_p_min"
    names(reference)[names(reference) == "p_max"] <- "reference_p_max"
    names(comparison)[names(comparison) == "cell_id"] <- "comparison_cell_id"
    names(comparison)[names(comparison) == "p_min"] <- "comparison_p_min"
    names(comparison)[names(comparison) == "p_max"] <- "comparison_p_max"
    joined <- merge(
      reference, comparison, by = pair_group_columns, all = FALSE, sort = FALSE
    )
    assert_true(
      nrow(joined) == nrow(reference) && nrow(joined) == nrow(comparison),
      sprintf(
        "Paired contrast %s does not pair every expected scientific group.",
        spec$contrast_id
      )
    )
    assert_true(
      all(joined$reference_p_min == joined$comparison_p_min) &&
        all(joined$reference_p_max == joined$comparison_p_max),
      sprintf("Paired contrast %s pairs unequal order intervals.",
              spec$contrast_id)
    )
    joined$contrast_id <- spec$contrast_id
    joined$contrast_definition <- spec$contrast_definition
    joined$reference_variant_id <- spec$reference_variant_id
    joined$comparison_variant_id <- spec$comparison_variant_id
    joined$p_min <- joined$reference_p_min
    joined$p_max <- joined$reference_p_max
    expected_pair_counter <- expected_pair_counter + 1L
    expected_pair_groups[[expected_pair_counter]] <- joined
  }
  expected_pair_groups <- do.call(rbind, expected_pair_groups)
  row.names(expected_pair_groups) <- NULL
  expected_pair_key <- do.call(paste, c(
    lapply(c("contrast_id", pair_group_columns), function(name) {
      as.character(expected_pair_groups[[name]])
    }),
    list(sep = "::")
  ))
  assert_true(!anyDuplicated(expected_pair_key),
              "The expected sensitivity-pair skeleton contains duplicates.")

  paired_metric_value <- function(paired, metric, side, default) {
    name <- paste0(metric, ".", side)
    if (name %in% names(paired)) paired[[name]] else rep(default, nrow(paired))
  }
  paired_rows <- list()
  paired_counter <- 0L
  enclosure_depth_audit_rows <- list()
  enclosure_depth_audit_counter <- 0L
  for (pair_index in seq_len(nrow(expected_pair_groups))) {
    pair <- expected_pair_groups[pair_index, , drop = FALSE]
    reference <- if (nrow(results) > 0L) {
      results[
        results$cell_id == pair$reference_cell_id &
          results$moment_type == pair$moment_type &
          results$interval_id == pair$interval_id,
        , drop = FALSE
      ]
    } else {
      data.frame(rep_id = integer(0L))
    }
    comparison <- if (nrow(results) > 0L) {
      results[
        results$cell_id == pair$comparison_cell_id &
          results$moment_type == pair$moment_type &
          results$interval_id == pair$interval_id,
        , drop = FALSE
      ]
    } else {
      data.frame(rep_id = integer(0L))
    }
    assert_true(!anyDuplicated(reference$rep_id) && !anyDuplicated(comparison$rep_id),
                "A paired variant has duplicate replication rows.")
    paired <- merge(
      reference, comparison, by = "rep_id", all = FALSE, sort = FALSE,
      suffixes = c(".reference", ".comparison")
    )
    if (nrow(paired) > 0L) {
      for (seed_name in c("seed_id", "bootstrap_seed_id")) {
        reference_seed <- paste0(seed_name, ".reference")
        comparison_seed <- paste0(seed_name, ".comparison")
        if (all(c(reference_seed, comparison_seed) %in% names(paired))) {
          assert_true(
            all(as.character(paired[[reference_seed]]) ==
                  as.character(paired[[comparison_seed]])),
            sprintf("Paired rows do not share %s.", seed_name)
          )
        }
      }
    }

    common <- data.frame(
      study_id = config$study_id,
      contrast_id = pair$contrast_id,
      contrast_definition = pair$contrast_definition,
      pairing_key = pair$pairing_key,
      cell_class = pair$cell_class,
      scenario_id = pair$scenario_id,
      n_x = pair$n_x, n_y = pair$n_y,
      moment_type = pair$moment_type, interval_id = pair$interval_id,
      p_min = pair$p_min, p_max = pair$p_max,
      reference_variant_id = pair$reference_variant_id,
      comparison_variant_id = pair$comparison_variant_id,
      reference_cell_id = pair$reference_cell_id,
      comparison_cell_id = pair$comparison_cell_id,
      expected_replications = as.integer(config$mc_reps_per_cell),
      reference_successful_replications = nrow(reference),
      comparison_successful_replications = nrow(comparison),
      paired_successful_replications = nrow(paired),
      stringsAsFactors = FALSE
    )

    reference_depth_name <- paste0(
      "configured_continuum_enclosure_maximum_levels.reference"
    )
    comparison_depth_name <- paste0(
      "configured_continuum_enclosure_maximum_levels.comparison"
    )
    if (nrow(paired) > 0L && all(c(
      reference_depth_name, comparison_depth_name
    ) %in% names(paired))) {
      reference_levels <- unique(suppressWarnings(as.integer(
        paired[[reference_depth_name]]
      )))
      comparison_levels <- unique(suppressWarnings(as.integer(
        paired[[comparison_depth_name]]
      )))
      reference_levels <- reference_levels[!is.na(reference_levels)]
      comparison_levels <- comparison_levels[!is.na(comparison_levels)]
      assert_true(
        length(reference_levels) == 1L &&
          length(comparison_levels) == 1L,
        "Paired variants do not persist one enclosure depth per group."
      )
      if (reference_levels[1L] != comparison_levels[1L]) {
        depth_audit <- paired_enclosure_depth_audit(
          paired, reference_levels[1L], comparison_levels[1L]
        )
        enclosure_depth_audit_counter <-
          enclosure_depth_audit_counter + 1L
        enclosure_depth_audit_rows[[enclosure_depth_audit_counter]] <-
          data.frame(
            study_id = config$study_id,
            contrast_id = pair$contrast_id,
            pairing_key = pair$pairing_key,
            cell_class = pair$cell_class,
            scenario_id = pair$scenario_id,
            n_x = pair$n_x, n_y = pair$n_y,
            moment_type = pair$moment_type,
            interval_id = pair$interval_id,
            reference_variant_id = pair$reference_variant_id,
            comparison_variant_id = pair$comparison_variant_id,
            reference_cell_id = pair$reference_cell_id,
            comparison_cell_id = pair$comparison_cell_id,
            reference_level = reference_levels[1L],
            comparison_level = comparison_levels[1L],
            shallow_side = depth_audit$shallow_side,
            deep_side = depth_audit$deep_side,
            shallow_level = depth_audit$shallow_level,
            deep_level = depth_audit$deep_level,
            paired_replications = depth_audit$paired_replications,
            subset_checks = depth_audit$subset_checks,
            subset_failures = depth_audit$subset_failures,
            root_retention_checks = depth_audit$root_retention_checks,
            root_losses_under_continuum_coverage =
              depth_audit$root_losses_under_continuum_coverage,
            length_bound_checks = depth_audit$length_bound_checks,
            length_bound_failures = depth_audit$length_bound_failures,
            maximum_length_bound_excess =
              depth_audit$maximum_length_bound_excess,
            minimum_length_reduction = depth_audit$minimum_length_reduction,
            decision_comparisons = depth_audit$decision_comparisons,
            decision_disagreements = depth_audit$decision_disagreements,
            decision_disagreement_replications =
              depth_audit$decision_disagreement_replications,
            decision_disagreement_detail =
              depth_audit$decision_disagreement_detail,
            stringsAsFactors = FALSE
          )
      }
    }

    for (metric in binary_metrics) {
      reference_values <- as.logical(paired_metric_value(
        paired, metric, "reference", NA
      ))
      comparison_values <- as.logical(paired_metric_value(
        paired, metric, "comparison", NA
      ))
      keep <- !is.na(reference_values) & !is.na(comparison_values)
      reference_values <- reference_values[keep]
      comparison_values <- comparison_values[keep]
      denominator <- length(reference_values)
      differences <- as.integer(comparison_values) - as.integer(reference_values)
      disagreements <- comparison_values != reference_values
      disagreement_count <- if (denominator) sum(disagreements) else 0L
      disagreement_interval <- wilson_interval(disagreement_count, denominator)
      paired_counter <- paired_counter + 1L
      paired_rows[[paired_counter]] <- cbind(
        common,
        data.frame(
          metric_type = "binary", metric = metric,
          metric_pair_denominator = denominator,
          reference_event_count = if (denominator) sum(reference_values) else 0L,
          comparison_event_count = if (denominator) sum(comparison_values) else 0L,
          reference_only_event_count = if (denominator) {
            sum(reference_values & !comparison_values)
          } else 0L,
          comparison_only_event_count = if (denominator) {
            sum(!reference_values & comparison_values)
          } else 0L,
          disagreement_count = disagreement_count,
          disagreement_rate = if (denominator) {
            disagreement_count / denominator
          } else NA_real_,
          disagreement_wilson_lower = disagreement_interval[["lower"]],
          disagreement_wilson_upper = disagreement_interval[["upper"]],
          reference_mean = if (denominator) mean(reference_values) else NA_real_,
          comparison_mean = if (denominator) mean(comparison_values) else NA_real_,
          mean_difference = if (denominator) mean(differences) else NA_real_,
          difference_sd = if (denominator >= 2L) stats::sd(differences) else NA_real_,
          paired_mcse = if (denominator >= 2L) {
            stats::sd(differences) / sqrt(denominator)
          } else NA_real_,
          difference_q05 = if (denominator) {
            unname(stats::quantile(differences, 0.05))
          } else NA_real_,
          difference_q50 = if (denominator) {
            unname(stats::quantile(differences, 0.50))
          } else NA_real_,
          difference_q95 = if (denominator) {
            unname(stats::quantile(differences, 0.95))
          } else NA_real_,
          stringsAsFactors = FALSE
        )
      )
    }

    for (metric in continuous_metrics) {
      reference_values <- suppressWarnings(as.numeric(paired_metric_value(
        paired, metric, "reference", NA_real_
      )))
      comparison_values <- suppressWarnings(as.numeric(paired_metric_value(
        paired, metric, "comparison", NA_real_
      )))
      keep <- is.finite(reference_values) & is.finite(comparison_values)
      reference_values <- reference_values[keep]
      comparison_values <- comparison_values[keep]
      denominator <- length(reference_values)
      differences <- comparison_values - reference_values
      paired_counter <- paired_counter + 1L
      paired_rows[[paired_counter]] <- cbind(
        common,
        data.frame(
          metric_type = "continuous", metric = metric,
          metric_pair_denominator = denominator,
          reference_event_count = NA_integer_, comparison_event_count = NA_integer_,
          reference_only_event_count = NA_integer_,
          comparison_only_event_count = NA_integer_,
          disagreement_count = NA_integer_, disagreement_rate = NA_real_,
          disagreement_wilson_lower = NA_real_,
          disagreement_wilson_upper = NA_real_,
          reference_mean = if (denominator) mean(reference_values) else NA_real_,
          comparison_mean = if (denominator) mean(comparison_values) else NA_real_,
          mean_difference = if (denominator) mean(differences) else NA_real_,
          difference_sd = if (denominator >= 2L) stats::sd(differences) else NA_real_,
          paired_mcse = if (denominator >= 2L) {
            stats::sd(differences) / sqrt(denominator)
          } else NA_real_,
          difference_q05 = if (denominator) {
            unname(stats::quantile(differences, 0.05))
          } else NA_real_,
          difference_q50 = if (denominator) {
            unname(stats::quantile(differences, 0.50))
          } else NA_real_,
          difference_q95 = if (denominator) {
            unname(stats::quantile(differences, 0.95))
          } else NA_real_,
          stringsAsFactors = FALSE
        )
      )
    }
  }
  paired_sensitivity_table <- do.call(rbind, paired_rows)
  row.names(paired_sensitivity_table) <- NULL
  paired_sensitivity_table <- paired_sensitivity_table[
    , names(paired_table_template), drop = FALSE
  ]
  paired_sensitivity_table <- paired_sensitivity_table[order(
    paired_sensitivity_table$scenario_id,
    paired_sensitivity_table$moment_type,
    paired_sensitivity_table$p_max,
    paired_sensitivity_table$contrast_id,
    paired_sensitivity_table$metric_type,
    paired_sensitivity_table$metric
  ), , drop = FALSE]
  if (length(enclosure_depth_audit_rows)) {
    enclosure_depth_audit_table <- do.call(
      rbind, enclosure_depth_audit_rows
    )
    row.names(enclosure_depth_audit_table) <- NULL
    enclosure_depth_audit_table <- enclosure_depth_audit_table[
      , names(enclosure_depth_audit_template), drop = FALSE
    ]
    enclosure_depth_audit_table <- enclosure_depth_audit_table[order(
      enclosure_depth_audit_table$scenario_id,
      enclosure_depth_audit_table$n_x,
      enclosure_depth_audit_table$n_y,
      enclosure_depth_audit_table$moment_type,
      enclosure_depth_audit_table$interval_id,
      enclosure_depth_audit_table$contrast_id
    ), , drop = FALSE]
  }
}

# Unit-level failure accounting is separate from scientific output rows,
# because one failed unit would otherwise be replicated across moment types
# and order intervals.  The frozen unit queue supplies the denominator.
queue <- readRDS(file.path(design_directory, "queue.rds"))
design_units <- queue$units
successful_unit_map <- if (nrow(results) > 0L) {
  unique(results[, c("cell_id", "unit_id"), drop = FALSE])
} else {
  data.frame(cell_id = character(0L), unit_id = integer(0L))
}
failed_unit_map <- if (nrow(numerical_failures) > 0L) {
  unique(numerical_failures[, c("cell_id", "unit_id"), drop = FALSE])
} else {
  data.frame(cell_id = character(0L), unit_id = integer(0L))
}
failure_rows <- lapply(seq_len(nrow(cell_metadata)), function(index) {
  metadata <- cell_metadata[index, , drop = FALSE]
  cell_id <- metadata$cell_id[1L]
  expected_units <- sum(design_units$cell_id == cell_id)
  successful_units <- sum(successful_unit_map$cell_id == cell_id)
  failed_units <- sum(failed_unit_map$cell_id == cell_id)
  assert_true(
    successful_units + failed_units == expected_units,
    sprintf("Successful and failed units do not partition cell %s.", cell_id)
  )
  data.frame(
    study_id = config$study_id,
    cell_id = cell_id,
    cell_class = metadata$cell_class[1L],
    variant_id = metadata$variant_id[1L],
    scenario_id = metadata$scenario_id[1L],
    n_x = metadata$n_x[1L], n_y = metadata$n_y[1L],
    expected_units = expected_units,
    successful_units = successful_units,
    numerical_failure_units = failed_units,
    numerical_failure_rate = failed_units / expected_units,
    stringsAsFactors = FALSE
  )
})
failure_summary <- do.call(rbind, failure_rows)
row.names(failure_summary) <- NULL

# Bridge.8-ready studies carry a frozen row contract.  Validate it before any
# summary is published so a missing metric family, accidental duplicated K
# event, or changed root denominator cannot yield a superficially complete
# output directory.
output_contract <- config$expected_output_contract %||% list()
if (length(output_contract)) {
  observed_contract <- list(
    cells = nrow(cell_metadata),
    units = validation$expected_units,
    tasks = validation$expected_tasks,
    independent_sample_streams = validation$independent_sample_streams,
    independent_bootstrap_streams = validation$independent_bootstrap_streams,
    scientific_groups = nrow(summary_table),
    replication_rows = nrow(results),
    binary_metric_rows = nrow(metric_table),
    continuous_metric_rows = nrow(continuous_table),
    paired_comparison_rows = nrow(paired_sensitivity_table),
    enclosure_depth_audit_rows = nrow(enclosure_depth_audit_table),
    across_interval_groups = nrow(across_interval_summary),
    across_interval_binary_metric_rows = nrow(across_interval_binary_table),
    root_metric_rows = nrow(root_metric_table),
    root_collection_deduplicated_rows =
      nrow(root_collection_summary_deduplicated),
    root_metric_deduplicated_rows =
      nrow(root_metrics_deduplicated_by_root)
  )
  checked_fields <- intersect(names(output_contract), names(observed_contract))
  assert_true(
    length(checked_fields) == length(output_contract),
    sprintf(
      "Study %s declares unknown expected-output fields: %s.",
      config$study_id,
      paste(setdiff(names(output_contract), names(observed_contract)),
            collapse = ", ")
    )
  )
  mismatches <- checked_fields[vapply(checked_fields, function(field) {
    !identical(
      as.integer(observed_contract[[field]]),
      as.integer(output_contract[[field]])
    )
  }, logical(1L))]
  assert_true(
    !length(mismatches),
    sprintf(
      "Study %s failed its frozen output contract: %s.",
      config$study_id,
      paste(vapply(mismatches, function(field) {
        sprintf(
          "%s observed=%s expected=%s", field,
          observed_contract[[field]], output_contract[[field]]
        )
      }, character(1L)), collapse = "; ")
    )
  )
}

summary_directory <- file.path(cmo_root(), "results", "summary", config$study_id)
dir.create(summary_directory, recursive = TRUE, showWarnings = FALSE)
atomic_write_csv(summary_table, file.path(summary_directory, "summary.csv"))
atomic_write_csv(metric_table, file.path(summary_directory, "binary_metrics.csv"))
atomic_write_csv(
  across_interval_summary,
  file.path(summary_directory, "across_interval_summary.csv")
)
atomic_write_csv(
  across_interval_binary_table,
  file.path(summary_directory, "across_interval_binary_metrics.csv")
)
atomic_write_csv(continuous_table, file.path(summary_directory, "continuous_metrics.csv"))
atomic_write_csv(
  root_metric_table,
  file.path(
    summary_directory, "root_metrics_conditional_on_joint_isolation.csv"
  )
)
atomic_write_csv(
  root_collection_summary_deduplicated,
  file.path(summary_directory, "root_collection_summary_deduplicated.csv")
)
atomic_write_csv(
  root_metrics_deduplicated_by_root,
  file.path(summary_directory, "root_metrics_deduplicated_by_root.csv")
)
atomic_write_csv(tail_table, file.path(summary_directory, "tail_stratified.csv"))
atomic_write_csv(
  paired_sensitivity_table,
  file.path(summary_directory, "sensitivity_paired_comparisons.csv")
)
atomic_write_csv(
  paired_sensitivity_table,
  file.path(summary_directory, "paired_comparisons.csv")
)
atomic_write_csv(
  enclosure_depth_audit_table,
  file.path(summary_directory, "enclosure_depth_audit.csv")
)
atomic_write_csv(numerical_failures, file.path(summary_directory, "numerical_failures.csv"))
atomic_write_csv(
  failure_summary,
  file.path(summary_directory, "numerical_failure_summary_by_cell.csv")
)
atomic_save_rds(results, file.path(summary_directory, "replication_results.rds"),
                compress = "xz")
atomic_save_rds(
  list(
    validation = validation, summary = summary_table,
    binary_metrics = metric_table, continuous_metrics = continuous_table,
    across_interval_summary = across_interval_summary,
    across_interval_binary_metrics = across_interval_binary_table,
    root_metrics_conditional_on_joint_isolation = root_metric_table,
    root_collection_summary_deduplicated =
      root_collection_summary_deduplicated,
    root_metrics_deduplicated_by_root = root_metrics_deduplicated_by_root,
    tail_stratified = tail_table,
    paired_comparisons = paired_sensitivity_table,
    sensitivity_paired_comparisons = paired_sensitivity_table,
    enclosure_depth_audit = enclosure_depth_audit_table,
    numerical_failures = numerical_failures,
    numerical_failure_summary = failure_summary,
    shard_paths = shard_paths,
    session = capture_session()
  ),
  file.path(summary_directory, "summary_bundle.rds"), compress = "xz"
)
atomic_write_lines(
  c(
    paste0("study_id=", config$study_id),
    paste0("schema_version=", config$schema_version),
    paste0("generated_utc=", timestamp_utc()),
    paste0("shards=", length(shard_paths)),
    paste0("replication_rows=", nrow(results)),
    paste0("scientific_groups=", nrow(summary_table)),
    paste0("across_interval_groups=", nrow(across_interval_summary)),
    paste0("across_interval_binary_metric_rows=",
           nrow(across_interval_binary_table)),
    paste0("root_metric_rows=", nrow(root_metric_table)),
    paste0(
      "root_collection_deduplicated_rows=",
      nrow(root_collection_summary_deduplicated)
    ),
    paste0(
      "root_metric_deduplicated_rows=",
      nrow(root_metrics_deduplicated_by_root)
    ),
    paste0("paired_comparison_rows=", nrow(paired_sensitivity_table)),
    paste0("paired_sensitivity_rows=", nrow(paired_sensitivity_table)),
    paste0("enclosure_depth_audit_rows=", nrow(enclosure_depth_audit_table)),
    "paired_difference_direction=comparison_variant_minus_reference_variant",
    paste0(
      "root_error_definition=root_bias_conditional_on_isolation and ",
      "root_rmse_conditional_on_isolation use only replications with a ",
      "successful prespecified root isolation and finite root estimate"
    ),
    paste0(
      "root_collection_error_definition=collection bias and RMSE pool every ",
      "root error only over replications with joint isolation of the complete ",
      "prespecified root collection; the replication and root-error ",
      "denominators are reported separately"
    ),
    paste0(
      "legacy_joint_scope=joint_* fields are joint across moment types within ",
      "one fixed interval and are retained for backward compatibility"
    ),
    paste0(
      "across_interval_scope=across_interval_joint_* fields aggregate all ",
      "configured intervals separately by moment type; the ALL_K tables ",
      "deduplicate to one event per successful complete replication"
    ),
    paste0(
      "across_interval_nominal_note=K-specific bands use separately ",
      "calibrated restricted suprema, so the empirical ALL_K rate is a joint ",
      "diagnostic and is not asserted to have new 95 percent nominal coverage"
    ),
    paste0(
      "root_interval_reporting=report_rate, conditional_coverage, and ",
      "report_and_cover_rate are distinct; vacuous_fallback_set_coverage ",
      "assigns the whole real line on nonreport and is not ordinary ",
      "finite-sample coverage"
    ),
    paste0(
      "root_K_deduplication=root-collection and root-specific publication ",
      "tables contain one inferential result per cell/root collection after ",
      "asserting exact agreement of copies repeated across enclosing K sets"
    ),
    paste0(
      "cell_class_policy=benchmark evaluates calibration; power evaluates ",
      "detection/localization; stress documents adverse finite-sample ",
      "behaviour and is not required to attain 95 percent coverage"
    ),
    paste0(
      "tangency_geometry=global K geometry is conditioned on certified ",
      "continuum-band coverage and reports its denominator explicitly; ",
      "legacy local geometry is truncated to a fixed neighbourhood and its ",
      "boundary-censoring rate must accompany it"
    ),
    paste0(
      "enclosure_limit_reporting=legacy enclosure_limit_hit is decomposed ",
      "into depth, node, variance, and statistical flags with the count and ",
      "total, proportional, and maximum width of retained limit cells"
    ),
    "peak_memory_metric=peak_r_heap_mb is the sum of R Ncells and Vcells maxima since gc reset; it is not process RSS",
    paste0("expected_units=", validation$expected_units),
    paste0("numerical_failures=", nrow(numerical_failures)),
    paste0("numerically_clean=", isTRUE(validation$numerically_clean)),
    paste0("cells_with_numerical_failures=", sum(
      failure_summary$numerical_failure_units > 0L
    )),
    "binomial_intervals=two-sided Wilson 95%",
    "zero-event_note=use the count and Wilson interval; plug-in MCSE=0 is not certainty"
  ),
  file.path(summary_directory, "summary_manifest.txt")
)

log_message(
  "Wrote v2 summaries: groups=", nrow(summary_table),
  ", binary_metrics=", nrow(metric_table),
  ", continuous_metrics=", nrow(continuous_table),
  ", across_interval_groups=", nrow(across_interval_summary),
  ", root_metric_rows=", nrow(root_metric_table),
  ", root_collection_deduplicated_rows=",
  nrow(root_collection_summary_deduplicated),
  ", root_metric_deduplicated_rows=",
  nrow(root_metrics_deduplicated_by_root),
  ", paired_comparison_rows=", nrow(paired_sensitivity_table),
  ", enclosure_depth_audit_rows=", nrow(enclosure_depth_audit_table),
  ", directory=", summary_directory
)
