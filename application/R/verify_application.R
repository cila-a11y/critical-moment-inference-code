#!/usr/bin/env Rscript

source(file.path(dirname(normalizePath(sub("^--file=", "", grep(
  "^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1L]),
  winslash = "/", mustWork = TRUE)), "application_common.R"))

arguments <- parse_named_args()
config <- load_hillstrom_config(arguments$config %||% "config/application_v001.R")
data_path <- arguments$data %||% Sys.getenv("HILL_DATA", unset = "")
run_id <- arguments[["run-id"]] %||% "primary"
assert_true(nzchar(data_path), "Provide --data=PATH or set HILL_DATA.")
assert_true(run_id %in% c("primary", "confirmation"),
            "run-id must be primary or confirmation.")

options(warn = 2, digits = 17)
Sys.setenv(TZ = "UTC")
set_single_thread_math()
assert_true(
  identical(paste(R.version$major, R.version$minor, sep = "."),
            config$compute$r_version),
  sprintf("Canonical verification requires R %s.", config$compute$r_version)
)
verify_design_lock()
verify_source_manifest()
source_hillstrom_engine(config)
validated <- read_and_validate_hillstrom(data_path, config, strict = TRUE)

run_directory <- file.path(hill_root, "results", config$application_id, run_id)
assert_true(dir.exists(run_directory),
            sprintf("Run directory not found: %s", run_directory))
required_files <- c(
  scientific_output_files(),
  "data_provenance.csv", "run_metadata.csv", "session_info.txt",
  "application_result_compact.rds", "table_01_hillstrom_application.tex",
  "figure_01_hillstrom_moment_contrast.pdf",
  "figure_01_hillstrom_moment_contrast.png", "MANIFEST_SHA256.txt"
)
missing <- required_files[!file.exists(file.path(run_directory, required_files))]
assert_true(
  !length(missing),
  sprintf("Required result files are missing: %s", paste(missing, collapse = ", "))
)
allowed_files <- c(required_files, "postrun_contract.txt")
unexpected <- setdiff(
  list.files(run_directory, all.files = TRUE, no.. = TRUE, recursive = FALSE),
  allowed_files
)
assert_true(
  !length(unexpected),
  sprintf("Unexpected result files are present: %s",
          paste(unexpected, collapse = ", "))
)
assert_true(verify_result_manifest(run_directory),
            "The result SHA-256 manifest failed.")

read_result <- function(filename) {
  read.csv(file.path(run_directory, filename), stringsAsFactors = FALSE,
           check.names = FALSE, na.strings = c("", "NA"))
}

summary <- read_result("application_summary.csv")
band <- read_result("moment_band.csv")
audit <- read_result("moment_band_audit.csv")
cells <- read_result("enclosure_cells.csv")
outer_components <- read_result("outer_set_components.csv")
sign_intervals <- read_result("certified_sign_intervals.csv")
brackets <- read_result("certified_brackets.csv")
anchors <- read_result("anchor_moments.csv")
descriptive <- read_result("descriptive_statistics.csv")
tail_diagnostics_output <- read_result("tail_diagnostics.csv")
table_csv <- read_result("table_01_hillstrom_application.csv")
bootstrap <- read_result("bootstrap_suprema.csv")
audit_checks <- read_result("audit_checks.csv")
rng_hashes <- read_result("rng_state_hashes.csv")
metadata <- read_result("run_metadata.csv")
validation_output <- read_result("data_validation.csv")
provenance <- read_result("data_provenance.csv")

frame_equal <- function(actual, expected, label, tolerance = 1e-11) {
  assert_true(identical(names(actual), names(expected)),
              sprintf("%s has the wrong columns.", label))
  assert_true(nrow(actual) == nrow(expected),
              sprintf("%s has the wrong number of rows.", label))
  for (column in names(expected)) {
    left <- actual[[column]]
    right <- expected[[column]]
    if (!length(right)) next
    assert_true(identical(is.na(left), is.na(right)),
                sprintf("%s has a different NA pattern in %s.", label, column))
    keep <- !is.na(right)
    if (is.numeric(right) || is.integer(right)) {
      left_numeric <- suppressWarnings(as.numeric(left[keep]))
      right_numeric <- as.numeric(right[keep])
      scale <- pmax(1, abs(left_numeric), abs(right_numeric))
      assert_true(
        all(is.finite(left_numeric)) && all(is.finite(right_numeric)) &&
          all(abs(left_numeric - right_numeric) <= tolerance * scale),
        sprintf("%s disagrees numerically in %s.", label, column)
      )
    } else if (is.logical(right)) {
      assert_true(identical(as.logical(left[keep]), right[keep]),
                  sprintf("%s disagrees logically in %s.", label, column))
    } else {
      assert_true(identical(as.character(left[keep]), as.character(right[keep])),
                  sprintf("%s disagrees textually in %s.", label, column))
    }
  }
  invisible(TRUE)
}

assert_true(nrow(summary) == 1L,
            "application_summary.csv must contain one row.")
analysis <- config$analysis
interval <- analysis$order_intervals$K1
bootstrap_orders <- regular_interval_grid(
  interval, spacing = analysis$bootstrap_grid_spacing
)
enclosure_orders <- regular_interval_grid(
  interval, spacing = analysis$enclosure_grid_spacing
)
audit_orders <- regular_interval_grid(
  interval, spacing = analysis$audit_grid_spacing
)
assert_true(
  nrow(band) == analysis$expected_bootstrap_grid_size &&
    nrow(audit) == analysis$expected_audit_grid_size &&
    nrow(bootstrap) == analysis$bootstrap_reps,
  "A locked grid or bootstrap output has the wrong size."
)
assert_true(
  identical(names(bootstrap), c("replication", "supremum")) &&
    identical(as.integer(bootstrap$replication),
            seq_len(analysis$bootstrap_reps)) &&
    all(is.finite(bootstrap$supremum)),
  "The bootstrap-supremum output is malformed."
)

quantile_index <- min(
  analysis$bootstrap_reps,
  ceiling((analysis$bootstrap_reps + 1L) * (1 - analysis$alpha))
)
recomputed_critical_value <- sort.int(
  bootstrap$supremum, partial = quantile_index
)[quantile_index]
critical_tolerance <- 64 * .Machine$double.eps *
  max(1, abs(summary$critical_value), abs(recomputed_critical_value))
assert_true(
  is.finite(summary$critical_value) && summary$critical_value > 0 &&
    abs(summary$critical_value - recomputed_critical_value) <= critical_tolerance,
  "The simultaneous critical value is not the locked empirical quantile."
)

# Reconstruct every aggregate, the variance guard, and the numerical
# continuum enclosure without rerunning the multiplier bootstrap.
fit <- build_feature_object(
  validated$x, validated$y, bootstrap_orders, analysis$moment_types
)
fit$interval_id <- "K1"
fit$interval <- interval
fit$critical_value <- as.numeric(recomputed_critical_value)
fit$guard_failed <- FALSE
state_cache <- make_state_cache(fit)
constant_cache <- new.env(parent = emptyenv(), hash = TRUE)
guard <- validate_relative_guard(
  object = fit, interval = interval, moment_types = analysis$moment_types,
  tolerance = analysis$variance_tolerance, initial_grid = bootstrap_orders,
  max_levels = analysis$guard_max_levels,
  max_nodes = analysis$enclosure_max_nodes,
  state_cache = state_cache, constant_cache = constant_cache
)
assert_true(isTRUE(guard$valid),
            "The verifier did not reproduce the relative-variance guard.")
fit$guard <- guard

columns <- column_block(fit, analysis$moment_types)
moment_x <- colMeans(fit$feature_x[[analysis$moment_types]])
moment_y <- colMeans(fit$feature_y[[analysis$moment_types]])
variance <- fit$variance[columns]
standard_error <- sqrt(pmax(variance, 0) / fit$n_eff)
half_width <- fit$critical_value * standard_error
expected_band <- data.frame(
  order = fit$orders,
  moment_x = moment_x,
  moment_y = moment_y,
  delta_hat = fit$delta[columns],
  variance_hat = variance,
  standard_error = standard_error,
  critical_value = fit$critical_value,
  half_width = half_width,
  lower = fit$delta[columns] - half_width,
  upper = fit$delta[columns] + half_width,
  studentized = fit$delta[columns] / standard_error,
  sign_status = ifelse(
    fit$delta[columns] - half_width > analysis$numerical_tolerance, "positive",
    ifelse(fit$delta[columns] + half_width < -analysis$numerical_tolerance,
           "negative", "unresolved")
  ),
  stringsAsFactors = FALSE
)
frame_equal(band, expected_band, "moment_band.csv")

evaluated_audit <- evaluate_fitted_band(
  fit, audit_orders, analysis$moment_types,
  chunk_size = analysis$evaluation_chunk_size
)
audit_standard_error <- sqrt(pmax(evaluated_audit$variance_hat, 0) / fit$n_eff)
expected_audit <- data.frame(
  order = evaluated_audit$orders,
  delta_hat = evaluated_audit$delta_hat,
  variance_hat = evaluated_audit$variance_hat,
  standard_error = audit_standard_error,
  critical_value = fit$critical_value,
  half_width = evaluated_audit$half_width,
  lower = evaluated_audit$lower,
  upper = evaluated_audit$upper,
  studentized = evaluated_audit$delta_hat / audit_standard_error,
  sign_status = ifelse(
    evaluated_audit$lower > analysis$numerical_tolerance, "positive",
    ifelse(evaluated_audit$upper < -analysis$numerical_tolerance,
           "negative", "unresolved")
  ),
  stringsAsFactors = FALSE
)
frame_equal(audit, expected_audit, "moment_band_audit.csv")

enclosure <- adaptive_continuum_enclosure(
  fit = fit, moment_type = analysis$moment_types,
  initial_grid = enclosure_orders,
  numerical_tolerance = analysis$numerical_tolerance,
  safety_margin = analysis$enclosure_safety_margin,
  max_levels = analysis$enclosure_max_levels,
  max_nodes = analysis$enclosure_max_nodes,
  refine_unresolved_levels = analysis$enclosure_refine_unresolved_levels,
  maximum_cell_half_width = analysis$maximum_enclosure_half_width,
  roundoff_inflation = analysis$roundoff_inflation,
  state_cache = state_cache, constant_cache = constant_cache
)
frame_equal(cells, enclosure$cells, "enclosure_cells.csv")
cells_expected <- enclosure$cells

partition_tolerance <- 1e-12
assert_true(nrow(cells_expected) >= 1L &&
              all(cells_expected$right > cells_expected$left),
            "The continuum-cell table is malformed.")
assert_true(
  abs(cells_expected$left[1L] - interval[1L]) <= partition_tolerance &&
    abs(tail(cells_expected$right, 1L) - interval[2L]) <=
      partition_tolerance,
  "The continuum cells do not cover K1."
)
if (nrow(cells_expected) > 1L) {
  assert_true(
    all(abs(head(cells_expected$right, -1L) -
              tail(cells_expected$left, -1L)) <= partition_tolerance),
    "The continuum cells contain a gap or overlap."
  )
}
assert_true(
  all(rowSums(cbind(cells_expected$positive, cells_expected$negative,
                    cells_expected$outer)) == 1L) &&
    all(cells_expected$lower_bound <= cells_expected$upper_bound) &&
    !any(cells_expected$variance_unresolved_at_limit),
  "The continuum-cell classifications or bounds are invalid."
)

audit_containment <- vapply(seq_len(nrow(expected_audit)), function(index) {
  order <- expected_audit$order[index]
  candidates <- which(cells_expected$left <= order + partition_tolerance &
                        cells_expected$right >= order - partition_tolerance)
  length(candidates) >= 1L &&
    expected_audit$lower[index] >=
      max(cells_expected$lower_bound[candidates]) - partition_tolerance &&
    expected_audit$upper[index] <=
      min(cells_expected$upper_bound[candidates]) + partition_tolerance
}, logical(1L))
audit_sign_consistent <- vapply(seq_len(nrow(expected_audit)), function(index) {
  order <- expected_audit$order[index]
  candidates <- which(cells_expected$left <= order + partition_tolerance &
                        cells_expected$right >= order - partition_tolerance)
  length(candidates) >= 1L && all(vapply(candidates, function(cell_index) {
    (!cells_expected$positive[cell_index] ||
       expected_audit$lower[index] > analysis$numerical_tolerance) &&
      (!cells_expected$negative[cell_index] ||
         expected_audit$upper[index] < -analysis$numerical_tolerance)
  }, logical(1L)))
}, logical(1L))
expected_audit_checks <- data.frame(
  check = c("audit_nodes", "whole_cell_containment_violations",
            "cell_sign_contradictions"),
  value = c(nrow(expected_audit), sum(!audit_containment),
            sum(!audit_sign_consistent)),
  status = c(
    if (nrow(expected_audit) == analysis$expected_audit_grid_size)
      "PASSED" else "FAILED",
    if (all(audit_containment)) "PASSED" else "FAILED",
    if (all(audit_sign_consistent)) "PASSED" else "FAILED"
  ),
  stringsAsFactors = FALSE
)
frame_equal(audit_checks, expected_audit_checks, "audit_checks.csv")
assert_true(all(expected_audit_checks$status == "PASSED"),
            "The separate fine audit grid contradicted the numerical enclosure.")

outer <- outer_set_summary(cells_expected, true_roots = numeric(0L))
expected_outer_components <- outer$components
expected_outer_components$component <- seq_len(nrow(expected_outer_components))
expected_outer_components$width <- expected_outer_components$right -
  expected_outer_components$left
expected_outer_components <- expected_outer_components[
  , c("component", "left", "right", "width"), drop = FALSE
]
frame_equal(outer_components, expected_outer_components,
            "outer_set_components.csv")

positive_intervals <- merge_adjacent_cells(
  cells_expected, cells_expected$positive, "positive"
)
negative_intervals <- merge_adjacent_cells(
  cells_expected, cells_expected$negative, "negative"
)
expected_sign_intervals <- rbind(positive_intervals, negative_intervals)
if (nrow(expected_sign_intervals)) {
  expected_sign_intervals <- expected_sign_intervals[
    order(expected_sign_intervals$left, expected_sign_intervals$right),
    , drop = FALSE
  ]
  expected_sign_intervals$width <- expected_sign_intervals$right -
    expected_sign_intervals$left
} else {
  expected_sign_intervals$width <- numeric(0L)
}
expected_brackets <- brackets_from_sign_intervals(expected_sign_intervals)
assert_bracket_invariants(
  expected_brackets, expected_sign_intervals, interval,
  tolerance = partition_tolerance
)
frame_equal(sign_intervals, expected_sign_intervals,
            "certified_sign_intervals.csv")
frame_equal(brackets, expected_brackets, "certified_brackets.csv")

describe_group <- function(values, segment, role) {
  positive <- values[values > 0]
  data.frame(
    role = role, segment = segment, n = length(values),
    n_positive = length(positive), positive_rate = mean(values > 0),
    mean_spend = mean(values), sd_spend = stats::sd(values),
    median_spend = stats::median(values),
    conditional_mean_positive = if (length(positive)) mean(positive) else NA_real_,
    maximum_spend = max(values), stringsAsFactors = FALSE
  )
}
expected_descriptive <- rbind(
  describe_group(validated$x, analysis$x_segment, "X"),
  describe_group(validated$y, analysis$y_segment, "Y")
)
frame_equal(descriptive, expected_descriptive, "descriptive_statistics.csv")
frame_equal(validation_output, validated$validation, "data_validation.csv")

expected_tail_diagnostics <- do.call(rbind, lapply(c("X", "Y"), function(role) {
  values <- if (identical(role, "X")) validated$x else validated$y
  segment <- if (identical(role, "X")) analysis$x_segment else
    analysis$y_segment
  do.call(rbind, lapply(analysis$tail_diagnostic_orders, function(order) {
    diagnostic_targets <- c("estimand_moment", "second_moment")
    diagnostic_exponents <- c(order, 2 * order)
    do.call(rbind, lapply(seq_along(diagnostic_targets), function(target_index) {
      diagnostic_target <- diagnostic_targets[target_index]
      power_exponent <- diagnostic_exponents[target_index]
      diagnostic <- moment_weight_diagnostic(
        values, order = power_exponent,
        top_fraction = analysis$tail_top_fraction
      )
      magnitude <- abs(values)
      log_weight <- power_exponent * log(magnitude)
      log_weight[magnitude == 0] <- -Inf
      maximum <- max(log_weight[is.finite(log_weight)])
      weights <- exp(log_weight - maximum)
      weights[!is.finite(weights)] <- 0
      shares <- sort(weights / sum(weights), decreasing = TRUE)
      n_positive <- sum(values > 0)
      share_top <- function(number) {
        sum(head(shares, min(number, length(shares))))
      }
      data.frame(
        role = role, segment = segment, order = order,
        diagnostic_target = diagnostic_target,
        power_exponent = power_exponent,
        n = length(values), n_positive = n_positive,
        ess = unname(diagnostic[["ess"]]),
        ess_ratio = unname(diagnostic[["ess_ratio"]]),
        ess_ratio_among_positive = unname(diagnostic[["ess"]]) / n_positive,
        maximum_share = unname(diagnostic[["maximum_share"]]),
        top_1pct_share = unname(diagnostic[["top_share"]]),
        top_1_contributor_share = share_top(1L),
        top_5_contributors_share = share_top(5L),
        top_10_contributors_share = share_top(10L),
        top_10pct_positive_share = share_top(
          max(1L, ceiling(0.10 * n_positive))
        ),
        stringsAsFactors = FALSE
      )
    }))
  }))
}))
frame_equal(tail_diagnostics_output, expected_tail_diagnostics,
            "tail_diagnostics.csv")

anchor_indices <- vapply(analysis$anchor_orders, function(order) {
  index <- which.min(abs(expected_band$order - order))
  assert_true(abs(expected_band$order[index] - order) <= 1e-12,
              "An anchor order is absent from the construction grid.")
  index
}, integer(1L))
expected_anchors <- do.call(rbind, lapply(
  seq_along(analysis$anchor_orders), function(index) {
    order <- analysis$anchor_orders[index]
    candidates <- which(cells_expected$left <= order + partition_tolerance &
                          cells_expected$right >= order - partition_tolerance)
    assert_true(length(candidates) >= 1L,
                "An anchor order is absent from the continuum partition.")
    endpoint_lower <- max(cells_expected$lower_bound[candidates])
    endpoint_upper <- min(cells_expected$upper_bound[candidates])
    assert_true(endpoint_lower <= endpoint_upper + partition_tolerance,
                "Whole-cell endpoint enclosures disagree at an anchor order.")
    status <- if (endpoint_lower > analysis$numerical_tolerance) {
      "positive"
    } else if (endpoint_upper < -analysis$numerical_tolerance) {
      "negative"
    } else {
      "outer"
    }
    data.frame(
      order = order,
      moment_x = expected_band$moment_x[anchor_indices[index]],
      moment_y = expected_band$moment_y[anchor_indices[index]],
      delta_hat = expected_band$delta_hat[anchor_indices[index]],
      implemented_lower = expected_band$lower[anchor_indices[index]],
      implemented_upper = expected_band$upper[anchor_indices[index]],
      endpoint_lower = endpoint_lower,
      endpoint_upper = endpoint_upper,
      anchor_sign_status = status,
      stringsAsFactors = FALSE
    )
  }
))
frame_equal(anchors, expected_anchors, "anchor_moments.csv")

expected_table_csv <- rbind(
  data.frame(
    panel = "A", role = expected_descriptive$role,
    segment = expected_descriptive$segment, order = NA_real_,
    n = expected_descriptive$n,
    n_positive = expected_descriptive$n_positive,
    positive_rate = expected_descriptive$positive_rate,
    mean_spend = expected_descriptive$mean_spend,
    sd_spend = expected_descriptive$sd_spend,
    conditional_mean_positive =
      expected_descriptive$conditional_mean_positive,
    moment_x = NA_real_, moment_y = NA_real_, delta_hat = NA_real_,
    implemented_lower = NA_real_, implemented_upper = NA_real_,
    endpoint_lower = NA_real_, endpoint_upper = NA_real_,
    anchor_sign_status = NA_character_, stringsAsFactors = FALSE
  ),
  data.frame(
    panel = "B", role = NA_character_, segment = NA_character_,
    order = expected_anchors$order, n = NA_integer_, n_positive = NA_integer_,
    positive_rate = NA_real_, mean_spend = NA_real_,
    sd_spend = NA_real_,
    conditional_mean_positive = NA_real_,
    moment_x = expected_anchors$moment_x,
    moment_y = expected_anchors$moment_y,
    delta_hat = expected_anchors$delta_hat,
    implemented_lower = expected_anchors$implemented_lower,
    implemented_upper = expected_anchors$implemented_upper,
    endpoint_lower = expected_anchors$endpoint_lower,
    endpoint_upper = expected_anchors$endpoint_upper,
    anchor_sign_status = expected_anchors$anchor_sign_status,
    stringsAsFactors = FALSE
  )
)
frame_equal(table_csv, expected_table_csv,
            "table_01_hillstrom_application.csv")

uniform_y <- nrow(cells_expected) > 0L && all(cells_expected$positive)
uniform_x <- nrow(cells_expected) > 0L && all(cells_expected$negative)
empty_outer <- nrow(cells_expected) > 0L && !any(cells_expected$outer)
assert_true(!empty_outer || xor(uniform_y, uniform_x),
            "An empty outer set lacks a uniform certified sign.")
no_root <- empty_outer && (uniform_y || uniform_x)
expected_conclusion <- if (uniform_y) {
  "UNIFORM_Y_DOMINANCE_NO_ROOT"
} else if (uniform_x) {
  "UNIFORM_X_DOMINANCE_NO_ROOT"
} else if (nrow(expected_brackets) > 0L) {
  "REVERSAL_AT_LEAST_ONE_ROOT"
} else {
  "INCONCLUSIVE"
}
expected_summary <- data.frame(
  application_id = config$application_id,
  x_segment = analysis$x_segment, y_segment = analysis$y_segment,
  n_x = length(validated$x), n_y = length(validated$y), n_eff = fit$n_eff,
  p_min = interval[1L], p_max = interval[2L],
  alpha = analysis$alpha, bootstrap_reps = analysis$bootstrap_reps,
  bootstrap_seed = analysis$bootstrap_seed,
  multiplier_distribution = analysis$multiplier_distribution,
  bootstrap_grid_size = length(bootstrap_orders),
  enclosure_initial_grid_size = length(enclosure_orders),
  audit_grid_size = length(audit_orders),
  critical_value = fit$critical_value,
  guard_valid = guard$valid,
  guard_minimum_certified_ratio = guard$minimum_certified_ratio,
  guard_nodes_used = guard$nodes_used,
  guard_maximum_level = guard$maximum_level,
  enclosure_nodes_used = enclosure$nodes_used,
  enclosure_maximum_level = enclosure$maximum_level,
  enclosure_limit_hit = enclosure$limit_hit,
  enclosure_variance_limit_hit = enclosure$variance_limit_hit,
  enclosure_statistical_limit_hit = enclosure$statistical_limit_hit,
  outer_component_count = outer$component_count,
  outer_total_length = outer$total_length,
  certified_alternations = nrow(expected_brackets),
  certified_bracket_count = nrow(expected_brackets),
  no_root_on_K_certified = no_root,
  uniform_y_dominance_certified = uniform_y,
  uniform_x_dominance_certified = uniform_x,
  conclusion = expected_conclusion,
  structural_crossing_budget_used = FALSE,
  stringsAsFactors = FALSE
)
frame_equal(summary, expected_summary, "application_summary.csv")

assert_true(
  identical(names(rng_hashes), c("state", "sha256")) &&
    nrow(rng_hashes) == 2L &&
    identical(as.character(rng_hashes$state), c("initial", "final")) &&
    all(grepl("^[0-9a-f]{64}$", as.character(rng_hashes$sha256))),
  "The RNG-state audit is malformed."
)

expected_provenance_columns <- c(
  "dataset", "source_mode", "source_used", "source_page", "accessed_utc",
  "raw_filename", "bytes", "sha256", "license_status", "redistributed"
)
assert_true(
  nrow(provenance) == 1L &&
    identical(names(provenance), expected_provenance_columns) &&
    identical(as.character(provenance$sha256), config$data$sha256) &&
    as.numeric(provenance$bytes) == as.numeric(config$data$bytes) &&
    identical(as.character(provenance$raw_filename),
              config$data$raw_filename) &&
    identical(toupper(as.character(provenance$redistributed)), "FALSE"),
  "The result data-provenance row violates the locked contract."
)

expected_metadata_fields <- c(
  "application_id", "application_version", "run_id", "started_utc",
  "finished_utc", "host", "slurm_job_id", "r_version", "platform",
  "engine_version", "engine_manifest_sha256", "config_sha256",
  "source_manifest_sha256", "design_lock_sha256", "data_sha256",
  "data_bytes", "rng_kind", "bootstrap_seed", "initial_rng_state_sha256",
  "final_rng_state_sha256"
)
assert_true(
  identical(names(metadata), c("field", "value")) &&
    identical(as.character(metadata$field), expected_metadata_fields) &&
    !anyDuplicated(metadata$field),
  "run_metadata.csv is malformed."
)
metadata_value <- function(field) {
  value <- metadata$value[match(field, metadata$field)]
  assert_true(length(value) == 1L && !is.na(value),
              sprintf("Metadata field is missing: %s", field))
  as.character(value)
}
assert_true(
  identical(metadata_value("application_id"), config$application_id) &&
    identical(metadata_value("application_version"),
              config$application_version) &&
    identical(metadata_value("run_id"), run_id) &&
    identical(metadata_value("r_version"), R.version.string) &&
    identical(metadata_value("engine_version"), config$engine$version) &&
    identical(metadata_value("engine_manifest_sha256"),
              config$engine$manifest_sha256) &&
    identical(metadata_value("config_sha256"),
              sha256_file(config$config_path)) &&
    identical(metadata_value("source_manifest_sha256"),
              sha256_file(file.path(hill_root,
                                    "SOURCE_MANIFEST_SHA256.txt"))) &&
    identical(metadata_value("design_lock_sha256"),
              sha256_file(file.path(hill_root, "DESIGN_LOCK_SHA256.txt"))) &&
    identical(metadata_value("data_sha256"), config$data$sha256) &&
    as.numeric(metadata_value("data_bytes")) == as.numeric(config$data$bytes) &&
    identical(metadata_value("rng_kind"),
              paste(analysis$rng_kind, collapse = ";")) &&
    as.integer(metadata_value("bootstrap_seed")) ==
      analysis$bootstrap_seed,
  "The run metadata disagrees with a locked source, data, RNG, or runtime field."
)
assert_true(
  identical(metadata_value("initial_rng_state_sha256"),
            as.character(rng_hashes$sha256[rng_hashes$state == "initial"])) &&
    identical(metadata_value("final_rng_state_sha256"),
              as.character(rng_hashes$sha256[rng_hashes$state == "final"])),
  "The metadata and RNG-state table disagree."
)

compact <- readRDS(file.path(run_directory,
                             "application_result_compact.rds"))
allowed_compact_names <- c(
  "config", "application_summary", "descriptive_statistics", "moment_band",
  "moment_band_audit", "enclosure_cells", "outer_set_components",
  "certified_sign_intervals", "certified_brackets", "tail_diagnostics",
  "audit_checks", "rng_state_hashes"
)
assert_true(identical(names(compact), allowed_compact_names),
            "The compact RDS has unexpected fields and may contain raw data.")
assert_true(identical(compact$config, config),
            "The compact RDS contains the wrong configuration.")
compact_map <- list(
  application_summary = summary,
  descriptive_statistics = descriptive,
  moment_band = band,
  moment_band_audit = audit,
  enclosure_cells = cells,
  outer_set_components = outer_components,
  certified_sign_intervals = sign_intervals,
  certified_brackets = brackets,
  tail_diagnostics = tail_diagnostics_output,
  audit_checks = audit_checks,
  rng_state_hashes = rng_hashes
)
for (name in names(compact_map)) {
  frame_equal(compact[[name]], compact_map[[name]],
              sprintf("compact RDS field %s", name))
}
contains_raw_schema <- function(object) {
  if (is.data.frame(object)) {
    return(all(config$data$expected_columns %in% names(object)))
  }
  if (is.list(object)) {
    return(any(vapply(object, contains_raw_schema, logical(1L))))
  }
  FALSE
}
assert_true(!contains_raw_schema(compact),
            "The compact RDS contains a data frame with the raw-data schema.")

article_files <- file.path(
  run_directory,
  c("table_01_hillstrom_application.tex",
    "figure_01_hillstrom_moment_contrast.pdf",
    "figure_01_hillstrom_moment_contrast.png")
)
assert_true(all(file.info(article_files)$size > 100L),
            "An article table or figure is empty.")
pdf_header <- readBin(article_files[2L], what = "raw", n = 4L)
png_header <- readBin(article_files[3L], what = "raw", n = 8L)
assert_true(
  identical(rawToChar(pdf_header), "%PDF") &&
    identical(as.integer(png_header),
              c(137L, 80L, 78L, 71L, 13L, 10L, 26L, 10L)),
  "A publication figure has the wrong file signature."
)
table_lines <- readLines(article_files[1L], warn = FALSE)
expected_table_lines <- render_hillstrom_table(
  expected_descriptive, expected_anchors, expected_summary,
  expected_outer_components, expected_brackets
)
assert_true(
  identical(table_lines, expected_table_lines) &&
    any(grepl("Panel A", table_lines, fixed = TRUE)) &&
    any(grepl("Panel B", table_lines, fixed = TRUE)) &&
    any(grepl("Panel C", table_lines, fixed = TRUE)),
  "The publication table is not the exact rendering of its audited inputs."
)

contract_lines <- c(
  sprintf("application_id=%s", config$application_id),
  sprintf("application_version=%s", config$application_version),
  sprintf("run_id=%s", run_id),
  sprintf("engine_version=%s", config$engine$version),
  sprintf("source_manifest_sha256=%s",
          sha256_file(file.path(hill_root, "SOURCE_MANIFEST_SHA256.txt"))),
  sprintf("design_lock_sha256=%s",
          sha256_file(file.path(hill_root, "DESIGN_LOCK_SHA256.txt"))),
  sprintf("data_sha256=%s", config$data$sha256),
  sprintf("n_x=%d", summary$n_x),
  sprintf("n_y=%d", summary$n_y),
  sprintf("bootstrap_reps=%d", summary$bootstrap_reps),
  sprintf("bootstrap_quantile_index=%d", quantile_index),
  sprintf("bootstrap_grid_size=%d", nrow(band)),
  sprintf("audit_grid_size=%d", nrow(audit)),
  sprintf("continuum_cells=%d", nrow(cells)),
  sprintf("guard_minimum_certified_ratio=%.17g",
          summary$guard_minimum_certified_ratio),
  sprintf("enclosure_statistical_limit_hit=%s",
          summary$enclosure_statistical_limit_hit),
  sprintf("scientific_conclusion=%s", summary$conclusion),
  "simultaneous_control=nominal_95_percent_asymptotic_under_article_conditions",
  "selected_contrast_post_selection_correction=FALSE",
  "finite_sample_guarantee=FALSE",
  "raw_or_row_level_data_in_outputs=FALSE",
  "HILLSTROM_APPLICATION_V001 POST-RUN CONTRACT: PASSED"
)
status_directory <- file.path(hill_root, "status", config$application_id,
                              run_id)
dir.create(status_directory, recursive = TRUE, showWarnings = FALSE)
atomic_write_lines(
  contract_lines, file.path(run_directory, "postrun_contract.txt"),
  overwrite = TRUE
)
write_result_manifest(run_directory)
assert_true(
  verify_result_manifest(run_directory) &&
    any(grepl("^[0-9a-f]{64}  \\./postrun_contract\\.txt$", readLines(
      file.path(run_directory, "MANIFEST_SHA256.txt"), warn = FALSE
    ))),
  "The final result manifest does not anchor the passed post-run contract."
)
atomic_write_lines(
  contract_lines, file.path(status_directory, "postrun_contract.txt"),
  overwrite = TRUE
)
assert_true(
  file.exists(file.path(status_directory, "postrun_contract.txt")) &&
    any(readLines(
      file.path(status_directory, "postrun_contract.txt"), warn = FALSE
    ) == "HILLSTROM_APPLICATION_V001 POST-RUN CONTRACT: PASSED"),
  "The external status contract was not published after result verification."
)

cat(paste(contract_lines, collapse = "\n"), "\n")
