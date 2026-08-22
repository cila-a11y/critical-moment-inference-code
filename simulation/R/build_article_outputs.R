#!/usr/bin/env Rscript

# Build a compact, publication-facing layer from the validated study summaries.
# This script deliberately never re-aggregates replication rows.  In
# particular, regular-root outputs are read from the K-deduplicated CSV files
# written by summarize_results.R, so repeated copies of one root calculation
# on K1/K2/K3 cannot be mistaken for independent Monte Carlo evidence.

source(file.path(Sys.getenv("CMO_ROOT", unset = "."), "R", "common.R"))
source_cmo("R/dgp.R")

arguments <- parse_named_args()
config_argument <- arguments$config %||% Sys.getenv("CMO_CONFIG", unset = "")
assert_true(nzchar(config_argument), "Supply --config=... or set CMO_CONFIG.")
config <- read_config(config_argument)

summary_directory <- file.path(
  cmo_root(), "results", "summary", config$study_id
)
article_directory <- file.path(summary_directory, "article")
validation_path <- file.path(cmo_root(), "status", config$study_id,
                             "validation.rds")
assert_true(file.exists(validation_path), "Run validate_results.R first.")
validation <- readRDS(validation_path)
assert_true(isTRUE(validation$complete),
            "Article outputs require a complete validated study.")

required_summary_files <- c(
  summary = "summary.csv",
  binary = "binary_metrics.csv",
  continuous = "continuous_metrics.csv",
  across = "across_interval_summary.csv",
  across_binary = "across_interval_binary_metrics.csv",
  root_collection = "root_collection_summary_deduplicated.csv",
  root_specific = "root_metrics_deduplicated_by_root.csv",
  failures = "numerical_failure_summary_by_cell.csv"
)
required_paths <- file.path(summary_directory, required_summary_files)
assert_true(
  all(file.exists(required_paths)),
  sprintf(
    "Missing validated summary files: %s.",
    paste(required_summary_files[!file.exists(required_paths)], collapse = ", ")
  )
)

read_summary_csv <- function(name) {
  read.csv(
    file.path(summary_directory, required_summary_files[[name]]),
    stringsAsFactors = FALSE, check.names = FALSE
  )
}

summary_table <- read_summary_csv("summary")
binary_table <- read_summary_csv("binary")
continuous_table <- read_summary_csv("continuous")
across_table <- read_summary_csv("across")
across_binary_table <- read_summary_csv("across_binary")
root_collection <- read_summary_csv("root_collection")
root_specific <- read_summary_csv("root_specific")
failure_table <- read_summary_csv("failures")

assert_true(nrow(summary_table) > 0L, "The scientific summary is empty.")
assert_true(nrow(binary_table) > 0L && nrow(continuous_table) > 0L,
            "Metric summaries are empty.")
assert_true(nrow(across_table) > 0L && nrow(across_binary_table) > 0L,
            "The across-interval summary is empty.")

required_summary_columns <- c(
  "study_id", "cell_id", "cell_class", "variant_id", "scenario_id",
  "n_x", "n_y", "moment_type", "interval_id",
  "expected_replications", "successful_replications",
  "continuum_band_coverage", "any_primary_false_report_rate",
  "exact_count_certification_rate", "reversal_certification_rate",
  "outer_set_truth_coverage", "mean_outer_set_length"
)
assert_true(
  all(required_summary_columns %in% names(summary_table)),
  sprintf(
    "The scientific summary lacks publication fields: %s.",
    paste(setdiff(required_summary_columns, names(summary_table)), collapse = ", ")
  )
)

# A single, explicit name is used for the familywise event in publication
# tables.  Keep the engine-level name as an exact alias for auditability.
if (!"familywise_primary_error_rate" %in% names(summary_table)) {
  summary_table$familywise_primary_error_rate <-
    summary_table$any_primary_false_report_rate
}
alias_comparable <- is.finite(summary_table$familywise_primary_error_rate) |
  is.finite(summary_table$any_primary_false_report_rate)
assert_true(
  all(
    abs(
      summary_table$familywise_primary_error_rate[alias_comparable] -
        summary_table$any_primary_false_report_rate[alias_comparable]
    ) <= 1e-14
  ),
  "The familywise-primary-error alias does not equal any-primary-false-report."
)

select_existing <- function(data, columns) {
  data[, intersect(columns, names(data)), drop = FALSE]
}

collapse_unique <- function(values, separator = "; ") {
  values <- unique(as.character(values))
  values <- values[!is.na(values) & nzchar(values)]
  paste(values, collapse = separator)
}

format_sample_pair <- function(cell) {
  sprintf("%d/%d", as.integer(cell$n_x), as.integer(cell$n_y))
}

count_root_intervals <- function(value) {
  if (is.null(value) || !length(value)) return(0L)
  if (is.numeric(value) && length(value) == 2L) return(1L)
  if (!is.list(value)) return(0L)
  sum(vapply(value, count_root_intervals, integer(1L)))
}

root_mode_for_cells <- function(cells) {
  scenario <- vapply(cells, `[[`, character(1L), "scenario_id")
  if (any(grepl("TANGENCY", scenario, fixed = TRUE))) {
    return("outer set only (tangential root)")
  }
  counts <- vapply(cells, function(cell) {
    count_root_intervals(cell$root_intervals %||% list())
  }, integer(1L))
  if (any(counts > 0L)) {
    return("prespecified simple-root interval(s)")
  }
  "no regular-root refinement"
}

format_population_truth <- function(scenario, moment_types) {
  paste(vapply(moment_types, function(moment_type) {
    if (isTRUE(scenario$identically_zero_by_type[[moment_type]])) {
      return(sprintf("%s: equality at every order", moment_type))
    }
    truth <- scenario$root_truth[[moment_type]]
    if (!nrow(truth)) return(sprintf("%s: no root", moment_type))
    roots <- vapply(seq_len(nrow(truth)), function(index) {
      sprintf(
        "p=%s (multiplicity %d, %s)",
        format(truth$order[index], trim = TRUE, digits = 6),
        truth$multiplicity[index], truth$direction[index]
      )
    }, character(1L))
    sprintf("%s: %s", moment_type, paste(roots, collapse = "; "))
  }, character(1L)), collapse = " | ")
}

format_structural_budgets <- function(scenario, moment_types) {
  paste(vapply(moment_types, function(moment_type) {
    budget <- scenario$structural_budgets[[moment_type]]
    sprintf("%s: %s", moment_type, if (length(budget) == 1L &&
                                          is.finite(budget)) {
      as.character(as.integer(budget))
    } else "unbounded")
  }, character(1L)), collapse = " | ")
}

# Table 1: one row per DGP and scientific class, preserving first appearance
# in the frozen config.  Keeping class in the key protects the table if a DGP
# is intentionally used as both a power design and a stress design.
scenario_values <- vapply(config$cells, `[[`, character(1L), "scenario_id")
class_values <- vapply(config$cells, function(cell) {
  cell$cell_class %||% "unclassified"
}, character(1L))
design_keys <- unique(paste(scenario_values, class_values, sep = "::"))
design_rows <- lapply(design_keys, function(design_key) {
  key_parts <- strsplit(design_key, "::", fixed = TRUE)[[1L]]
  scenario_id <- key_parts[1L]
  scientific_class <- key_parts[2L]
  cells <- config$cells[vapply(config$cells, function(cell) {
    identical(cell$scenario_id, scenario_id) &&
      identical(cell$cell_class %||% "unclassified", scientific_class)
  }, logical(1L))]
  sample_pairs <- unique(vapply(cells, format_sample_pair, character(1L)))
  interval_ids <- unique(unlist(lapply(cells, `[[`, "interval_ids"),
                                use.names = FALSE))
  moment_types <- unique(unlist(lapply(cells, `[[`, "moment_types"),
                                 use.names = FALSE))
  scenario <- get_scenario(scenario_id)
  cell_classes <- unique(vapply(cells, function(cell) {
    cell$cell_class %||% "unclassified"
  }, character(1L)))
  data.frame(
    scenario_id = scenario_id,
    scenario_label = scenario$label,
    distribution_family = scenario$family,
    scientific_class = collapse_unique(cell_classes),
    sample_sizes_nx_over_ny = paste(sample_pairs, collapse = "; "),
    moment_types = paste(moment_types, collapse = ";"),
    order_intervals = paste(interval_ids, collapse = ";"),
    population_root_truth = format_population_truth(scenario, moment_types),
    structural_multiplicity_budget =
      format_structural_budgets(scenario, moment_types),
    root_inference = root_mode_for_cells(cells),
    replications_per_cell = as.integer(config$mc_reps_per_cell),
    bootstrap_repetitions = as.integer(config$bootstrap_reps),
    stringsAsFactors = FALSE
  )
})
table_01_design <- do.call(rbind, design_rows)
row.names(table_01_design) <- NULL

across_key <- function(data) paste(data$cell_id, data$moment_type, sep = "::")

metric_details <- function(metric) {
  rows <- across_binary_table[across_binary_table$metric == metric, , drop = FALSE]
  keys <- across_key(rows)
  assert_true(!anyDuplicated(keys),
              sprintf("Across-interval metric %s is duplicated.", metric))
  positions <- match(across_key(table_02_source), keys)
  assert_true(!anyNA(positions), sprintf(
    "Across-interval metric %s is missing an editorial row.", metric
  ))
  data.frame(
    successes = rows$successes[positions],
    denominator = rows$denominator[positions],
    estimate = rows$estimate[positions],
    mcse = rows$mcse[positions],
    lower = rows$wilson_lower[positions],
    upper = rows$wilson_upper[positions]
  )
}

interval_metric_columns <- function(field, prefix) {
  keys <- paste(summary_table$cell_id, summary_table$moment_type,
                summary_table$interval_id, sep = "::")
  assert_true(!anyDuplicated(keys),
              sprintf("Interval metric %s has duplicated rows.", field))
  output <- list()
  for (interval_id in names(config$order_intervals)) {
    target <- paste(table_02_source$cell_id, table_02_source$moment_type,
                    interval_id, sep = "::")
    output[[paste0(interval_id, "_", prefix)]] <-
      summary_table[[field]][match(target, keys)]
  }
  as.data.frame(output, stringsAsFactors = FALSE, check.names = FALSE)
}

# Table 2: benchmark rows are simultaneous over all configured K intervals.
# The smoke has diagnostic rather than benchmark classes; use its diagnostic
# rows so that the publication pipeline itself is exercised end to end.
table_02_source <- across_table[across_table$cell_class == "benchmark", ,
                                drop = FALSE]
if (!nrow(table_02_source) &&
    identical(config$article_output_profile %||% "", "smoke_bridge8")) {
  table_02_source <- across_table
}
continuum_details <- metric_details(
  "across_interval_joint_continuum_band_covers_truth"
)
familywise_details <- metric_details(
  "across_interval_joint_any_primary_false_report"
)
table_02_benchmark <- data.frame(
  select_existing(table_02_source, c(
    "study_id", "cell_id", "cell_class", "scenario_id", "n_x", "n_y",
    "moment_type", "configured_interval_ids", "expected_replications",
    "successful_replications", "complete_interval_replications"
  )),
  simultaneous_continuum_band_successes = continuum_details$successes,
  simultaneous_continuum_band_denominator = continuum_details$denominator,
  simultaneous_continuum_band_coverage =
    continuum_details$estimate,
  simultaneous_continuum_band_coverage_mcse = continuum_details$mcse,
  simultaneous_continuum_band_coverage_wilson_lower = continuum_details$lower,
  simultaneous_continuum_band_coverage_wilson_upper = continuum_details$upper,
  familywise_primary_error_count = familywise_details$successes,
  familywise_primary_error_denominator = familywise_details$denominator,
  any_primary_false_report_rate =
    familywise_details$estimate,
  familywise_primary_error_rate = familywise_details$estimate,
  familywise_primary_error_mcse = familywise_details$mcse,
  familywise_primary_error_wilson_lower = familywise_details$lower,
  familywise_primary_error_wilson_upper = familywise_details$upper,
  interval_metric_columns(
    "mean_outer_set_retained_proportion", "outer_set_retained_proportion"
  ),
  interval_metric_columns(
    "continuum_no_root_certification_rate", "no_root_certification_rate"
  ),
  interval_metric_columns(
    "exact_count_certification_rate", "exact_count_certification_rate"
  ),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
assert_true(
  isTRUE(all.equal(
    table_02_benchmark$any_primary_false_report_rate,
    table_02_benchmark$familywise_primary_error_rate,
    tolerance = 0
  )),
  "The editorial familywise-error alias changed the underlying event."
)

# Tables 3 and 3b come only from K-deduplicated sources.  Table 3 is long:
# each root collection contributes one Wald row and one multiplier-bootstrap
# row, so report rate, conditional coverage, and report-and-cover are visually
# distinct and cannot be confused by sharing a wide block of columns.
root_collection_columns <- c(
  "study_id", "cell_id", "cell_class", "variant_id", "scenario_id",
  "n_x", "n_y", "moment_type",
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
  "root_multiplier_report_and_cover_count",
  "source_interval_ids", "source_interval_count"
)
assert_true(
  all(root_collection_columns %in% names(root_collection)),
  "The K-deduplicated root-collection table is incomplete."
)
root_identity_columns <- c(
  "study_id", "cell_id", "cell_class", "variant_id", "scenario_id",
  "n_x", "n_y", "moment_type",
  "root_collection_applicable_replications",
  "root_collection_joint_isolation_reports",
  "root_collection_joint_isolation_report_rate",
  "root_collection_conditional_replication_denominator",
  "root_collection_error_denominator",
  "root_collection_bias_conditional_on_joint_isolation",
  "root_collection_rmse_conditional_on_joint_isolation",
  "source_interval_ids", "source_interval_count"
)
root_collection_continuous_details <- function(metric) {
  rows <- continuous_table[continuous_table$metric == metric, , drop = FALSE]
  details <- lapply(seq_len(nrow(root_collection)), function(index) {
    selected <- rows[
      rows$cell_id == root_collection$cell_id[index] &
        rows$variant_id == root_collection$variant_id[index] &
        rows$moment_type == root_collection$moment_type[index],
      , drop = FALSE
    ]
    means <- selected$mean[is.finite(selected$mean)]
    mcses <- selected$mean_mcse[is.finite(selected$mean_mcse)]
    if (length(means) > 1L) {
      assert_true(
        max(abs(means - means[1L])) <= 1e-12 * max(1, abs(means)),
        sprintf("Root-length metric %s disagrees across K for %s.",
                metric, root_collection$cell_id[index])
      )
    }
    if (length(mcses) > 1L) {
      assert_true(
        max(abs(mcses - mcses[1L])) <= 1e-12 * max(1, abs(mcses)),
        sprintf("Root-length MCSE %s disagrees across K for %s.",
                metric, root_collection$cell_id[index])
      )
    }
    c(
      mean = if (length(means)) means[1L] else NA_real_,
      mean_mcse = if (length(mcses)) mcses[1L] else NA_real_
    )
  })
  as.data.frame(do.call(rbind, details), stringsAsFactors = FALSE)
}

wald_length_details <- root_collection_continuous_details(
  "wald_collection_total_length"
)
root_multiplier_length_details <- root_collection_continuous_details(
  "root_multiplier_total_length"
)

binomial_uncertainty <- function(successes, denominator, confidence = 0.95) {
  successes <- as.numeric(successes)
  denominator <- as.numeric(denominator)
  valid <- is.finite(successes) & is.finite(denominator) & denominator > 0
  estimate <- mcse <- lower <- upper <- rep(NA_real_, length(denominator))
  estimate[valid] <- successes[valid] / denominator[valid]
  mcse[valid] <- sqrt(
    estimate[valid] * (1 - estimate[valid]) / denominator[valid]
  )
  z <- stats::qnorm(1 - (1 - confidence) / 2)
  centre <- (estimate[valid] + z^2 / (2 * denominator[valid])) /
    (1 + z^2 / denominator[valid])
  half <- z * sqrt(
    estimate[valid] * (1 - estimate[valid]) / denominator[valid] +
      z^2 / (4 * denominator[valid]^2)
  ) / (1 + z^2 / denominator[valid])
  lower[valid] <- pmax(0, centre - half)
  upper[valid] <- pmin(1, centre + half)
  data.frame(
    mcse = mcse,
    wilson_lower = lower,
    wilson_upper = upper
  )
}

make_method_rows <- function(method, prefix, length_details) {
  report_count <- root_collection[[paste0(prefix, "_report_count")]]
  report_denominator <-
    root_collection[[paste0(prefix, "_applicable_replications")]]
  conditional_count <- root_collection[[paste0(
    prefix, "_conditional_coverage_count"
  )]]
  conditional_denominator <- root_collection[[paste0(
    prefix, "_conditional_coverage_denominator"
  )]]
  report_cover_count <- root_collection[[paste0(
    prefix, "_report_and_cover_count"
  )]]
  report_cover_denominator <- root_collection[[paste0(
    prefix, "_report_and_cover_denominator"
  )]]
  report_uncertainty <- binomial_uncertainty(
    report_count, report_denominator
  )
  conditional_uncertainty <- binomial_uncertainty(
    conditional_count, conditional_denominator
  )
  report_cover_uncertainty <- binomial_uncertainty(
    report_cover_count, report_cover_denominator
  )
  data.frame(
    root_collection[, root_identity_columns, drop = FALSE],
    interval_method = method,
    report_rate = root_collection[[paste0(prefix, "_report_rate")]],
    applicable_replications = report_denominator,
    report_count = report_count,
    report_rate_mcse = report_uncertainty$mcse,
    report_rate_wilson_lower = report_uncertainty$wilson_lower,
    report_rate_wilson_upper = report_uncertainty$wilson_upper,
    conditional_coverage =
      root_collection[[paste0(prefix, "_conditional_coverage")]],
    conditional_coverage_denominator = conditional_denominator,
    conditional_coverage_count = conditional_count,
    conditional_coverage_mcse = conditional_uncertainty$mcse,
    conditional_coverage_wilson_lower =
      conditional_uncertainty$wilson_lower,
    conditional_coverage_wilson_upper =
      conditional_uncertainty$wilson_upper,
    report_and_cover_rate =
      root_collection[[paste0(prefix, "_report_and_cover_rate")]],
    report_and_cover_denominator = report_cover_denominator,
    report_and_cover_count = report_cover_count,
    report_and_cover_mcse = report_cover_uncertainty$mcse,
    report_and_cover_wilson_lower =
      report_cover_uncertainty$wilson_lower,
    report_and_cover_wilson_upper =
      report_cover_uncertainty$wilson_upper,
    mean_total_length_conditional_on_report = length_details$mean,
    mean_total_length_mcse_conditional_on_report = length_details$mean_mcse,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}
table_03_root_inference <- rbind(
  make_method_rows("Wald", "wald", wald_length_details),
  make_method_rows(
    "Multiplier bootstrap", "root_multiplier",
    root_multiplier_length_details
  )
)
table_03_root_inference <- table_03_root_inference[order(
  table_03_root_inference$cell_class,
  table_03_root_inference$scenario_id,
  table_03_root_inference$n_x,
  table_03_root_inference$n_y,
  table_03_root_inference$moment_type,
  table_03_root_inference$interval_method
), , drop = FALSE]
row.names(table_03_root_inference) <- NULL
assert_true(
  nrow(table_03_root_inference) == 2L * nrow(root_collection) &&
    !anyDuplicated(table_03_root_inference[, c(
      "cell_id", "variant_id", "moment_type", "interval_method"
    )]),
  "Table 3 is not one row per deduplicated root collection and method."
)

root_accuracy_columns <- c(
  "study_id", "cell_id", "cell_class", "variant_id", "scenario_id",
  "n_x", "n_y", "moment_type", "root_id", "true_root",
  "applicable_successful_replications", "joint_isolation_replications",
  "finite_error_replications", "joint_isolation_report_rate",
  "bias_conditional_on_joint_isolation",
  "bias_mcse_conditional_on_joint_isolation",
  "rmse_conditional_on_joint_isolation",
  "rmse_mcse_conditional_on_joint_isolation",
  "error_sd_conditional_on_joint_isolation",
  "error_q05_conditional_on_joint_isolation",
  "error_q50_conditional_on_joint_isolation",
  "error_q95_conditional_on_joint_isolation",
  "source_interval_ids", "source_interval_count"
)
assert_true(
  all(root_accuracy_columns %in% names(root_specific)),
  "The K-deduplicated root-specific table is incomplete."
)
table_03b_root_accuracy <- root_specific[, root_accuracy_columns, drop = FALSE]
table_03b_root_accuracy <- table_03b_root_accuracy[order(
  table_03b_root_accuracy$cell_class,
  table_03b_root_accuracy$scenario_id,
  table_03b_root_accuracy$n_x,
  table_03b_root_accuracy$n_y,
  table_03b_root_accuracy$moment_type,
  table_03b_root_accuracy$true_root
), , drop = FALSE]
row.names(table_03b_root_accuracy) <- NULL

# Table 4: the six alternatives selected before the confirmatory Monte Carlo.
# Each selected DGP is configured on exactly one K, hence no K aggregation is
# needed.  Joint-root reporting fields are joined from the deduplicated table.
special_scenarios <- c("TANGENCY_SV_POWER", "TWO_ROOT_SV_POWER")
special_source <- summary_table[
  summary_table$scenario_id %in% special_scenarios &
    summary_table$moment_type == "absolute",
  , drop = FALSE
]
assert_true(
  !anyDuplicated(special_source[, c("cell_id", "moment_type")]),
  "A selected-power cell has more than one K publication row."
)
special_source$.article_order <- seq_len(nrow(special_source))
special_join_columns <- c(
  "cell_id", "variant_id", "moment_type",
  "root_collection_applicable_replications",
  "root_collection_joint_isolation_reports",
  "root_collection_joint_isolation_report_rate",
  "wald_report_rate", "wald_report_count",
  "wald_conditional_coverage", "wald_conditional_coverage_denominator",
  "wald_report_and_cover_rate", "wald_report_and_cover_count",
  "root_multiplier_report_rate", "root_multiplier_report_count",
  "root_multiplier_conditional_coverage",
  "root_multiplier_conditional_coverage_denominator",
  "root_multiplier_report_and_cover_rate",
  "root_multiplier_report_and_cover_count"
)
special_inference_columns <- setdiff(
  special_join_columns, c("cell_id", "variant_id", "moment_type")
)
# Drop the interval-specific copies before joining.  The columns in Table 4
# below therefore have only one possible source: the K-deduplicated table.
special_source <- special_source[, setdiff(
  names(special_source), special_inference_columns
), drop = FALSE]
special_join <- root_collection[, intersect(
  special_join_columns, names(root_collection)
), drop = FALSE]
special_source <- merge(
  special_source, special_join,
  by = c("cell_id", "variant_id", "moment_type"),
  all.x = TRUE, sort = FALSE
)
special_source <- special_source[order(special_source$.article_order), ,
                                 drop = FALSE]
special_source$n_eff <- with(
  special_source, as.numeric(n_x) * as.numeric(n_y) /
    (as.numeric(n_x) + as.numeric(n_y))
)
special_columns <- c(
  "study_id", "cell_id", "cell_class", "scenario_id", "n_x", "n_y",
  "n_eff", "moment_type", "interval_id", "expected_replications",
  "successful_replications", "continuum_band_coverage",
  "outer_set_truth_coverage", "reversal_certification_rate",
  "exact_count_certification_rate",
  "root_collection_applicable_replications",
  "root_collection_joint_isolation_reports",
  "root_collection_joint_isolation_report_rate",
  "wald_report_rate", "wald_report_count", "wald_conditional_coverage",
  "wald_conditional_coverage_denominator", "wald_report_and_cover_rate",
  "wald_report_and_cover_count", "root_multiplier_report_rate",
  "root_multiplier_report_count", "root_multiplier_conditional_coverage",
  "root_multiplier_conditional_coverage_denominator",
  "root_multiplier_report_and_cover_rate",
  "root_multiplier_report_and_cover_count", "mean_outer_set_length",
  "tangency_global_geometry_applicable_replications",
  "tangency_global_geometry_continuum_coverage_conditioned_replications",
  "tangency_global_geometry_conditioning_rate",
  "tangency_global_outer_nonempty_rate_conditional_on_continuum_band_coverage",
  "tangency_global_outer_contains_root_rate_conditional_on_continuum_band_coverage",
  "tangency_global_outer_total_length_finite_count",
  "tangency_global_outer_total_length_missing_count",
  "tangency_global_outer_total_length_positive_infinite_count",
  "tangency_global_outer_hausdorff_finite_count",
  "tangency_global_outer_hausdorff_missing_count",
  "tangency_global_outer_hausdorff_positive_infinite_count",
  "mean_tangency_global_outer_total_length_conditional_on_continuum_band_coverage",
  "mean_tangency_global_outer_hausdorff_conditional_on_continuum_band_coverage",
  "mean_tangency_global_scaled_outer_total_length_conditional_on_continuum_band_coverage",
  "mean_tangency_global_scaled_outer_hausdorff_conditional_on_continuum_band_coverage"
)
assert_true(
  all(special_columns %in% names(special_source)),
  sprintf(
    "The selected-power table lacks publication fields: %s.",
    paste(setdiff(special_columns, names(special_source)), collapse = ", ")
  )
)
table_04_special_power <- select_existing(special_source, special_columns)
table_04_special_power <- table_04_special_power[order(
  match(table_04_special_power$scenario_id, special_scenarios),
  table_04_special_power$n_x, table_04_special_power$n_y
), , drop = FALSE]
row.names(table_04_special_power) <- NULL
tangency_rows_index <- table_04_special_power$scenario_id ==
  "TANGENCY_SV_POWER"
for (field in c(
    "reversal_certification_rate", "exact_count_certification_rate",
    "root_collection_joint_isolation_report_rate", "wald_report_rate",
    "wald_conditional_coverage", "wald_report_and_cover_rate",
    "root_multiplier_report_rate", "root_multiplier_conditional_coverage",
    "root_multiplier_report_and_cover_rate")) {
  table_04_special_power[[field]][tangency_rows_index] <- NA_real_
}

continuous_details_for_rows <- function(data, metric) {
  rows <- continuous_table[continuous_table$metric == metric, , drop = FALSE]
  key_columns <- c("cell_id", "moment_type", "interval_id")
  keys <- do.call(paste, c(rows[key_columns], list(sep = "::")))
  targets <- do.call(paste, c(data[key_columns], list(sep = "::")))
  assert_true(!anyDuplicated(keys),
              sprintf("Continuous publication metric %s is duplicated.", metric))
  positions <- match(targets, keys)
  data.frame(
    mean = rows$mean[positions],
    mean_mcse = rows$mean_mcse[positions]
  )
}

tangency_metric_map <- c(
  mean_tangency_global_outer_total_length_conditional_on_continuum_band_coverage =
    "tangency_global_outer_total_length_conditional_on_continuum_band_coverage",
  mean_tangency_global_outer_hausdorff_conditional_on_continuum_band_coverage =
    "tangency_global_outer_hausdorff_conditional_on_continuum_band_coverage",
  mean_tangency_global_scaled_outer_total_length_conditional_on_continuum_band_coverage =
    "tangency_global_scaled_outer_total_length_conditional_on_continuum_band_coverage",
  mean_tangency_global_scaled_outer_hausdorff_conditional_on_continuum_band_coverage =
    "tangency_global_scaled_outer_hausdorff_conditional_on_continuum_band_coverage"
)
for (output_name in names(tangency_metric_map)) {
  details <- continuous_details_for_rows(
    table_04_special_power, tangency_metric_map[[output_name]]
  )
  comparable <- is.finite(table_04_special_power[[output_name]]) &
    is.finite(details$mean)
  assert_true(
    all(abs(table_04_special_power[[output_name]][comparable] -
              details$mean[comparable]) <= 1e-12),
    sprintf("Wide and long summaries disagree for %s.", output_name)
  )
  table_04_special_power[[paste0(output_name, "_mcse")]] <-
    details$mean_mcse
}

# Supplementary tables retain all interval-specific results and the complete
# stress class.  No nominal-coverage pass/fail rule is imposed on stress rows.
supplement_all_intervals <- summary_table
supplement_all_intervals <- supplement_all_intervals[order(
  supplement_all_intervals$cell_class,
  supplement_all_intervals$scenario_id,
  supplement_all_intervals$n_x,
  supplement_all_intervals$n_y,
  supplement_all_intervals$moment_type,
  supplement_all_intervals$interval_id
), , drop = FALSE]
row.names(supplement_all_intervals) <- NULL
stress_source <- supplement_all_intervals[
  supplement_all_intervals$cell_class == "stress", , drop = FALSE
]
stress_columns <- c(
  "study_id", "cell_id", "cell_class", "scenario_id", "n_x", "n_y",
  "moment_type", "interval_id", "p_min", "p_max", "bootstrap_reps",
  "expected_replications", "successful_replications",
  "replication_failure_rate", "continuum_band_coverage",
  "continuum_band_coverage_wilson_lower",
  "continuum_band_coverage_wilson_upper",
  "any_primary_false_report_rate", "familywise_primary_error_rate",
  "familywise_primary_error_wilson_lower",
  "familywise_primary_error_wilson_upper",
  "reversal_certification_rate", "exact_count_certification_rate",
  "continuum_no_root_certification_rate", "outer_set_truth_coverage",
  "root_isolation_report_rate", "root_bias_conditional_on_isolation",
  "root_rmse_conditional_on_isolation", "wald_report_rate",
  "wald_conditional_coverage", "wald_report_and_cover_rate",
  "root_multiplier_report_rate", "root_multiplier_conditional_coverage",
  "root_multiplier_report_and_cover_rate", "guard_failure_rate",
  "enclosure_limit_rate", "enclosure_variance_limit_rate",
  "mean_maximum_half_width", "mean_outer_set_length",
  "mean_y_ess_ratio_upper", "mean_y_max_share_upper",
  "mean_elapsed_seconds_per_replication"
)
assert_true(all(stress_columns %in% names(stress_source)), sprintf(
  "The stress supplement lacks fields: %s.",
  paste(setdiff(stress_columns, names(stress_source)), collapse = ", ")
))
supplement_stress <- stress_source[, stress_columns, drop = FALSE]

continuous_means <- function(metric) {
  rows <- continuous_table[continuous_table$metric == metric, , drop = FALSE]
  key_columns <- c("cell_id", "variant_id", "moment_type", "interval_id")
  keys <- do.call(paste, c(rows[key_columns], list(sep = "::")))
  assert_true(!anyDuplicated(keys),
              sprintf("Continuous metric %s is duplicated.", metric))
  target_keys <- do.call(paste, c(summary_table[key_columns], list(sep = "::")))
  rows$mean[match(target_keys, keys)]
}

failure_columns <- c(
  "cell_id", "expected_units", "successful_units",
  "numerical_failure_units", "numerical_failure_rate"
)
assert_true(all(failure_columns %in% names(failure_table)),
            "The numerical-failure cell summary is incomplete.")
assert_true(!anyDuplicated(failure_table$cell_id),
            "The numerical-failure summary duplicates a cell.")
numerical_source <- summary_table
numerical_source$.article_order <- seq_len(nrow(numerical_source))
numerical_source <- merge(
  numerical_source,
  failure_table[, failure_columns, drop = FALSE],
  by = "cell_id", all.x = TRUE, sort = FALSE
)
numerical_source <- numerical_source[order(numerical_source$.article_order), ,
                                     drop = FALSE]
for (metric in c(
    "enclosure_nodes_used", "enclosure_maximum_level",
    "enclosure_unresolved_cells", "elapsed_seconds", "peak_r_heap_mb")) {
  numerical_source[[paste0("mean_", metric)]] <- continuous_means(metric)
}
numerical_columns <- c(
  "study_id", "cell_id", "cell_class", "scenario_id", "n_x", "n_y",
  "moment_type", "interval_id", "expected_replications",
  "successful_replications", "replication_failure_rate", "expected_units",
  "successful_units", "numerical_failure_units", "numerical_failure_rate",
  "guard_failure_rate", "truth_score_limit_rate", "enclosure_limit_rate",
  "enclosure_depth_limit_rate", "enclosure_node_limit_rate",
  "enclosure_variance_limit_rate", "enclosure_statistical_limit_rate",
  "mean_enclosure_limit_hit_cells", "mean_enclosure_limit_hit_total_width",
  "maximum_enclosure_limit_hit_cell_width", "mean_enclosure_nodes_used",
  "mean_enclosure_maximum_level", "mean_enclosure_unresolved_cells",
  "mean_elapsed_seconds", "mean_peak_r_heap_mb"
)
supplement_numerical_audit <- select_existing(numerical_source,
                                              numerical_columns)

# Enforce the frozen publication contract for the confirmatory main.  The
# smoke intentionally has fewer and differently classified rows, but must
# still exercise every builder and both graphics devices.
article_contract <- config$article_output_contract %||% list()
if (length(article_contract)) {
  observed_article_contract <- list(
    design_rows = nrow(table_01_design),
    benchmark_editorial_rows = nrow(table_02_benchmark),
    # The frozen contract counts K-deduplicated source collections.  The
    # editorial table deliberately expands each source to two method rows.
    root_collection_rows = nrow(root_collection),
    root_specific_rows = nrow(table_03b_root_accuracy),
    special_power_rows = nrow(table_04_special_power),
    stress_scientific_groups = nrow(supplement_stress)
  )
  unknown <- setdiff(names(article_contract), names(observed_article_contract))
  assert_true(!length(unknown), sprintf(
    "Unknown article-output contract fields: %s.", paste(unknown, collapse = ", ")
  ))
  mismatches <- names(article_contract)[vapply(names(article_contract), function(name) {
    as.integer(article_contract[[name]]) !=
      as.integer(observed_article_contract[[name]])
  }, logical(1L))]
  assert_true(!length(mismatches), sprintf(
    "Article-output contract mismatch: %s.",
    paste(vapply(mismatches, function(name) sprintf(
      "%s observed=%d expected=%d", name,
      observed_article_contract[[name]], article_contract[[name]]
    ), character(1L)), collapse = "; ")
  ))
}

dir.create(article_directory, recursive = TRUE, showWarnings = FALSE)

write_article_csv <- function(data, filename) {
  atomic_write_csv(data, file.path(article_directory, filename))
}

latex_escape <- function(values) {
  values <- as.character(values)
  values[is.na(values)] <- ""
  replacements <- c(
    "\\" = "\\textbackslash{}", "&" = "\\&", "%" = "\\%",
    "$" = "\\$", "#" = "\\#", "_" = "\\_", "{" = "\\{",
    "}" = "\\}", "~" = "\\textasciitilde{}", "^" = "\\textasciicircum{}"
  )
  for (pattern in names(replacements)) {
    values <- gsub(pattern, replacements[[pattern]], values,
                   fixed = TRUE, useBytes = TRUE)
  }
  values
}

write_article_tex <- function(data, path) {
  display <- data
  for (name in names(display)) {
    if (is.numeric(display[[name]])) {
      display[[name]] <- format(display[[name]], trim = TRUE, scientific = FALSE)
    }
    display[[name]] <- latex_escape(display[[name]])
  }
  column_count <- ncol(display)
  lines <- c(
    "% Generated from the frozen simulation configuration.",
    sprintf("\\begin{tabular}{%s}", paste(rep("l", column_count), collapse = "")),
    "\\hline",
    paste(latex_escape(names(display)), collapse = " & "),
    "\\\\",
    "\\hline"
  )
  if (nrow(display)) {
    body <- apply(display, 1L, function(row) {
      paste0(paste(row, collapse = " & "), " \\\\")
    })
    lines <- c(lines, body)
  }
  lines <- c(lines, "\\hline", "\\end{tabular}")
  atomic_write_lines(lines, path)
}

format_probability <- function(estimate, successes = NA_real_,
                               denominator = NA_real_, lower = NA_real_,
                               upper = NA_real_) {
  if (!length(estimate)) return(character(0L))
  size <- max(length(estimate), length(successes), length(denominator),
              length(lower), length(upper))
  estimate <- rep_len(as.numeric(estimate), size)
  successes <- rep_len(as.numeric(successes), size)
  denominator <- rep_len(as.numeric(denominator), size)
  lower <- rep_len(as.numeric(lower), size)
  upper <- rep_len(as.numeric(upper), size)
  vapply(seq_len(size), function(index) {
    if (!is.finite(estimate[index])) return("--")
    value <- sprintf("%.3f", estimate[index])
    if (is.finite(successes[index]) && is.finite(denominator[index])) {
      value <- paste0(value, sprintf(
        " (%d/%d)", as.integer(successes[index]),
        as.integer(denominator[index])
      ))
    }
    if (is.finite(lower[index]) && is.finite(upper[index])) {
      value <- paste0(value, sprintf(" [%.3f, %.3f]",
                                     lower[index], upper[index]))
    }
    value
  }, character(1L))
}

format_number_mcse <- function(estimate, mcse) {
  estimate <- as.numeric(estimate)
  mcse <- rep_len(as.numeric(mcse), length(estimate))
  vapply(seq_along(estimate), function(index) {
    if (!is.finite(estimate[index])) return("--")
    if (is.finite(mcse[index])) {
      sprintf("%.4f (MCSE %.4f)", estimate[index], mcse[index])
    } else sprintf("%.4f", estimate[index])
  }, character(1L))
}

sample_labels <- function(data) ifelse(
  data$n_x == data$n_y, as.character(data$n_x),
  paste0(data$n_x, "/", data$n_y)
)

format_k_diagnostic <- function(data, interval_id) {
  retained <- data[[paste0(interval_id, "_outer_set_retained_proportion")]]
  no_root <- data[[paste0(interval_id, "_no_root_certification_rate")]]
  exact <- data[[paste0(interval_id, "_exact_count_certification_rate")]]
  vapply(seq_len(nrow(data)), function(index) {
    pieces <- character(0L)
    if (is.finite(retained[index])) {
      pieces <- c(pieces, sprintf("retained %.3f", retained[index]))
    }
    if (is.finite(no_root[index])) {
      pieces <- c(pieces, sprintf("no-root %.3f", no_root[index]))
    }
    if (is.finite(exact[index])) {
      pieces <- c(pieces, sprintf("exact %.3f", exact[index]))
    }
    if (length(pieces)) paste(pieces, collapse = "; ") else "--"
  }, character(1L))
}

# CSV files retain every audit column.  The LaTeX counterparts are deliberately
# compact editorial views so the manuscript need not expose the full machine
# contract in the main text.
table_01_tex <- data.frame(
  Design = table_01_design$scenario_label,
  Class = table_01_design$scientific_class,
  `Sample sizes` = table_01_design$sample_sizes_nx_over_ny,
  Moments = table_01_design$moment_types,
  `Order intervals` = table_01_design$order_intervals,
  `Population roots` = table_01_design$population_root_truth,
  `Structural budget` = table_01_design$structural_multiplicity_budget,
  `Regular refinement` = table_01_design$root_inference,
  check.names = FALSE, stringsAsFactors = FALSE
)
table_02_tex <- data.frame(
  Design = table_02_benchmark$scenario_id,
  n = sample_labels(table_02_benchmark),
  Moment = table_02_benchmark$moment_type,
  `Simultaneous band coverage` = format_probability(
    table_02_benchmark$simultaneous_continuum_band_coverage,
    table_02_benchmark$simultaneous_continuum_band_successes,
    table_02_benchmark$simultaneous_continuum_band_denominator,
    table_02_benchmark$simultaneous_continuum_band_coverage_wilson_lower,
    table_02_benchmark$simultaneous_continuum_band_coverage_wilson_upper
  ),
  `Familywise primary error` = format_probability(
    table_02_benchmark$familywise_primary_error_rate,
    table_02_benchmark$familywise_primary_error_count,
    table_02_benchmark$familywise_primary_error_denominator,
    table_02_benchmark$familywise_primary_error_wilson_lower,
    table_02_benchmark$familywise_primary_error_wilson_upper
  ),
  K1 = format_k_diagnostic(table_02_benchmark, "K1"),
  K2 = format_k_diagnostic(table_02_benchmark, "K2"),
  K3 = format_k_diagnostic(table_02_benchmark, "K3"),
  check.names = FALSE, stringsAsFactors = FALSE
)
table_03_tex <- data.frame(
  Design = table_03_root_inference$scenario_id,
  n = sample_labels(table_03_root_inference),
  Moment = table_03_root_inference$moment_type,
  Method = table_03_root_inference$interval_method,
  `Joint isolation` = sprintf(
    "%.3f (%d/%d)",
    table_03_root_inference$root_collection_joint_isolation_report_rate,
    table_03_root_inference$root_collection_joint_isolation_reports,
    table_03_root_inference$root_collection_applicable_replications
  ),
  `Report rate` = format_probability(
    table_03_root_inference$report_rate,
    table_03_root_inference$report_count,
    table_03_root_inference$applicable_replications,
    table_03_root_inference$report_rate_wilson_lower,
    table_03_root_inference$report_rate_wilson_upper
  ),
  `Coverage given report` = format_probability(
    table_03_root_inference$conditional_coverage,
    table_03_root_inference$conditional_coverage_count,
    table_03_root_inference$conditional_coverage_denominator,
    table_03_root_inference$conditional_coverage_wilson_lower,
    table_03_root_inference$conditional_coverage_wilson_upper
  ),
  `Report and cover` = format_probability(
    table_03_root_inference$report_and_cover_rate,
    table_03_root_inference$report_and_cover_count,
    table_03_root_inference$report_and_cover_denominator,
    table_03_root_inference$report_and_cover_wilson_lower,
    table_03_root_inference$report_and_cover_wilson_upper
  ),
  `Length given report` = format_number_mcse(
    table_03_root_inference$mean_total_length_conditional_on_report,
    table_03_root_inference$mean_total_length_mcse_conditional_on_report
  ),
  check.names = FALSE, stringsAsFactors = FALSE
)
table_03b_tex <- data.frame(
  Design = table_03b_root_accuracy$scenario_id,
  n = sample_labels(table_03b_root_accuracy),
  Moment = table_03b_root_accuracy$moment_type,
  Root = table_03b_root_accuracy$true_root,
  `Joint-isolation denominator` =
    table_03b_root_accuracy$joint_isolation_replications,
  `Conditional bias` = format_number_mcse(
    table_03b_root_accuracy$bias_conditional_on_joint_isolation,
    table_03b_root_accuracy$bias_mcse_conditional_on_joint_isolation
  ),
  `Conditional RMSE` = format_number_mcse(
    table_03b_root_accuracy$rmse_conditional_on_joint_isolation,
    table_03b_root_accuracy$rmse_mcse_conditional_on_joint_isolation
  ),
  check.names = FALSE, stringsAsFactors = FALSE
)
table_04_tex <- data.frame(
  Design = table_04_special_power$scenario_id,
  n = sample_labels(table_04_special_power),
  `Band coverage` = format_probability(
    table_04_special_power$continuum_band_coverage
  ),
  `Outer-set root coverage` = format_probability(
    table_04_special_power$outer_set_truth_coverage
  ),
  `Joint isolation` = format_probability(
    table_04_special_power$root_collection_joint_isolation_report_rate
  ),
  `Exact certification` = format_probability(
    table_04_special_power$exact_count_certification_rate
  ),
  `Multiplier report and cover` = format_probability(
    table_04_special_power$root_multiplier_report_and_cover_rate,
    table_04_special_power$root_multiplier_report_and_cover_count,
    table_04_special_power$root_collection_applicable_replications
  ),
  `Global outer length` = format_number_mcse(
    table_04_special_power$mean_tangency_global_outer_total_length_conditional_on_continuum_band_coverage,
    table_04_special_power$mean_tangency_global_outer_total_length_conditional_on_continuum_band_coverage_mcse
  ),
  `Global Hausdorff` = format_number_mcse(
    table_04_special_power$mean_tangency_global_outer_hausdorff_conditional_on_continuum_band_coverage,
    table_04_special_power$mean_tangency_global_outer_hausdorff_conditional_on_continuum_band_coverage_mcse
  ),
  check.names = FALSE, stringsAsFactors = FALSE
)
supplement_stress_tex <- data.frame(
  Design = supplement_stress$scenario_id,
  n = sample_labels(supplement_stress),
  Moment = supplement_stress$moment_type,
  K = supplement_stress$interval_id,
  `Band coverage` = format_probability(supplement_stress$continuum_band_coverage),
  `Familywise error` = format_probability(
    supplement_stress$familywise_primary_error_rate
  ),
  `No-root certificate` = format_probability(
    supplement_stress$continuum_no_root_certification_rate
  ),
  `Exact certificate` = format_probability(
    supplement_stress$exact_count_certification_rate
  ),
  `Outer-set root coverage` = format_probability(
    supplement_stress$outer_set_truth_coverage
  ),
  check.names = FALSE, stringsAsFactors = FALSE
)

write_article_csv(table_01_design, "table_01_design.csv")
write_article_tex(table_01_tex,
                  file.path(article_directory, "table_01_design.tex"))
write_article_csv(table_02_benchmark, "table_02_benchmark.csv")
write_article_tex(table_02_tex,
                  file.path(article_directory, "table_02_benchmark.tex"))
write_article_csv(table_03_root_inference, "table_03_root_inference.csv")
write_article_tex(table_03_tex,
                  file.path(article_directory, "table_03_root_inference.tex"))
write_article_csv(table_03b_root_accuracy, "table_03b_root_accuracy.csv")
write_article_tex(table_03b_tex,
                  file.path(article_directory, "table_03b_root_accuracy.tex"))
write_article_csv(table_04_special_power, "table_04_special_power.csv")
write_article_tex(table_04_tex,
                  file.path(article_directory, "table_04_special_power.tex"))
write_article_csv(supplement_stress, "supplement_stress.csv")
write_article_tex(supplement_stress_tex,
                  file.path(article_directory, "supplement_stress.tex"))
write_article_csv(supplement_all_intervals, "supplement_all_intervals.csv")
write_article_csv(supplement_numerical_audit,
                  "supplement_numerical_audit.csv")

safe_plot_range <- function(values, log_scale = FALSE) {
  values <- as.numeric(values)
  values <- values[is.finite(values)]
  if (log_scale) values <- values[values > 0]
  if (!length(values)) return(if (log_scale) c(1, 2) else c(0, 1))
  range_values <- range(values)
  if (range_values[1L] == range_values[2L]) {
    if (log_scale) {
      range_values * c(0.8, 1.25)
    } else {
      delta <- max(0.05, abs(range_values[1L]) * 0.08)
      range_values + c(-delta, delta)
    }
  } else range_values
}

plot_lines_or_message <- function(x, series, labels, colours, ylab,
                                  log_axes = "", ylim = NULL) {
  keep_series <- vapply(series, function(values) {
    any(is.finite(as.numeric(values)))
  }, logical(1L))
  series <- series[keep_series]
  labels <- labels[keep_series]
  colours <- colours[keep_series]
  if (!length(series) || !any(is.finite(x))) {
    plot.new()
    title(ylab = ylab)
    text(0.5, 0.5, "Metric not applicable in this smoke/design")
    return(invisible(NULL))
  }
  y_values <- unlist(series, use.names = FALSE)
  if (is.null(ylim)) ylim <- safe_plot_range(y_values, grepl("y", log_axes))
  xlim <- safe_plot_range(x, grepl("x", log_axes))
  point_symbols <- seq.int(16L, length.out = length(series))
  plot(
    x, series[[1L]], type = "b", pch = point_symbols[1L], lwd = 2,
    col = colours[1L], xlab = "Effective sample size",
    ylab = ylab, xlim = xlim, ylim = ylim, log = log_axes
  )
  if (length(series) > 1L) {
    for (index in 2:length(series)) {
      lines(x, series[[index]], type = "b", pch = point_symbols[index],
            lwd = 2, col = colours[index])
    }
  }
  grid(col = "grey88")
  legend("bottomright", legend = labels, col = colours, lwd = 2,
         pch = point_symbols, bty = "n",
         cex = 0.85)
  invisible(NULL)
}

render_atomic_figure <- function(filename, renderer) {
  final_path <- file.path(article_directory, filename)
  extension <- tolower(tools::file_ext(filename))
  temporary <- tempfile(
    pattern = paste0(".", tools::file_path_sans_ext(basename(filename)), "."),
    tmpdir = article_directory, fileext = paste0(".", extension)
  )
  on.exit(unlink(temporary), add = TRUE)
  if (extension == "pdf") {
    grDevices::pdf(temporary, width = 9.0, height = 7.0,
                   useDingbats = FALSE, onefile = TRUE)
  } else if (extension == "png") {
    png_arguments <- list(filename = temporary, width = 1800L, height = 1400L,
                          res = 180L)
    if (isTRUE(capabilities("cairo"))) png_arguments$type <- "cairo-png"
    do.call(grDevices::png, png_arguments)
  } else {
    stopf("Unsupported article figure extension: %s", extension)
  }
  device_open <- TRUE
  tryCatch(
    renderer(),
    finally = {
      if (device_open) grDevices::dev.off()
      device_open <- FALSE
    }
  )
  assert_true(file.exists(temporary) && file.info(temporary)$size > 0L,
              sprintf("Figure %s is empty.", filename))
  atomic_publish(temporary, final_path, overwrite = TRUE)
}

tangency_rows <- table_04_special_power[
  table_04_special_power$scenario_id == "TANGENCY_SV_POWER", , drop = FALSE
]
tangency_rows <- tangency_rows[order(tangency_rows$n_eff), , drop = FALSE]

plot_tangency_metric <- function(x, y, mcse, title_text, y_label,
                                 scaled = FALSE) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  mcse <- as.numeric(mcse)
  keep <- is.finite(x) & x > 0 & is.finite(y) & y > 0
  if (!any(keep)) {
    plot.new()
    title(main = title_text, ylab = y_label)
    text(0.5, 0.5, "Metric not applicable in this smoke/design")
    return(invisible(NULL))
  }
  x <- x[keep]
  y <- y[keep]
  mcse <- mcse[keep]
  mcse[!is.finite(mcse)] <- 0
  lower <- pmax(y - 1.96 * mcse, min(y) / 20)
  upper <- y + 1.96 * mcse
  plot(
    x, y, type = "b", pch = 16L, lwd = 2, col = "#1b6ca8",
    log = if (scaled) "x" else "xy",
    xlim = safe_plot_range(x, TRUE),
    ylim = safe_plot_range(c(lower, upper), !scaled),
    xlab = "Effective sample size", ylab = y_label,
    main = title_text
  )
  arrows(x, lower, x, upper, angle = 90, code = 3,
         length = 0.035, col = "#1b6ca8")
  if (!scaled && length(x) >= 2L) {
    reference <- y[1L] * (x / x[1L])^(-1 / 4)
    lines(x, reference, lty = 2L, lwd = 1.6, col = "#d1495b")
    legend("topright", c("Monte Carlo mean", "slope -1/4"),
           col = c("#1b6ca8", "#d1495b"), lwd = c(2, 1.6),
           lty = c(1, 2), pch = c(16, NA), bty = "n", cex = 0.78)
  }
  grid(col = "grey88")
  invisible(NULL)
}

render_tangency <- function() {
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)
  par(mfrow = c(2, 2), mar = c(4.1, 4.2, 2.5, 1.0),
      oma = c(0, 0, 2.0, 0), las = 1)
  plot_tangency_metric(
    tangency_rows$n_eff,
    tangency_rows$mean_tangency_global_outer_total_length_conditional_on_continuum_band_coverage,
    tangency_rows$mean_tangency_global_outer_total_length_conditional_on_continuum_band_coverage_mcse,
    "Outer-set length", "Conditional mean"
  )
  plot_tangency_metric(
    tangency_rows$n_eff,
    tangency_rows$mean_tangency_global_outer_hausdorff_conditional_on_continuum_band_coverage,
    tangency_rows$mean_tangency_global_outer_hausdorff_conditional_on_continuum_band_coverage_mcse,
    "Hausdorff distance", "Conditional mean"
  )
  plot_tangency_metric(
    tangency_rows$n_eff,
    tangency_rows$mean_tangency_global_scaled_outer_total_length_conditional_on_continuum_band_coverage,
    tangency_rows$mean_tangency_global_scaled_outer_total_length_conditional_on_continuum_band_coverage_mcse,
    "Fourth-root-scaled length", expression(N[eff]^{1/4} * L),
    scaled = TRUE
  )
  plot_tangency_metric(
    tangency_rows$n_eff,
    tangency_rows$mean_tangency_global_scaled_outer_hausdorff_conditional_on_continuum_band_coverage,
    tangency_rows$mean_tangency_global_scaled_outer_hausdorff_conditional_on_continuum_band_coverage_mcse,
    "Fourth-root-scaled Hausdorff", expression(N[eff]^{1/4} * d[H]),
    scaled = TRUE
  )
  mtext("Conditioned on simultaneous continuum-band coverage",
        side = 3, outer = TRUE, line = 0.4, cex = 0.88)
}
render_atomic_figure("figure_01_tangency_contraction.pdf", render_tangency)
render_atomic_figure("figure_01_tangency_contraction.png", render_tangency)

two_root_rows <- table_04_special_power[
  table_04_special_power$scenario_id == "TWO_ROOT_SV_POWER", , drop = FALSE
]
two_root_rows <- two_root_rows[order(two_root_rows$n_eff), , drop = FALSE]
render_two_root <- function() {
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)
  par(mfrow = c(1, 2), mar = c(4.2, 4.3, 2.4, 1.0),
      oma = c(0, 0, 1.8, 0), las = 1)
  plot_lines_or_message(
    two_root_rows$n_eff,
    list(
      two_root_rows$root_collection_joint_isolation_report_rate,
      two_root_rows$exact_count_certification_rate,
      two_root_rows$root_multiplier_report_and_cover_rate
    ),
    c("Joint isolation", "Exact-count certificate", "Report and cover"),
    c("#1b6ca8", "#238b45", "#8c5da5"),
    "Probability", log_axes = "x", ylim = c(0, 1)
  )
  title("Two roots: structural conclusions")
  plot_lines_or_message(
    two_root_rows$n_eff,
    list(
      two_root_rows$wald_conditional_coverage,
      two_root_rows$root_multiplier_conditional_coverage
    ),
    c("Wald: coverage | report", "Multiplier: coverage | report"),
    c("#1b6ca8", "#d1495b"),
    "Probability", log_axes = "x", ylim = c(0, 1)
  )
  title("Conditional coverage by method")
  mtext("Prespecified isolating intervals; no K-copy averaging",
        side = 3, outer = TRUE, line = 0.25, cex = 0.82)
}
render_atomic_figure("figure_02_two_root_power.pdf", render_two_root)
render_atomic_figure("figure_02_two_root_power.png", render_two_root)

# The manifest is generated last from a fixed allowlist and covers every
# required publication artifact other than itself.  A stale file from an
# earlier smoke can therefore never enter the confirmatory manifest.  Base R's
# MD5 is portable across the supported compute environments.
expected_artifact_files <- c(
  "table_01_design.csv", "table_01_design.tex",
  "table_02_benchmark.csv", "table_02_benchmark.tex",
  "table_03_root_inference.csv", "table_03_root_inference.tex",
  "table_03b_root_accuracy.csv", "table_03b_root_accuracy.tex",
  "table_04_special_power.csv", "table_04_special_power.tex",
  "figure_01_tangency_contraction.pdf",
  "figure_01_tangency_contraction.png",
  "figure_02_two_root_power.pdf", "figure_02_two_root_power.png",
  "supplement_stress.csv", "supplement_stress.tex",
  "supplement_all_intervals.csv", "supplement_numerical_audit.csv"
)
artifact_paths <- file.path(article_directory, expected_artifact_files)
assert_true(
  all(file.exists(artifact_paths)) && all(file.info(artifact_paths)$size > 0L),
  "At least one required article output is absent or empty."
)
artifact_info <- file.info(artifact_paths)
manifest <- data.frame(
  study_id = config$study_id,
  article_output_profile = config$article_output_profile %||% "unspecified",
  generated_utc = timestamp_utc(),
  file = basename(artifact_paths),
  bytes = as.numeric(artifact_info$size),
  md5 = unname(tools::md5sum(artifact_paths)),
  stringsAsFactors = FALSE
)
assert_true(nrow(manifest) == length(expected_artifact_files) &&
              all(manifest$bytes > 0L),
            "The article artifact set is incomplete or contains empty files.")
atomic_write_csv(manifest, file.path(article_directory, "manifest.csv"))
atomic_write_lines(
  c(
    paste0("study_id=", config$study_id),
    paste0("article_output_profile=",
           config$article_output_profile %||% "unspecified"),
    paste0("validation_complete=", isTRUE(validation$complete)),
    paste0("generated_utc=", manifest$generated_utc[1L]),
    paste0("artifact_count=", nrow(manifest)),
    vapply(seq_len(nrow(manifest)), function(index) sprintf(
      "%s bytes=%s md5=%s", manifest$file[index], manifest$bytes[index],
      manifest$md5[index]
    ), character(1L))
  ),
  file.path(article_directory, "article_outputs_manifest.txt")
)

log_message(
  "Wrote publication outputs for study=", config$study_id,
  ": design_rows=", nrow(table_01_design),
  ", benchmark_rows=", nrow(table_02_benchmark),
  ", root_collection_source_rows=", nrow(root_collection),
  ", root_method_rows=", nrow(table_03_root_inference),
  ", root_specific_rows=", nrow(table_03b_root_accuracy),
  ", special_power_rows=", nrow(table_04_special_power),
  ", stress_rows=", nrow(supplement_stress),
  ", directory=", article_directory
)
