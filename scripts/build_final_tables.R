#!/usr/bin/env Rscript

arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 3L) {
  stop(
    paste(
      "Usage: Rscript scripts/build_final_tables.R",
      "SUMMARY_DIRECTORY ARTICLE_DIRECTORY CONFIG"
    ),
    call. = FALSE
  )
}

options(digits = 17)
summary_directory <- normalizePath(arguments[[1L]], mustWork = TRUE)
article_directory <- normalizePath(arguments[[2L]], mustWork = TRUE)
config_path <- normalizePath(arguments[[3L]], mustWork = TRUE)
config <- dget(config_path)

read_csv <- function(directory, filename) {
  read.csv(
    file.path(directory, filename), stringsAsFactors = FALSE,
    check.names = FALSE, na.strings = c("", "NA")
  )
}
format_three <- function(value) sprintf("%.3f", as.numeric(value))
write_table <- function(value, filename) {
  encode_field <- function(field) {
    if (length(field) == 0L || is.na(field)) return("")
    field <- as.character(field)
    escaped <- gsub('"', '""', field, fixed = TRUE)
    if (grepl('[,"\\r\\n]', field)) paste0('"', escaped, '"') else escaped
  }
  rows <- vapply(seq_len(nrow(value)), function(index) {
    paste(vapply(value[index, , drop = FALSE], encode_field, character(1L)),
          collapse = ",")
  }, character(1L))
  writeLines(
    c(paste(names(value), collapse = ","), rows),
    file.path(article_directory, filename), useBytes = TRUE
  )
}

cell_classes <- table(vapply(config$cells, `[[`, character(1L), "cell_class"))
stopifnot(
  length(config$cells) == 54L,
  identical(
    as.integer(cell_classes[c("benchmark", "power", "stress")]),
    c(10L, 18L, 26L)
  ),
  identical(as.integer(config$mc_reps_per_cell), 2032L),
  identical(as.integer(config$bootstrap_reps), 999L),
  identical(names(config$order_intervals), c("K1", "K2", "K3"))
)

table_01 <- data.frame(
  role = c(
    "Benchmark", "Benchmark", "Benchmark", "Power", "Power", "Power",
    "Power", "Challenging"
  ),
  mechanism = c(
    "Lognormal equality", "Lognormal no root", "Mixed sign",
    "Lognormal upcrossing", "Rare-tail upcrossing", "Selected tangency",
    "Selected two-root design",
    "High variance, sparse tails, Pareto, weak separation"
  ),
  population_feature = c(
    "identical absolute-moment curves",
    "no positive root; budget 0",
    "absolute equality; signed root 1 (up); budget 1",
    "root 4/3 (up); moderate or weak signal; budget 1",
    "root 1 (up); mixture probability 0.05; budget 1",
    "root 0.6, multiplicity 2, no sign change; budget 2",
    "roots 0.75 (down) and 1.75 (up); budget 2",
    "equality, no root, one root, tangency, or two roots"
  ),
  sample_sizes_nx_equals_ny = c(
    "50, 100, 250, 500", "50, 100, 250, 500", "100, 250",
    "100, 250, 500, 1000", "100, 250, 500, 1000",
    "2000, 4000, 8000", "2000, 4000, 8000",
    "100--1000; some unequal sizes"
  ),
  interval_numbers = c("1--3", "1--3", "1--3", "1--3", "1--3", "1", "3", "1--3"),
  stringsAsFactors = FALSE
)

summary <- read_csv(summary_directory, "summary.csv")
stopifnot(nrow(summary) == 140L)
binary <- read_csv(summary_directory, "binary_metrics.csv")
stopifnot(nrow(binary) == 14980L)
benchmark_spec <- data.frame(
  scenario_id = c("LN_EQUAL", "LN_NO_ROOT", "MIXED_SIGN"),
  design = c("Lognormal equality", "Lognormal no root", "Mixed sign, joint"),
  stringsAsFactors = FALSE
)
table_02_rows <- list()
counter <- 0L
for (spec_index in seq_len(nrow(benchmark_spec))) {
  scenario <- benchmark_spec$scenario_id[spec_index]
  scenario_rows <- summary[
    summary$scenario_id == scenario & summary$moment_type == "absolute",
    , drop = FALSE
  ]
  for (sample_size in sort(unique(scenario_rows$n_x))) {
    counter <- counter + 1L
    rows <- scenario_rows[scenario_rows$n_x == sample_size, , drop = FALSE]
    positions <- match(c("K1", "K2", "K3"), rows$interval_id)
    stopifnot(nrow(rows) == 3L, !anyNA(positions))
    coverage_field <- if (scenario == "MIXED_SIGN") {
      "joint_continuum_band_coverage"
    } else {
      "continuum_band_coverage"
    }
    false_field <- if (scenario == "MIXED_SIGN") {
      "joint_any_primary_false_report_rate"
    } else {
      "any_primary_false_report_rate"
    }
    table_02_rows[[counter]] <- data.frame(
      design = benchmark_spec$design[spec_index],
      n_x_equals_n_y = as.integer(sample_size),
      K1_band_coverage = format_three(rows[[coverage_field]][positions[1L]]),
      K2_band_coverage = format_three(rows[[coverage_field]][positions[2L]]),
      K3_band_coverage = format_three(rows[[coverage_field]][positions[3L]]),
      K1_false_report_rate = format_three(rows[[false_field]][positions[1L]]),
      K2_false_report_rate = format_three(rows[[false_field]][positions[2L]]),
      K3_false_report_rate = format_three(rows[[false_field]][positions[3L]]),
      stringsAsFactors = FALSE
    )
  }
}
table_02 <- do.call(rbind, table_02_rows)
row.names(table_02) <- NULL
stopifnot(nrow(table_02) == 10L)

root_collection <- read_csv(
  summary_directory, "root_collection_summary_deduplicated.csv"
)
source_table_03 <- read_csv(article_directory, "table_03_root_inference.csv")
scenario_order <- c("LN_MODERATE", "LN_WEAK", "RARE_CROSSING_P05")
selected <- root_collection[root_collection$scenario_id %in% scenario_order, , drop = FALSE]
selected <- selected[order(match(selected$scenario_id, scenario_order), selected$n_x), , drop = FALSE]
multiplier <- source_table_03[
  source_table_03$interval_method == "Multiplier bootstrap" &
    source_table_03$cell_id %in% selected$cell_id,
  , drop = FALSE
]
positions <- match(selected$cell_id, multiplier$cell_id)
stopifnot(
  nrow(selected) == 12L, !anyNA(positions),
  all(selected$root_collection_joint_isolation_reports ==
        selected$root_multiplier_report_count),
  all(selected$root_collection_error_denominator ==
        selected$root_multiplier_report_count)
)
design_labels <- c(
  LN_MODERATE = "Moderate lognormal",
  LN_WEAK = "Weak lognormal",
  RARE_CROSSING_P05 = "Rare-tail, pi=.05"
)
table_03 <- data.frame(
  design = unname(design_labels[selected$scenario_id]),
  n_x_equals_n_y = as.integer(selected$n_x),
  reports = as.integer(selected$root_multiplier_report_count),
  report_rate = format_three(selected$root_multiplier_report_rate),
  coverage_conditional_on_report =
    format_three(selected$root_multiplier_conditional_coverage),
  report_and_cover_rate =
    format_three(selected$root_multiplier_report_and_cover_rate),
  bias_conditional_on_isolation =
    format_three(selected$root_collection_bias_conditional_on_joint_isolation),
  rmse_conditional_on_isolation =
    format_three(selected$root_collection_rmse_conditional_on_joint_isolation),
  length_conditional_on_report = format_three(
    multiplier$mean_total_length_conditional_on_report[positions]
  ),
  stringsAsFactors = FALSE
)

source_table_04 <- read_csv(article_directory, "table_04_special_power.csv")
source_table_04 <- source_table_04[order(
  match(source_table_04$scenario_id,
        c("TANGENCY_SV_POWER", "TWO_ROOT_SV_POWER")),
  source_table_04$n_x
), , drop = FALSE]
stopifnot(
  nrow(source_table_04) == 6L,
  identical(
    as.character(source_table_04$scenario_id),
    rep(c("TANGENCY_SV_POWER", "TWO_ROOT_SV_POWER"), each = 3L)
  )
)
tangency <- source_table_04[seq_len(3L), , drop = FALSE]
two_root <- source_table_04[4:6, , drop = FALSE]
empty <- rep(NA_character_, 3L)
table_04 <- rbind(
  data.frame(
    panel = "A", design = "Selected tangency",
    n_x_equals_n_y = as.integer(tangency$n_x),
    band_coverage = format_three(tangency$continuum_band_coverage),
    root_in_outer_set = format_three(tangency$outer_set_truth_coverage),
    length_conditional_on_band = format_three(
      tangency$mean_tangency_global_outer_total_length_conditional_on_continuum_band_coverage
    ),
    hausdorff_conditional_on_band = format_three(
      tangency$mean_tangency_global_outer_hausdorff_conditional_on_continuum_band_coverage
    ),
    scaled_length = format_three(
      tangency$mean_tangency_global_scaled_outer_total_length_conditional_on_continuum_band_coverage
    ),
    scaled_hausdorff = format_three(
      tangency$mean_tangency_global_scaled_outer_hausdorff_conditional_on_continuum_band_coverage
    ),
    outer_set_cover = empty, joint_isolation = empty, reversal = empty,
    exact_count = empty, coverage_conditional_on_report = empty,
    report_and_cover_rate = empty, stringsAsFactors = FALSE
  ),
  data.frame(
    panel = "B", design = "Selected two-root design",
    n_x_equals_n_y = as.integer(two_root$n_x),
    band_coverage = format_three(two_root$continuum_band_coverage),
    root_in_outer_set = empty, length_conditional_on_band = empty,
    hausdorff_conditional_on_band = empty, scaled_length = empty,
    scaled_hausdorff = empty,
    outer_set_cover = format_three(two_root$outer_set_truth_coverage),
    joint_isolation = format_three(
      two_root$root_collection_joint_isolation_report_rate
    ),
    reversal = format_three(two_root$reversal_certification_rate),
    exact_count = format_three(two_root$exact_count_certification_rate),
    coverage_conditional_on_report = format_three(
      two_root$root_multiplier_conditional_coverage
    ),
    report_and_cover_rate = format_three(
      two_root$root_multiplier_report_and_cover_rate
    ),
    stringsAsFactors = FALSE
  )
)
row.names(table_04) <- NULL

tail <- read_csv(summary_directory, "tail_stratified.csv")
tail <- tail[
  tail$cell_id == "M019_RARE_CROSSING_P05_n100" & tail$interval_id == "K3",
  , drop = FALSE
]
tail_order <- match(tail$tail_count_stratum, c("0", "1", "2", "3-5", "6+"))
stopifnot(nrow(tail) == 5L, !anyNA(tail_order), !anyDuplicated(tail_order))
tail <- tail[order(tail_order), , drop = FALSE]
table_05 <- data.frame(
  observed_rare_draws = c("0", "1", "2", "3--5", "6 or more"),
  replications = as.integer(tail$replications),
  band_coverage = format_three(tail$continuum_band_coverage),
  false_report_rate = format_three(tail$any_primary_false_report_rate),
  stringsAsFactors = FALSE
)

ablation_metric <- function(cell_id, metric) {
  rows <- binary[
    binary$cell_id == cell_id &
      binary$moment_type == "absolute" &
      binary$interval_id == "K3" &
      binary$metric == metric,
    , drop = FALSE
  ]
  stopifnot(
    nrow(rows) == 1L,
    as.integer(rows$denominator) == 2032L,
    abs(
      as.numeric(rows$estimate) -
        as.integer(rows$successes) / as.integer(rows$denominator)
    ) <= 1e-14
  )
  rows
}

ablation_spec <- data.frame(
  panel = c("A", "A", "A", "B", "B"),
  cell_id = c(
    "M004_LN_EQUAL_n500",
    "M008_LN_NO_ROOT_n500",
    "M052_TWO_ROOT_SV_POWER_n2000",
    "M008_LN_NO_ROOT_n500",
    "M052_TWO_ROOT_SV_POWER_n2000"
  ),
  design = c(
    "LN-EQUAL", "LN-NO-ROOT", "SELECTED-TWO-ROOT",
    "LN-NO-ROOT", "SELECTED-TWO-ROOT"
  ),
  n_x_equals_n_y = c(500L, 500L, 2000L, 500L, 2000L),
  quantity = c(
    "Curve coverage", "Curve coverage", "Curve coverage",
    "No root", "Reversal"
  ),
  continuum_metric = c(
    rep("continuum_band_covers_truth", 3L),
    "continuum_no_root_certified", "reversal_certified"
  ),
  simultaneous_grid_metric = c(
    rep("band_covers_truth_on_grid", 3L),
    "grid_only_no_root", "grid_only_reversal_certified"
  ),
  pointwise_grid_metric = c(
    rep("pointwise_grid_curve_coverage", 3L),
    "pointwise_no_root_report", "pointwise_reversal_report"
  ),
  stringsAsFactors = FALSE
)

ablation_rows <- lapply(seq_len(nrow(ablation_spec)), function(index) {
  spec <- ablation_spec[index, , drop = FALSE]
  continuum <- ablation_metric(spec$cell_id, spec$continuum_metric)
  simultaneous_grid <- ablation_metric(
    spec$cell_id, spec$simultaneous_grid_metric
  )
  pointwise_grid <- ablation_metric(spec$cell_id, spec$pointwise_grid_metric)
  stopifnot(
    identical(continuum$denominator, simultaneous_grid$denominator),
    identical(continuum$denominator, pointwise_grid$denominator)
  )
  data.frame(
    panel = spec$panel,
    cell_id = spec$cell_id,
    design = spec$design,
    n_x_equals_n_y = spec$n_x_equals_n_y,
    quantity = spec$quantity,
    denominator = as.integer(continuum$denominator),
    continuum_successes = as.integer(continuum$successes),
    continuum_simultaneous = format_three(continuum$estimate),
    grid_simultaneous_successes = as.integer(simultaneous_grid$successes),
    grid_simultaneous = format_three(simultaneous_grid$estimate),
    grid_pointwise_successes = as.integer(pointwise_grid$successes),
    grid_pointwise = format_three(pointwise_grid$estimate),
    stringsAsFactors = FALSE
  )
})
table_07 <- do.call(rbind, ablation_rows)
row.names(table_07) <- NULL

write_table(table_01, "table_01_design.csv")
write_table(table_02, "table_02_benchmark.csv")
write_table(table_03, "table_03_root_inference.csv")
write_table(table_04, "table_04_special_power.csv")
write_table(table_05, "table_05_tail_strata.csv")
write_table(table_07, "table_07_ablation.csv")
cat("FINAL PUBLIC SIMULATION TABLE LAYER: PASSED\n")
