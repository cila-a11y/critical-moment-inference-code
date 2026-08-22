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
assert_true(identical(paste(R.version$major, R.version$minor, sep = "."),
                      config$compute$r_version),
            sprintf("Canonical execution requires R %s.", config$compute$r_version))
assert_true(isTRUE(capabilities("cairo")),
            "Canonical execution requires Cairo support for headless PNG rendering.")
verify_design_lock()
verify_source_manifest()
source_hillstrom_engine(config)
source_manifest_sha256 <- sha256_file(file.path(
  hill_root, "SOURCE_MANIFEST_SHA256.txt"
))
design_lock_sha256 <- sha256_file(file.path(
  hill_root, "DESIGN_LOCK_SHA256.txt"
))

analysis <- config$analysis
assert_true(analysis$alpha > 0 && analysis$alpha < 1,
            "alpha must lie in (0,1).")
assert_true(analysis$bootstrap_reps >= 2L,
            "At least two multiplier replications are required.")
assert_true(analysis$bootstrap_batch_size >= 1L,
            "bootstrap_batch_size must be positive.")

result_parent <- file.path(hill_root, "results", config$application_id)
final_directory <- file.path(result_parent, run_id)
dir.create(result_parent, recursive = TRUE, showWarnings = FALSE)
assert_true(!file.exists(final_directory),
            sprintf("Refusing to overwrite existing run directory: %s", final_directory))
temporary_directory <- tempfile(pattern = paste0(".", run_id, "."),
                                tmpdir = result_parent)
dir.create(temporary_directory, recursive = FALSE, showWarnings = FALSE,
           mode = "0750")

started <- Sys.time()
log_message("Validating locked Hillstrom data contract.")
validated <- read_and_validate_hillstrom(data_path, config, strict = TRUE)
x <- validated$x
y <- validated$y

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
assert_true(length(bootstrap_orders) == analysis$expected_bootstrap_grid_size,
            "Unexpected bootstrap-grid size.")
assert_true(length(enclosure_orders) == analysis$expected_enclosure_grid_size,
            "Unexpected enclosure-grid size.")
assert_true(length(audit_orders) == analysis$expected_audit_grid_size,
            "Unexpected audit-grid size.")

do.call(RNGkind, as.list(analysis$rng_kind))
set.seed(analysis$bootstrap_seed)

seed_hash <- function(seed) {
  temporary <- tempfile(tmpdir = temporary_directory, fileext = ".rds")
  on.exit(unlink(temporary), add = TRUE)
  saveRDS(seed, temporary, version = 3L, compress = FALSE)
  sha256_file(temporary)
}
initial_seed_hash <- seed_hash(.Random.seed)

log_message("Computing the locked B=", analysis$bootstrap_reps,
            " simultaneous multiplier band.")
nested <- nested_multiplier_bands(
  x = x,
  y = y,
  orders = bootstrap_orders,
  moment_types = analysis$moment_types,
  order_intervals = analysis$order_intervals,
  bootstrap_reps = analysis$bootstrap_reps,
  alpha = analysis$alpha,
  multiplier_distribution = analysis$multiplier_distribution,
  bootstrap_batch_size = analysis$bootstrap_batch_size,
  variance_tolerance = analysis$variance_tolerance,
  guard_max_levels = analysis$guard_max_levels,
  guard_max_nodes = analysis$enclosure_max_nodes
)
final_seed_hash <- seed_hash(.Random.seed)
fit <- nested$fits$K1
band <- extract_band(fit, analysis$moment_types)

assert_true(!fit$guard_failed && isTRUE(fit$guard$valid),
            "The continuum relative-variance guard failed.")
assert_true(is.finite(fit$critical_value) && fit$critical_value > 0,
            "The simultaneous critical value is not positive and finite.")
assert_true(all(is.finite(unlist(band))),
            "The construction-grid band contains a nonfinite value.")

moment_x <- colMeans(fit$feature_x[[analysis$moment_types]])
moment_y <- colMeans(fit$feature_y[[analysis$moment_types]])
standard_error <- sqrt(pmax(band$variance_hat, 0) / fit$n_eff)
studentized <- band$delta_hat / standard_error
moment_band <- data.frame(
  order = band$orders,
  moment_x = moment_x,
  moment_y = moment_y,
  delta_hat = band$delta_hat,
  variance_hat = band$variance_hat,
  standard_error = standard_error,
  critical_value = fit$critical_value,
  half_width = band$half_width,
  lower = band$lower,
  upper = band$upper,
  studentized = studentized,
  sign_status = ifelse(
    band$lower > analysis$numerical_tolerance, "positive",
    ifelse(band$upper < -analysis$numerical_tolerance, "negative", "unresolved")
  ),
  stringsAsFactors = FALSE
)

log_message("Evaluating the separate fine audit grid and continuum enclosure.")
audit <- evaluate_fitted_band(
  fit, audit_orders, analysis$moment_types,
  chunk_size = analysis$evaluation_chunk_size
)
audit_standard_error <- sqrt(pmax(audit$variance_hat, 0) / fit$n_eff)
moment_band_audit <- data.frame(
  order = audit$orders,
  delta_hat = audit$delta_hat,
  variance_hat = audit$variance_hat,
  standard_error = audit_standard_error,
  critical_value = fit$critical_value,
  half_width = audit$half_width,
  lower = audit$lower,
  upper = audit$upper,
  studentized = audit$delta_hat / audit_standard_error,
  sign_status = ifelse(
    audit$lower > analysis$numerical_tolerance, "positive",
    ifelse(audit$upper < -analysis$numerical_tolerance, "negative", "unresolved")
  ),
  stringsAsFactors = FALSE
)
audit_numeric <- moment_band_audit[, setdiff(
  names(moment_band_audit), "sign_status"
), drop = FALSE]
assert_true(all(vapply(audit_numeric, function(value) all(is.finite(value)),
                       logical(1L))),
            "The audit-grid band contains a nonfinite value.")

enclosure <- adaptive_continuum_enclosure(
  fit = fit,
  moment_type = analysis$moment_types,
  initial_grid = enclosure_orders,
  numerical_tolerance = analysis$numerical_tolerance,
  safety_margin = analysis$enclosure_safety_margin,
  max_levels = analysis$enclosure_max_levels,
  max_nodes = analysis$enclosure_max_nodes,
  refine_unresolved_levels = analysis$enclosure_refine_unresolved_levels,
  maximum_cell_half_width = analysis$maximum_enclosure_half_width,
  roundoff_inflation = analysis$roundoff_inflation,
  state_cache = nested$state_cache,
  constant_cache = nested$constant_cache
)
cells <- enclosure$cells

partition_tolerance <- 1e-12
assert_true(nrow(cells) >= 1L && all(cells$right > cells$left),
            "The continuum enclosure contains a malformed cell.")
assert_true(abs(cells$left[1L] - interval[1L]) <= partition_tolerance &&
              abs(tail(cells$right, 1L) - interval[2L]) <= partition_tolerance,
            "The continuum enclosure does not cover the locked interval.")
if (nrow(cells) > 1L) {
  assert_true(all(abs(head(cells$right, -1L) - tail(cells$left, -1L)) <=
                    partition_tolerance),
              "The continuum enclosure has a gap or overlap.")
}
assert_true(all(rowSums(cbind(cells$positive, cells$negative, cells$outer)) == 1L),
            "A continuum cell has an invalid classification.")
assert_true(all(cells$lower_bound <= cells$upper_bound),
            "A continuum cell has reversed deterministic bounds.")
assert_true(!any(cells$variance_unresolved_at_limit),
            "The variance enclosure remained unresolved at a numerical limit.")

audit_containment <- vapply(seq_len(nrow(moment_band_audit)), function(index) {
  order <- moment_band_audit$order[index]
  candidates <- which(cells$left <= order + partition_tolerance &
                        cells$right >= order - partition_tolerance)
  length(candidates) >= 1L &&
    moment_band_audit$lower[index] >=
      max(cells$lower_bound[candidates]) - partition_tolerance &&
    moment_band_audit$upper[index] <=
      min(cells$upper_bound[candidates]) + partition_tolerance
}, logical(1L))
audit_sign_consistent <- vapply(seq_len(nrow(moment_band_audit)), function(index) {
  order <- moment_band_audit$order[index]
  candidates <- which(cells$left <= order + partition_tolerance &
                        cells$right >= order - partition_tolerance)
  length(candidates) >= 1L && all(vapply(candidates, function(cell_index) {
    (!cells$positive[cell_index] ||
       moment_band_audit$lower[index] > analysis$numerical_tolerance) &&
      (!cells$negative[cell_index] ||
         moment_band_audit$upper[index] < -analysis$numerical_tolerance)
  }, logical(1L)))
}, logical(1L))
audit_checks <- data.frame(
  check = c("audit_nodes", "whole_cell_containment_violations",
            "cell_sign_contradictions"),
  value = c(nrow(moment_band_audit), sum(!audit_containment),
            sum(!audit_sign_consistent)),
  status = c(
    if (nrow(moment_band_audit) == analysis$expected_audit_grid_size) "PASSED" else "FAILED",
    if (all(audit_containment)) "PASSED" else "FAILED",
    if (all(audit_sign_consistent)) "PASSED" else "FAILED"
  ),
  stringsAsFactors = FALSE
)
assert_true(all(audit_checks$status == "PASSED"),
            "The separate fine audit grid contradicted the continuum enclosure.")

outer <- outer_set_summary(cells, true_roots = numeric(0L))
outer_components <- outer$components
outer_components$component <- seq_len(nrow(outer_components))
outer_components$width <- outer_components$right - outer_components$left
outer_components <- outer_components[, c("component", "left", "right", "width"),
                                     drop = FALSE]

positive_intervals <- merge_adjacent_cells(cells, cells$positive, "positive")
negative_intervals <- merge_adjacent_cells(cells, cells$negative, "negative")
certified_sign_intervals <- rbind(positive_intervals, negative_intervals)
if (nrow(certified_sign_intervals)) {
  certified_sign_intervals <- certified_sign_intervals[
    order(certified_sign_intervals$left, certified_sign_intervals$right), , drop = FALSE
  ]
  certified_sign_intervals$width <- certified_sign_intervals$right -
    certified_sign_intervals$left
} else {
  certified_sign_intervals$width <- numeric(0L)
}

certified_brackets <- brackets_from_sign_intervals(certified_sign_intervals)
assert_bracket_invariants(
  certified_brackets, certified_sign_intervals, interval,
  tolerance = partition_tolerance
)

describe_group <- function(values, segment, role) {
  positive <- values[values > 0]
  data.frame(
    role = role,
    segment = segment,
    n = length(values),
    n_positive = length(positive),
    positive_rate = mean(values > 0),
    mean_spend = mean(values),
    sd_spend = stats::sd(values),
    median_spend = stats::median(values),
    conditional_mean_positive = if (length(positive)) mean(positive) else NA_real_,
    maximum_spend = max(values),
    stringsAsFactors = FALSE
  )
}
descriptive <- rbind(
  describe_group(x, analysis$x_segment, "X"),
  describe_group(y, analysis$y_segment, "Y")
)

tail_diagnostics_table <- do.call(rbind, lapply(c("X", "Y"), function(role) {
  values <- if (identical(role, "X")) x else y
  segment <- if (identical(role, "X")) analysis$x_segment else analysis$y_segment
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
      }
    ))
  }))
}))

anchor_moment_indices <- vapply(analysis$anchor_orders, function(order) {
  index <- which.min(abs(moment_band$order - order))
  assert_true(abs(moment_band$order[index] - order) <= 1e-12,
              "An anchor order is absent from the construction grid.")
  index
}, integer(1L))
anchor_moments <- do.call(rbind, lapply(seq_along(analysis$anchor_orders), function(index) {
  order <- analysis$anchor_orders[index]
  candidates <- which(cells$left <= order + partition_tolerance &
                        cells$right >= order - partition_tolerance)
  assert_true(length(candidates) >= 1L,
              "An anchor order is absent from the continuum partition.")
  endpoint_lower <- max(cells$lower_bound[candidates])
  endpoint_upper <- min(cells$upper_bound[candidates])
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
    moment_x = moment_band$moment_x[anchor_moment_indices[index]],
    moment_y = moment_band$moment_y[anchor_moment_indices[index]],
    delta_hat = moment_band$delta_hat[anchor_moment_indices[index]],
    implemented_lower = moment_band$lower[anchor_moment_indices[index]],
    implemented_upper = moment_band$upper[anchor_moment_indices[index]],
    endpoint_lower = endpoint_lower,
    endpoint_upper = endpoint_upper,
    anchor_sign_status = status,
    stringsAsFactors = FALSE
  )
}))

uniform_y <- nrow(cells) > 0L && all(cells$positive)
uniform_x <- nrow(cells) > 0L && all(cells$negative)
empty_outer <- nrow(cells) > 0L && !any(cells$outer)
assert_true(
  !empty_outer || xor(uniform_y, uniform_x),
  paste(
    "An empty outer set was obtained without a uniform certified sign;",
    "the continuum classification is internally inconsistent."
  )
)
no_root <- empty_outer && (uniform_y || uniform_x)
conclusion <- if (uniform_y) {
  "UNIFORM_Y_DOMINANCE_NO_ROOT"
} else if (uniform_x) {
  "UNIFORM_X_DOMINANCE_NO_ROOT"
} else if (nrow(certified_brackets) > 0L) {
  "REVERSAL_AT_LEAST_ONE_ROOT"
} else {
  "INCONCLUSIVE"
}

application_summary <- data.frame(
  application_id = config$application_id,
  x_segment = analysis$x_segment,
  y_segment = analysis$y_segment,
  n_x = length(x), n_y = length(y), n_eff = fit$n_eff,
  p_min = interval[1L], p_max = interval[2L],
  alpha = analysis$alpha,
  bootstrap_reps = analysis$bootstrap_reps,
  bootstrap_seed = analysis$bootstrap_seed,
  multiplier_distribution = analysis$multiplier_distribution,
  bootstrap_grid_size = length(bootstrap_orders),
  enclosure_initial_grid_size = length(enclosure_orders),
  audit_grid_size = length(audit_orders),
  critical_value = fit$critical_value,
  guard_valid = fit$guard$valid,
  guard_minimum_certified_ratio = fit$guard$minimum_certified_ratio,
  guard_nodes_used = fit$guard$nodes_used,
  guard_maximum_level = fit$guard$maximum_level,
  enclosure_nodes_used = enclosure$nodes_used,
  enclosure_maximum_level = enclosure$maximum_level,
  enclosure_limit_hit = enclosure$limit_hit,
  enclosure_variance_limit_hit = enclosure$variance_limit_hit,
  enclosure_statistical_limit_hit = enclosure$statistical_limit_hit,
  outer_component_count = outer$component_count,
  outer_total_length = outer$total_length,
  certified_alternations = nrow(certified_brackets),
  certified_bracket_count = nrow(certified_brackets),
  no_root_on_K_certified = no_root,
  uniform_y_dominance_certified = uniform_y,
  uniform_x_dominance_certified = uniform_x,
  conclusion = conclusion,
  structural_crossing_budget_used = FALSE,
  stringsAsFactors = FALSE
)

bootstrap_suprema <- data.frame(
  replication = seq_len(analysis$bootstrap_reps),
  supremum = as.numeric(nested$bootstrap_suprema[, "K1"])
)
assert_true(all(is.finite(bootstrap_suprema$supremum)),
            "A multiplier-bootstrap supremum is nonfinite.")
quantile_index <- min(
  analysis$bootstrap_reps,
  ceiling((analysis$bootstrap_reps + 1L) * (1 - analysis$alpha))
)
recomputed_critical_value <- sort.int(
  bootstrap_suprema$supremum, partial = quantile_index
)[quantile_index]
critical_tolerance <- 64 * .Machine$double.eps *
  max(1, abs(fit$critical_value), abs(recomputed_critical_value))
assert_true(
  abs(fit$critical_value - recomputed_critical_value) <= critical_tolerance,
  "The stored critical value is not the locked empirical bootstrap quantile."
)

rng_state_hashes <- data.frame(
  state = c("initial", "final"),
  sha256 = c(initial_seed_hash, final_seed_hash),
  stringsAsFactors = FALSE
)

provenance_path <- file.path(hill_root, "status", config$application_id,
                             "data_provenance.csv")
assert_true(file.exists(provenance_path),
            "Data provenance is missing; run bin/acquire_data.sh first.")
data_provenance <- read.csv(provenance_path, stringsAsFactors = FALSE,
                            check.names = FALSE)
expected_provenance_columns <- c(
  "dataset", "source_mode", "source_used", "source_page", "accessed_utc",
  "raw_filename", "bytes", "sha256", "license_status", "redistributed"
)
assert_true(
  nrow(data_provenance) == 1L &&
    identical(names(data_provenance), expected_provenance_columns),
  "data_provenance.csv must contain exactly one locked provenance row."
)
assert_true(
  identical(as.character(data_provenance$sha256), config$data$sha256) &&
    as.numeric(data_provenance$bytes) == as.numeric(config$data$bytes) &&
    identical(as.character(data_provenance$raw_filename),
              config$data$raw_filename) &&
    identical(toupper(as.character(data_provenance$redistributed)), "FALSE"),
  "The data-provenance row disagrees with the locked data contract."
)

table_csv <- rbind(
  data.frame(
    panel = "A", role = descriptive$role, segment = descriptive$segment,
    order = NA_real_, n = descriptive$n,
    n_positive = descriptive$n_positive,
    positive_rate = descriptive$positive_rate,
    mean_spend = descriptive$mean_spend,
    sd_spend = descriptive$sd_spend,
    conditional_mean_positive = descriptive$conditional_mean_positive,
    moment_x = NA_real_, moment_y = NA_real_, delta_hat = NA_real_,
    implemented_lower = NA_real_, implemented_upper = NA_real_,
    endpoint_lower = NA_real_, endpoint_upper = NA_real_,
    anchor_sign_status = NA_character_, stringsAsFactors = FALSE
  ),
  data.frame(
    panel = "B", role = NA_character_, segment = NA_character_,
    order = anchor_moments$order, n = NA_integer_, n_positive = NA_integer_,
    positive_rate = NA_real_, mean_spend = NA_real_,
    sd_spend = NA_real_,
    conditional_mean_positive = NA_real_,
    moment_x = anchor_moments$moment_x,
    moment_y = anchor_moments$moment_y,
    delta_hat = anchor_moments$delta_hat,
    implemented_lower = anchor_moments$implemented_lower,
    implemented_upper = anchor_moments$implemented_upper,
    endpoint_lower = anchor_moments$endpoint_lower,
    endpoint_upper = anchor_moments$endpoint_upper,
    anchor_sign_status = anchor_moments$anchor_sign_status,
    stringsAsFactors = FALSE
  )
)

atomic_write_csv(validated$validation,
                 file.path(temporary_directory, "data_validation.csv"))
atomic_write_csv(data_provenance,
                 file.path(temporary_directory, "data_provenance.csv"))
atomic_write_csv(descriptive,
                 file.path(temporary_directory, "descriptive_statistics.csv"))
atomic_write_csv(moment_band,
                 file.path(temporary_directory, "moment_band.csv"))
atomic_write_csv(moment_band_audit,
                 file.path(temporary_directory, "moment_band_audit.csv"))
atomic_write_csv(anchor_moments,
                 file.path(temporary_directory, "anchor_moments.csv"))
atomic_write_csv(table_csv,
                 file.path(temporary_directory,
                           "table_01_hillstrom_application.csv"))
atomic_write_csv(cells,
                 file.path(temporary_directory, "enclosure_cells.csv"))
atomic_write_csv(outer_components,
                 file.path(temporary_directory, "outer_set_components.csv"))
atomic_write_csv(certified_sign_intervals,
                 file.path(temporary_directory, "certified_sign_intervals.csv"))
atomic_write_csv(certified_brackets,
                 file.path(temporary_directory, "certified_brackets.csv"))
atomic_write_csv(tail_diagnostics_table,
                 file.path(temporary_directory, "tail_diagnostics.csv"))
atomic_write_csv(audit_checks,
                 file.path(temporary_directory, "audit_checks.csv"))
atomic_write_csv(application_summary,
                 file.path(temporary_directory, "application_summary.csv"))
atomic_write_csv(bootstrap_suprema,
                 file.path(temporary_directory, "bootstrap_suprema.csv"))
atomic_write_csv(rng_state_hashes,
                 file.path(temporary_directory, "rng_state_hashes.csv"))

table_lines <- render_hillstrom_table(
  descriptive, anchor_moments, application_summary, outer_components,
  certified_brackets
)
atomic_write_lines(table_lines,
                   file.path(temporary_directory, "table_01_hillstrom_application.tex"))

draw_figure <- function(device) {
  device()
  old <- graphics::par(no.readonly = TRUE)
  on.exit({
    try(graphics::par(old), silent = TRUE)
    grDevices::dev.off()
  }, add = TRUE)
  graphics::par(mfrow = c(1, 2), mar = c(4.2, 4.5, 2.1, 0.9),
                mgp = c(2.5, 0.7, 0), las = 1, lend = "butt")
  y_range <- range(c(moment_band_audit$lower, moment_band_audit$upper, 0))
  graphics::plot(moment_band_audit$order, moment_band_audit$delta_hat,
                 type = "n", xlab = "Moment order p",
                 ylab = expression(hat(Delta)(p)), ylim = y_range)
  cell_colours <- c(
    positive = "#E5F5E0", negative = "#FDE0DD", outer = "#EEEEEE"
  )
  graphics::rect(
    min(moment_band_audit$order), y_range[1L],
    max(moment_band_audit$order), y_range[2L],
    col = cell_colours[["outer"]], border = NA
  )
  if (nrow(certified_sign_intervals)) {
    for (index in seq_len(nrow(certified_sign_intervals))) {
      interval_colour <- cell_colours[[certified_sign_intervals$sign[index]]]
      graphics::rect(
        certified_sign_intervals$left[index], y_range[1L],
        certified_sign_intervals$right[index], y_range[2L],
        col = interval_colour, border = NA
      )
    }
  }
  band_colour <- grDevices::adjustcolor("#377EB8", alpha.f = 0.30)
  graphics::polygon(
    c(moment_band_audit$order, rev(moment_band_audit$order)),
    c(moment_band_audit$lower, rev(moment_band_audit$upper)),
    col = band_colour, border = NA
  )
  graphics::lines(moment_band_audit$order, moment_band_audit$delta_hat,
                  col = "#084594", lwd = 2)
  if (nrow(certified_brackets)) {
    graphics::segments(
      certified_brackets$lower, 0, certified_brackets$upper, 0,
      col = "#984EA3", lwd = 5, lend = 1
    )
  }
  graphics::abline(h = 0, col = "#333333", lty = 2)
  graphics::abline(v = 1, col = "#777777", lty = 3)
  legend_labels <- c("Estimated contrast", "Nominal 95% simultaneous band")
  legend_colours <- c("#084594", "#8DBBDD")
  legend_lines <- c(1, NA)
  legend_widths <- c(2, NA)
  legend_points <- c(NA, 15)
  add_region_legend <- function(label, colour) {
    legend_labels <<- c(legend_labels, label)
    legend_colours <<- c(legend_colours, colour)
    legend_lines <<- c(legend_lines, NA)
    legend_widths <<- c(legend_widths, NA)
    legend_points <<- c(legend_points, 15)
  }
  if (any(cells$positive)) {
    add_region_legend("Certified positive region", cell_colours[["positive"]])
  }
  if (any(cells$negative)) {
    add_region_legend("Certified negative region", cell_colours[["negative"]])
  }
  if (any(cells$outer)) {
    add_region_legend("Unresolved outer region", cell_colours[["outer"]])
  }
  legend_labels <- c(legend_labels, "p = 1 reference")
  legend_colours <- c(legend_colours, "#777777")
  legend_lines <- c(legend_lines, 3)
  legend_widths <- c(legend_widths, 1)
  legend_points <- c(legend_points, NA)
  if (nrow(certified_brackets)) {
    legend_labels <- c(legend_labels, "certified root bracket")
    legend_colours <- c(legend_colours, "#984EA3")
    legend_lines <- c(legend_lines, 1)
    legend_widths <- c(legend_widths, 5)
    legend_points <- c(legend_points, NA)
  }
  graphics::legend(
    "topleft", legend = legend_labels, col = legend_colours,
    lty = legend_lines, lwd = legend_widths,
    pch = legend_points, pt.cex = 1.1, bty = "n", cex = 0.61,
    y.intersp = 0.9
  )
  graphics::mtext("(a) Moment contrast", side = 3, line = 0.35,
                  adj = 0, font = 2, cex = 0.90)

  z_range <- range(c(moment_band_audit$studentized,
                     -fit$critical_value, fit$critical_value, 0))
  graphics::plot(moment_band_audit$order, moment_band_audit$studentized,
                 type = "l", lwd = 2, col = "#084594",
                 xlab = "Moment order p", ylab = "Studentized contrast",
                 ylim = z_range)
  graphics::abline(h = c(-fit$critical_value, fit$critical_value),
                   col = "#B2182B", lty = 2)
  graphics::abline(h = 0, col = "#333333")
  graphics::abline(v = 1, col = "#777777", lty = 3)
  graphics::legend(
    "bottomleft",
    legend = c(
      "Studentized contrast",
      sprintf("Simultaneous thresholds: +/- %.3f", fit$critical_value),
      "p = 1 reference"
    ),
    col = c("#084594", "#B2182B", "#777777"),
    lty = c(1, 2, 3), lwd = c(2, 1, 1),
    bty = "n", cex = 0.66
  )
  graphics::mtext("(b) Studentized contrast", side = 3, line = 0.35,
                  adj = 0, font = 2, cex = 0.90)
}
draw_figure(function() grDevices::cairo_pdf(
  file.path(temporary_directory, "figure_01_hillstrom_moment_contrast.pdf"),
  width = 7.2, height = 3.5, family = "sans", bg = "white"
))
draw_figure(function() grDevices::png(
  file.path(temporary_directory, "figure_01_hillstrom_moment_contrast.png"),
  width = 2160, height = 1050, res = 300, type = "cairo-png", bg = "white"
))

run_metadata <- data.frame(
  field = c(
    "application_id", "application_version", "run_id", "started_utc",
    "finished_utc", "host", "slurm_job_id", "r_version", "platform",
    "engine_version", "engine_manifest_sha256", "config_sha256",
    "source_manifest_sha256", "design_lock_sha256",
    "data_sha256", "data_bytes", "rng_kind", "bootstrap_seed",
    "initial_rng_state_sha256", "final_rng_state_sha256"
  ),
  value = c(
    config$application_id, config$application_version, run_id,
    format(started, tz = "UTC", format = "%Y-%m-%dT%H:%M:%SZ"),
    timestamp_utc(), Sys.info()[["nodename"]],
    Sys.getenv("SLURM_JOB_ID", unset = NA_character_), R.version.string,
    R.version$platform, config$engine$version, config$engine$manifest_sha256,
    sha256_file(config$config_path), source_manifest_sha256, design_lock_sha256,
    config$data$sha256, config$data$bytes,
    paste(analysis$rng_kind, collapse = ";"), analysis$bootstrap_seed,
    initial_seed_hash, final_seed_hash
  ),
  stringsAsFactors = FALSE
)
atomic_write_csv(run_metadata,
                 file.path(temporary_directory, "run_metadata.csv"))
graphics_session <- c(
  utils::capture.output(sessionInfo()),
  "",
  sprintf("bitmap_type=%s", getOption("bitmapType", default = NA_character_)),
  sprintf("capabilities_cairo=%s", capabilities("cairo")),
  utils::capture.output(grDevices::grSoftVersion())
)
atomic_write_lines(graphics_session,
                   file.path(temporary_directory, "session_info.txt"))

compact_result <- list(
  config = config,
  application_summary = application_summary,
  descriptive_statistics = descriptive,
  moment_band = moment_band,
  moment_band_audit = moment_band_audit,
  enclosure_cells = cells,
  outer_set_components = outer_components,
  certified_sign_intervals = certified_sign_intervals,
  certified_brackets = certified_brackets,
  tail_diagnostics = tail_diagnostics_table,
  audit_checks = audit_checks,
  rng_state_hashes = rng_state_hashes
)
atomic_save_rds(compact_result,
                file.path(temporary_directory, "application_result_compact.rds"),
                compress = "xz")

write_result_manifest(temporary_directory)
assert_true(verify_result_manifest(temporary_directory),
            "The result manifest failed immediately after creation.")
assert_true(file.rename(temporary_directory, final_directory),
            "Could not atomically publish the completed run directory.")

rm(nested, fit, validated, x, y)
invisible(gc())
log_message("Published immutable run directory: ", final_directory)
cat(sprintf("SCIENTIFIC_CONCLUSION=%s\n", conclusion))
cat("HILLSTROM APPLICATION COMPUTATION: COMPLETED\n")
