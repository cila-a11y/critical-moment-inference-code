# Pure current-release reporting helpers.  Keeping these functions free of filesystem
# side effects lets the node audit exercise the exact aggregation code used by
# summarize_results.R.

split_semicolon_character <- function(value) {
  value <- as.character(value %||% "")
  if (length(value) != 1L || is.na(value) || !nzchar(value)) {
    return(character(0L))
  }
  strsplit(value, ";", fixed = TRUE)[[1L]]
}

split_semicolon_numeric <- function(value) {
  values <- split_semicolon_character(value)
  if (!length(values)) return(numeric(0L))
  parsed <- suppressWarnings(as.numeric(values))
  assert_true(
    all(is.finite(parsed)),
    "A persisted root-collection vector contains a nonfinite value."
  )
  parsed
}

root_collection_details_for_group <- function(group) {
  empty <- list(
    applicable = 0L, joint = 0L, root_count = 0L,
    root_ids = character(0L), true_roots = numeric(0L),
    applicable_mask = rep(FALSE, nrow(group)),
    joint_mask = rep(FALSE, nrow(group)),
    error_matrix = matrix(numeric(0L), nrow = 0L, ncol = 0L),
    errors = numeric(0L)
  )
  required_fields <- c(
    "root_collection_applicable",
    "root_collection_required_root_count",
    "root_collection_joint_isolation_success",
    "root_collection_root_ids", "root_collection_true_roots",
    "root_collection_root_hats", "root_collection_errors"
  )
  if (!nrow(group) || !all(required_fields %in% names(group))) return(empty)

  applicable_flag <- as.logical(group$root_collection_applicable)
  assert_true(
    !anyNA(applicable_flag),
    "Root-collection applicability is missing from a persisted replication."
  )
  applicable <- applicable_flag
  if (!any(applicable)) {
    empty$applicable_mask <- applicable
    return(empty)
  }
  applicable_indices <- which(applicable)
  required_counts <- suppressWarnings(as.integer(
    group$root_collection_required_root_count[applicable]
  ))
  assert_true(
    length(unique(required_counts)) == 1L && required_counts[1L] >= 1L,
    "The required root count changes within a scientific group."
  )
  root_count <- required_counts[1L]
  root_ids <- split_semicolon_character(
    group$root_collection_root_ids[applicable_indices[1L]]
  )
  true_roots <- split_semicolon_numeric(
    group$root_collection_true_roots[applicable_indices[1L]]
  )
  assert_true(
    length(root_ids) == root_count &&
      length(true_roots) == root_count &&
      !anyDuplicated(root_ids),
    "The reference root collection is incomplete or duplicated."
  )
  for (index in applicable_indices) {
    ids_now <- split_semicolon_character(
      group$root_collection_root_ids[index]
    )
    truth_now <- split_semicolon_numeric(
      group$root_collection_true_roots[index]
    )
    assert_true(
      identical(ids_now, root_ids) &&
        length(truth_now) == root_count &&
        max(abs(truth_now - true_roots)) <= 1e-10,
      "Root identifiers or truths change within a scientific group."
    )
  }

  joint_flag <- as.logical(
    group$root_collection_joint_isolation_success
  )
  assert_true(
    !anyNA(joint_flag[applicable]),
    paste(
      "Joint-isolation status is missing from an applicable root",
      "collection."
    )
  )
  joint <- applicable & !is.na(joint_flag) & joint_flag
  for (index in which(applicable & !joint)) {
    assert_true(
      !length(split_semicolon_numeric(
        group$root_collection_root_hats[index]
      )) &&
        !length(split_semicolon_numeric(
          group$root_collection_errors[index]
        )),
      paste(
        "A replication without joint root isolation contains persisted",
        "root estimates or errors."
      )
    )
  }
  error_matrix <- matrix(
    numeric(0L), nrow = 0L, ncol = root_count,
    dimnames = list(NULL, root_ids)
  )
  if (any(joint)) {
    error_matrix <- do.call(rbind, lapply(which(joint), function(index) {
      root_hats <- split_semicolon_numeric(
        group$root_collection_root_hats[index]
      )
      errors <- split_semicolon_numeric(
        group$root_collection_errors[index]
      )
      assert_true(
        length(root_hats) == root_count &&
          length(errors) == root_count &&
          max(abs(root_hats - true_roots - errors)) <= 1e-9,
        "A jointly isolated root collection is incomplete or misaligned."
      )
      errors
    }))
    colnames(error_matrix) <- root_ids
  }
  assert_true(
    nrow(error_matrix) == sum(joint) &&
      length(error_matrix) == sum(joint) * root_count &&
      all(is.finite(error_matrix)),
    "Joint root errors were silently dropped from the conditional denominator."
  )
  list(
    applicable = sum(applicable), joint = sum(joint),
    root_count = root_count, root_ids = root_ids, true_roots = true_roots,
    applicable_mask = applicable, joint_mask = joint,
    error_matrix = error_matrix,
    errors = as.numeric(t(error_matrix))
  )
}

root_collection_errors_for_group <- function(group) {
  details <- root_collection_details_for_group(group)
  list(
    applicable = details$applicable,
    joint = details$joint,
    root_count = details$root_count,
    errors = details$errors
  )
}

across_interval_metric_specification <- c(
  across_interval_joint_band_covers_truth_on_grid =
    "band_covers_truth_on_grid",
  across_interval_joint_band_covers_truth_on_audit_grid =
    "band_covers_truth_on_audit_grid",
  across_interval_joint_continuum_band_covers_truth =
    "continuum_band_covers_truth",
  across_interval_joint_any_primary_false_report =
    "any_primary_false_report"
)

deduplicate_across_interval_events <- function(
    observed, expected_interval_ids,
    metric_specification = across_interval_metric_specification) {
  expected_interval_ids <- sort(unique(as.character(expected_interval_ids)))
  assert_true(
    length(expected_interval_ids) >= 1L,
    "ALL_K aggregation requires at least one configured interval."
  )
  if (!nrow(observed)) {
    template <- data.frame(
      unit_id = integer(0L), rep_id = integer(0L),
      seed_id = character(0L), bootstrap_seed_id = character(0L),
      stringsAsFactors = FALSE
    )
    for (metric in names(metric_specification)) template[[metric]] <- logical(0L)
    return(list(events = template, successful_units = 0L,
                incomplete_units = 0L))
  }
  assert_true(
    all(c("unit_id", "interval_id") %in% names(observed)),
    "ALL_K aggregation requires unit and interval identifiers."
  )
  unit_groups <- if (nrow(observed)) {
    split(observed, observed$unit_id)
  } else list()
  complete_units <- list()
  incomplete_units <- 0L
  identity_fields <- intersect(
    c("unit_id", "rep_id", "seed_id", "bootstrap_seed_id"), names(observed)
  )
  for (unit in unit_groups) {
    observed_intervals <- sort(unique(as.character(unit$interval_id)))
    if (!identical(observed_intervals, expected_interval_ids) ||
        nrow(unit) != length(expected_interval_ids)) {
      incomplete_units <- incomplete_units + 1L
      next
    }
    unit_row <- unit[1L, identity_fields, drop = FALSE]
    for (metric in names(metric_specification)) {
      primitive <- unname(metric_specification[[metric]])
      assert_true(
        primitive %in% names(unit) && all(!is.na(unit[[primitive]])),
        sprintf("The primitive ALL_K event %s is absent.", primitive)
      )
      expected_value <- if (
          identical(metric,
                    "across_interval_joint_any_primary_false_report")) {
        any(as.logical(unit[[primitive]]))
      } else {
        all(as.logical(unit[[primitive]]))
      }
      if (metric %in% names(unit)) {
        assert_true(
          length(unique(unit[[metric]])) == 1L &&
            !is.na(unit[[metric]][1L]) &&
            identical(as.logical(unit[[metric]][1L]), expected_value),
          sprintf("The annotated ALL_K event %s is inconsistent.", metric)
        )
      }
      unit_row[[metric]] <- expected_value
    }
    complete_units[[length(complete_units) + 1L]] <- unit_row
  }
  deduplicated <- if (length(complete_units)) {
    do.call(rbind, complete_units)
  } else {
    template <- observed[FALSE, identity_fields, drop = FALSE]
    for (metric in names(metric_specification)) template[[metric]] <- logical(0L)
    template
  }
  list(
    events = deduplicated,
    successful_units = length(unit_groups),
    incomplete_units = incomplete_units
  )
}

# Decode the compact interval-union representation persisted by method.R.
# The empty set has the unambiguous encoding lefts="", rights="".  This
# helper is deliberately strict because paired depth diagnostics use these
# endpoints to verify nesting replication by replication.
decode_interval_union <- function(lefts, rights, tolerance = 1e-12) {
  left <- split_semicolon_numeric(lefts)
  right <- split_semicolon_numeric(rights)
  assert_true(
    length(left) == length(right),
    "A persisted interval union has unequal left/right endpoint counts."
  )
  if (!length(left)) {
    return(data.frame(left = numeric(0L), right = numeric(0L)))
  }
  assert_true(
    all(left <= right + tolerance),
    "A persisted interval union contains a reversed component."
  )
  if (length(left) >= 2L) {
    assert_true(
      all(diff(left) >= -tolerance) &&
        all(head(right, -1L) < tail(left, -1L) - tolerance),
      paste(
        "Persisted outer components must be ordered, disjoint, and",
        "already merged."
      )
    )
  }
  data.frame(left = left, right = right)
}

interval_union_contains_points <- function(union, points, tolerance = 1e-12) {
  points <- as.numeric(points)
  if (!length(points)) return(logical(0L))
  assert_true(all(is.finite(points)),
              "Interval-union containment needs finite points.")
  vapply(points, function(point) {
    nrow(union) > 0L && any(
      union$left <= point + tolerance &
        union$right >= point - tolerance
    )
  }, logical(1L))
}

interval_union_root_components <- function(
    union, roots, tolerance = 1e-12) {
  roots <- as.numeric(roots)
  assert_true(all(is.finite(roots)),
              "Root-component extraction needs finite roots.")
  if (!length(roots)) {
    return(data.frame(
      root = numeric(0L), present = logical(0L),
      left = numeric(0L), right = numeric(0L)
    ))
  }
  rows <- lapply(roots, function(root) {
    indices <- if (nrow(union)) which(
      union$left <= root + tolerance &
        union$right >= root - tolerance
    ) else integer(0L)
    assert_true(length(indices) <= 1L,
                "More than one disjoint component contains one root.")
    data.frame(
      root = root, present = length(indices) == 1L,
      left = if (length(indices)) union$left[indices] else NA_real_,
      right = if (length(indices)) union$right[indices] else NA_real_
    )
  })
  do.call(rbind, rows)
}

interval_union_is_subset <- function(
    inner, outer, tolerance = 1e-12) {
  if (!nrow(inner)) return(TRUE)
  if (!nrow(outer)) return(FALSE)
  all(vapply(seq_len(nrow(inner)), function(index) {
    any(
      outer$left <= inner$left[index] + tolerance &
        outer$right >= inner$right[index] - tolerance
    )
  }, logical(1L)))
}

persisted_interval_union_is_subset <- function(
    inner_lefts, inner_rights, outer_lefts, outer_rights,
    tolerance = 1e-12) {
  inner <- decode_interval_union(inner_lefts, inner_rights, tolerance)
  outer <- decode_interval_union(outer_lefts, outer_rights, tolerance)
  interval_union_is_subset(inner, outer, tolerance)
}

tangency_geometry_condition_masks <- function(group) {
  rows <- nrow(group)
  applicable <- if (
    "tangency_global_geometry_applicable" %in% names(group)
  ) {
    as.logical(group$tangency_global_geometry_applicable)
  } else if ("tangency_localization_applicable" %in% names(group)) {
    as.logical(group$tangency_localization_applicable)
  } else rep(FALSE, rows)
  covered <- if ("continuum_band_covers_truth" %in% names(group)) {
    as.logical(group$continuum_band_covers_truth)
  } else rep(NA, rows)
  persisted_condition <- if (
    "tangency_global_geometry_conditioned_on_continuum_band_coverage" %in%
      names(group)
  ) {
    as.logical(
      group$tangency_global_geometry_conditioned_on_continuum_band_coverage
    )
  } else applicable & !is.na(covered) & covered
  assert_true(!anyNA(applicable),
              "Tangency-geometry applicability contains a missing value.")
  expected_condition <- applicable & !is.na(covered) & covered
  assert_true(
    !anyNA(persisted_condition) &&
      identical(persisted_condition, expected_condition),
    paste(
      "The persisted tangency coverage condition disagrees with",
      "applicability and continuum-band coverage."
    )
  )
  list(
    applicable = applicable,
    continuum_covered = covered,
    condition = persisted_condition,
    applicable_denominator = sum(applicable),
    coverage_conditioned_denominator = sum(persisted_condition)
  )
}

paired_enclosure_depth_audit <- function(
    paired, reference_level, comparison_level, tolerance = 1e-10,
    decision_fields = c(
      "reversal_certified", "exact_count_certified",
      "root_isolation_success",
      "wald_reported", "wald_conditional_cover", "wald_report_and_cover",
      "root_multiplier_reported", "root_multiplier_conditional_cover",
      "root_multiplier_report_and_cover"
    )) {
  reference_level <- as.integer(reference_level)
  comparison_level <- as.integer(comparison_level)
  assert_true(
    length(reference_level) == 1L && length(comparison_level) == 1L &&
      !is.na(reference_level) && !is.na(comparison_level) &&
      reference_level != comparison_level,
    "A depth audit needs two distinct finite enclosure depths."
  )
  shallow_side <- if (reference_level < comparison_level) {
    "reference"
  } else "comparison"
  deep_side <- if (identical(shallow_side, "reference")) {
    "comparison"
  } else "reference"
  value <- function(field, side, default = NA) {
    name <- paste0(field, ".", side)
    if (name %in% names(paired)) paired[[name]] else rep(default, nrow(paired))
  }
  shallow_lefts <- value(
    "numerical_outer_component_lefts", shallow_side, ""
  )
  shallow_rights <- value(
    "numerical_outer_component_rights", shallow_side, ""
  )
  deep_lefts <- value("numerical_outer_component_lefts", deep_side, "")
  deep_rights <- value("numerical_outer_component_rights", deep_side, "")
  subset_ok <- vapply(seq_len(nrow(paired)), function(index) {
    persisted_interval_union_is_subset(
      deep_lefts[index], deep_rights[index],
      shallow_lefts[index], shallow_rights[index], tolerance
    )
  }, logical(1L))

  shallow_coverage <- as.logical(value(
    "continuum_band_covers_truth", shallow_side, NA
  ))
  deep_coverage <- as.logical(value(
    "continuum_band_covers_truth", deep_side, NA
  ))
  truth_strings <- as.character(value("true_roots", deep_side, ""))
  root_loss_eligible <- !is.na(shallow_coverage) & shallow_coverage &
    !is.na(deep_coverage) & deep_coverage & nzchar(truth_strings)
  root_lost <- rep(FALSE, nrow(paired))
  for (index in which(root_loss_eligible)) {
    roots <- split_semicolon_numeric(truth_strings[index])
    deep_union <- decode_interval_union(
      deep_lefts[index], deep_rights[index], tolerance
    )
    root_lost[index] <- !all(interval_union_contains_points(
      deep_union, roots, tolerance
    ))
  }

  shallow_length <- suppressWarnings(as.numeric(value(
    "numerical_outer_total_length", shallow_side, NA_real_
  )))
  deep_length <- suppressWarnings(as.numeric(value(
    "numerical_outer_total_length", deep_side, NA_real_
  )))
  shallow_depth_width <- suppressWarnings(as.numeric(value(
    "enclosure_depth_limit_total_width", shallow_side, NA_real_
  )))
  bound_eligible <- is.finite(shallow_length) & is.finite(deep_length) &
    is.finite(shallow_depth_width)
  length_reduction <- shallow_length - deep_length
  bound_excess <- length_reduction - shallow_depth_width
  bound_ok <- rep(NA, nrow(paired))
  bound_ok[bound_eligible] <-
    length_reduction[bound_eligible] >= -tolerance &
    bound_excess[bound_eligible] <= tolerance

  decision_comparisons <- 0L
  decision_disagreements <- 0L
  decision_disagreement_rows <- rep(FALSE, nrow(paired))
  per_decision <- character(0L)
  for (field in decision_fields) {
    shallow <- as.logical(value(field, shallow_side, NA))
    deep <- as.logical(value(field, deep_side, NA))
    keep <- !is.na(shallow) & !is.na(deep)
    disagreements <- keep & shallow != deep
    decision_comparisons <- decision_comparisons + sum(keep)
    decision_disagreements <- decision_disagreements + sum(disagreements)
    decision_disagreement_rows <- decision_disagreement_rows | disagreements
    per_decision <- c(
      per_decision,
      sprintf("%s=%d/%d", field, sum(disagreements), sum(keep))
    )
  }
  list(
    shallow_side = shallow_side, deep_side = deep_side,
    shallow_level = min(reference_level, comparison_level),
    deep_level = max(reference_level, comparison_level),
    paired_replications = nrow(paired),
    subset_checks = length(subset_ok), subset_failures = sum(!subset_ok),
    root_retention_checks = sum(root_loss_eligible),
    root_losses_under_continuum_coverage = sum(root_lost[root_loss_eligible]),
    length_bound_checks = sum(bound_eligible),
    length_bound_failures = sum(!bound_ok[bound_eligible]),
    maximum_length_bound_excess = if (any(bound_eligible)) {
      max(bound_excess[bound_eligible])
    } else NA_real_,
    minimum_length_reduction = if (any(bound_eligible)) {
      min(length_reduction[bound_eligible])
    } else NA_real_,
    decision_comparisons = decision_comparisons,
    decision_disagreements = decision_disagreements,
    decision_disagreement_replications = sum(decision_disagreement_rows),
    decision_disagreement_detail = paste(per_decision, collapse = ";")
  )
}
