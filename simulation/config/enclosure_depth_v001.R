local({
  intervals <- list(
    K1 = c(0.25, 1.50),
    K2 = c(0.25, 2.00),
    K3 = c(0.25, 2.50)
  )
  q_two_root <- list(absolute = list(
    Q1 = c(0.45, 1.05),
    Q2 = c(1.40, 2.10)
  ))
  no_regular_root <- list(absolute = list())

  # Both variants reproduce the bridge.6 selected scientific configuration.
  # The sole changed field is the dedicated maximum continuum-enclosure
  # depth.  The relative-variance guard remains fixed at legacy depth six.
  variants <- list(
    depth6_reference = list(
      bootstrap_reps = 999L,
      bootstrap_grid_spacing = 0.005,
      grid_size = 451L,
      enclosure_grid_spacing = 0.0025,
      enclosure_grid_size = 901L,
      maximum_enclosure_half_width = 0.0025,
      maximum_enclosure_bisections = 6L,
      maximum_continuum_enclosure_bisections = 6L,
      maximum_enclosure_nodes = 57601L
    ),
    depth8_refined = list(
      bootstrap_reps = 999L,
      bootstrap_grid_spacing = 0.005,
      grid_size = 451L,
      enclosure_grid_spacing = 0.0025,
      enclosure_grid_size = 901L,
      maximum_enclosure_half_width = 0.0025,
      maximum_enclosure_bisections = 6L,
      maximum_continuum_enclosure_bisections = 8L,
      maximum_enclosure_nodes = 57601L
    )
  )

  make_spec <- function(
      spec_id, scenario_id, n, interval_id, root_intervals,
      run_wald_refinement, run_root_multiplier_refinement) {
    list(
      spec_id = spec_id,
      scenario_id = scenario_id,
      n = as.integer(n),
      interval_id = interval_id,
      root_intervals = root_intervals,
      run_wald_refinement = run_wald_refinement,
      run_root_multiplier_refinement = run_root_multiplier_refinement
    )
  }

  specs <- list()
  for (n in c(2000L, 4000L, 8000L)) {
    specs[[length(specs) + 1L]] <- make_spec(
      sprintf("TANGENCY_SV_POWER_n%d", n),
      "TANGENCY_SV_POWER", n, "K1", no_regular_root, FALSE, FALSE
    )
  }
  for (n in c(2000L, 4000L, 8000L)) {
    specs[[length(specs) + 1L]] <- make_spec(
      sprintf("TWO_ROOT_SV_POWER_n%d", n),
      "TWO_ROOT_SV_POWER", n, "K3", q_two_root, TRUE, TRUE
    )
  }

  cells <- list()
  for (spec in specs) {
    pair_key <- sprintf("ED1_%s", spec$spec_id)
    for (variant_id in names(variants)) {
      cells[[length(cells) + 1L]] <- list(
        cell_id = sprintf(
          "ED1_%03d_%s_%s", length(cells) + 1L,
          spec$spec_id, variant_id
        ),
        scenario_id = spec$scenario_id,
        n_x = spec$n,
        n_y = spec$n,
        moment_types = "absolute",
        interval_ids = spec$interval_id,
        variant_id = variant_id,
        cell_class = "diagnostic",
        pairing_key = pair_key,
        multiplier_pairing_key = pair_key,
        root_intervals = spec$root_intervals,
        run_wald_refinement = spec$run_wald_refinement,
        run_root_multiplier_refinement =
          spec$run_root_multiplier_refinement,
        run_derivative_assisted_ablation = FALSE
      )
    }
  }

  variant_ids <- vapply(cells, `[[`, character(1L), "variant_id")
  pairing_keys <- vapply(cells, `[[`, character(1L), "pairing_key")
  tangency_cells <- cells[
    vapply(cells, `[[`, character(1L), "scenario_id") ==
      "TANGENCY_SV_POWER"
  ]
  two_root_cells <- cells[
    vapply(cells, `[[`, character(1L), "scenario_id") ==
      "TWO_ROOT_SV_POWER"
  ]
  stopifnot(
    length(specs) == 6L,
    length(cells) == 12L,
    length(unique(vapply(cells, `[[`, character(1L), "cell_id"))) == 12L,
    identical(
      as.integer(table(factor(variant_ids, levels = names(variants)))),
      c(6L, 6L)
    ),
    length(unique(pairing_keys)) == 6L,
    all(table(pairing_keys) == 2L),
    identical(
      pairing_keys,
      vapply(cells, `[[`, character(1L), "multiplier_pairing_key")
    ),
    identical(
      sort(unique(vapply(tangency_cells, `[[`, integer(1L), "n_x"))),
      c(2000L, 4000L, 8000L)
    ),
    identical(
      sort(unique(vapply(two_root_cells, `[[`, integer(1L), "n_x"))),
      c(2000L, 4000L, 8000L)
    ),
    all(vapply(tangency_cells, `[[`, character(1L), "interval_ids") ==
          "K1"),
    all(vapply(two_root_cells, `[[`, character(1L), "interval_ids") ==
          "K3"),
    all(vapply(
      tangency_cells,
      function(cell) !cell$run_wald_refinement &&
        !cell$run_root_multiplier_refinement &&
        identical(cell$root_intervals, no_regular_root),
      logical(1L)
    )),
    all(vapply(
      two_root_cells,
      function(cell) cell$run_wald_refinement &&
        cell$run_root_multiplier_refinement &&
        identical(cell$root_intervals, q_two_root),
      logical(1L)
    )),
    variants$depth6_reference$maximum_enclosure_bisections == 6L,
    variants$depth8_refined$maximum_enclosure_bisections == 6L,
    variants$depth6_reference$maximum_continuum_enclosure_bisections == 6L,
    variants$depth8_refined$maximum_continuum_enclosure_bisections == 8L,
    identical(
      variants$depth6_reference[names(variants$depth6_reference) !=
        "maximum_continuum_enclosure_bisections"],
      variants$depth8_refined[names(variants$depth8_refined) !=
        "maximum_continuum_enclosure_bisections"]
    ),
    !any(vapply(
      cells, `[[`, logical(1L), "run_derivative_assisted_ablation"
    )),
    length(cells) * 254L == 3048L
  )

  paired_contrasts <- data.frame(
    contrast_id = "depth8_refined_minus_depth6_reference",
    contrast_definition = paste(
      "maximum continuum-enclosure depth 8 minus 6, with guard depth fixed",
      "at 6 and with identical samples,",
      "multipliers, B=999, bootstrap grid 0.005, enclosure grid 0.0025,",
      "node ceiling 57601, and every remaining inferential setting"
    ),
    reference_variant_id = "depth6_reference",
    comparison_variant_id = "depth8_refined",
    stringsAsFactors = FALSE
  )

  list(
    schema_version = "2.0.0",
    study_id = "enclosure_depth_v001",
    release_status = "ready",
    design_locked = TRUE,
    study_role = "fully_paired_enclosure_depth_diagnostic",
    purpose = paste(
      "Fully paired one-factor diagnostic of maximum continuum-enclosure",
      "depth 6 versus 8 for both signal-variance-selected DGPs and all",
      "calibrated sample sizes; the diagnostic measures the severity, not",
      "merely the incidence, of statistically unresolved cells at the",
      "configured maximum level"
    ),
    master_seed = 2026082012L,
    cells = cells,
    variants = variants,
    paired_contrasts = paired_contrasts,
    diagnostic_targets = list(
      tangency = c(
        "continuum_band_coverage", "outer_set_truth_coverage",
        "tangency_global_outer_total_length_conditional_on_continuum_band_coverage",
        "tangency_global_outer_hausdorff_conditional_on_continuum_band_coverage",
        "tangency_global_outer_component_count_conditional_on_continuum_band_coverage",
        "tangency_global_root_component_length_conditional_on_continuum_band_coverage"
      ),
      two_root = c(
        "joint_root_isolation", "exact_count_certification",
        "report_rate", "conditional_coverage", "report_and_cover",
        "global_outer_total_length"
      ),
      enclosure_severity = c(
        "enclosure_depth_limit_cells",
        "enclosure_depth_limit_total_width",
        "enclosure_depth_limit_width_proportion",
        "enclosure_depth_limit_maximum_cell_width",
        "enclosure_variance_limit_cells", "enclosure_nodes_used",
        "enclosure_maximum_level"
      )
    ),
    interpretation = list(
      tangency_local_window = paste(
        "the fixed [0.35,0.85] summaries are truncated diagnostics and are",
        "not the primary contraction estimands"
      ),
      tangency_global_outer = paste(
        "global K1 outer length, Hausdorff distance and root-component",
        "geometry are the primary n^{-1/4} diagnostics"
      ),
      limit_hit = paste(
        "an any-cell limit indicator is an incidence flag; scientific",
        "severity is determined by unresolved-cell counts and widths"
      ),
      decision = paste(
        "manual scientific and numerical audit is required before freezing",
        "or submitting main_final_v003"
      )
    ),
    expected_output_contract = list(
      cells = 12L,
      units = 3048L,
      tasks = 12L,
      independent_sample_streams = 1524L,
      independent_bootstrap_streams = 1524L,
      scientific_groups = 12L,
      replication_rows = 3048L,
      binary_metric_rows = 1284L,
      continuous_metric_rows = 1464L,
      paired_comparison_rows = 1374L,
      enclosure_depth_audit_rows = 6L,
      across_interval_groups = 12L,
      across_interval_binary_metric_rows = 48L,
      root_metric_rows = 12L
    ),
    scenario_ids = c("TANGENCY_SV_POWER", "TWO_ROOT_SV_POWER"),
    sample_sizes = list(
      c(n_x = 2000L, n_y = 2000L),
      c(n_x = 4000L, n_y = 4000L),
      c(n_x = 8000L, n_y = 8000L)
    ),
    moment_types = "absolute",
    order_intervals = intervals,
    interval_calibration = "separate_restricted_suprema_from_common_process",
    root_interval_policy = list(
      lognormal_p4_over_3 = c(0.75, 1.90),
      selected_tangency_root = 0.60,
      selected_two_root = list(c(0.45, 1.05), c(1.40, 2.10))
    ),
    root_reporting_definitions = list(
      root_isolation_report_rate =
        "fraction_of_all_applicable_successful_replications_with_prespecified_root_isolation",
      report_rate =
        "fraction_of_all_applicable_successful_replications_with_required_intervals_reported",
      conditional_coverage =
        "coverage_conditional_on_required_intervals_reported",
      report_and_cover =
        "fraction_of_all_applicable_successful_replications_reporting_and_covering",
      root_bias_conditional_on_isolation =
        "conditional_on_root_isolation",
      root_rmse_conditional_on_isolation =
        "conditional_on_root_isolation",
      root_collection_bias_conditional_on_joint_isolation =
        "pooled_root_errors_only_after_joint_isolation_of_every_required_root",
      root_collection_rmse_conditional_on_joint_isolation =
        "sqrt_of_pooled_squared_root_errors_after_joint_isolation"
    ),
    p_min = 0.25,
    p_max = 2.50,
    bootstrap_grid_spacing = 0.005,
    enclosure_grid_spacing = 0.0025,
    grid_size = 451L,
    enclosure_grid_size = 901L,
    maximum_enclosure_half_width = 0.0025,
    audit_grid_spacing = 0.00125,
    audit_grid_size = 1801L,
    maximum_enclosure_bisections = 6L,
    maximum_continuum_enclosure_bisections = 6L,
    maximum_enclosure_nodes = 57601L,
    maximum_truth_bisection_levels = 20L,
    maximum_truth_subintervals = 50000L,
    enclosure_absolute_tolerance = 1e-12,
    roundoff_inflation_required = TRUE,
    root_interval = c(0.75, 1.90),
    tangency_neighborhood = c(0.35, 0.85),
    two_root_intervals = list(c(0.45, 1.05), c(1.40, 2.10)),
    root_scan_size = 2001L,
    maximum_root_bisections = 24L,
    maximum_root_evaluations = 50000L,
    root_absolute_tolerance = 1e-10,
    alpha = 0.05,
    bootstrap_reps = 999L,
    multiplier_distribution = "rademacher",
    bootstrap_batch_size = 50L,
    variance_tolerance = 1e-12,
    slope_tolerance = 1e-10,
    contrast_variance_tolerance = 1e-12,
    root_variance_tolerance = 1e-12,
    augmented_variance_tolerance = 1e-12,
    numerical_tolerance = 1e-12,
    run_derivative_assisted_ablation = FALSE,
    run_grid_only_ablation = FALSE,
    run_pointwise_ablation = FALSE,
    run_wald_refinement = FALSE,
    run_root_multiplier_refinement = FALSE,
    formal_wald_guard_set = "R",
    tail_diagnostic_orders = c(0.60, 0.75, 1.50, 1.75, 2.50),
    tail_diagnostic_probabilities = c(0.95, 0.99),
    tail_top_fraction = 0.01,
    mc_reps_per_cell = 254L,
    task_packing = "global",
    total_units = 3048L,
    units_per_task = 254L,
    units_per_worker_chunk = 1L,
    maximum_workers = 127L,
    checkpoint_reserve_seconds = 180L,
    internal_wall_budget_seconds = 13500L,
    confirmatory_main_seed_reserved = 2026082008L
  )
})
