local({
  intervals <- list(
    K1 = c(0.25, 1.50),
    K2 = c(0.25, 2.00),
    K3 = c(0.25, 2.50)
  )
  q_lognormal <- list(absolute = list(Q_LN = c(0.75, 1.90)))
  q_p1 <- list(absolute = list(Q_P1 = c(0.50, 1.45)))
  q_two_root <- list(absolute = list(
    Q1 = c(0.50, 1.45),
    Q2 = c(1.65, 2.35)
  ))
  q_two_root_selected <- list(absolute = list(
    Q1 = c(0.45, 1.05),
    Q2 = c(1.40, 2.10)
  ))
  no_regular_root <- list(absolute = list())
  mixed_sign_roots <- list(
    absolute = list(),
    signed = list(Q_P1 = c(0.50, 1.45))
  )

  make_cell <- function(
      cell_id, scenario_id, n_x, n_y = n_x,
      cell_class, interval_ids = names(intervals),
      moment_types = "absolute", root_intervals = no_regular_root,
      run_wald_refinement = FALSE,
      run_root_multiplier_refinement = FALSE, ...) {
    base <- list(
      cell_id = cell_id,
      scenario_id = scenario_id,
      n_x = as.integer(n_x),
      n_y = as.integer(n_y),
      moment_types = moment_types,
      interval_ids = interval_ids,
      variant_id = "main_primary",
      cell_class = cell_class,
      pairing_key = sprintf("%s_%s", cell_id, scenario_id),
      multiplier_pairing_key = sprintf("%s_%s", cell_id, scenario_id),
      root_intervals = root_intervals,
      run_wald_refinement = run_wald_refinement,
      run_root_multiplier_refinement = run_root_multiplier_refinement,
      run_derivative_assisted_ablation = FALSE
    )
    overrides <- list(...)
    stopifnot(
      is.list(overrides),
      is.null(names(overrides)) || all(nzchar(names(overrides))),
      !any(names(overrides) %in% names(base))
    )
    c(base, overrides)
  }

  cells <- list(
    # Benchmark cells: null/equality calibration and signed-moment interface.
    make_cell("M001_LN_EQUAL_n50", "LN_EQUAL", 50L,
              cell_class = "benchmark"),
    make_cell("M002_LN_EQUAL_n100", "LN_EQUAL", 100L,
              cell_class = "benchmark"),
    make_cell("M003_LN_EQUAL_n250", "LN_EQUAL", 250L,
              cell_class = "benchmark"),
    make_cell("M004_LN_EQUAL_n500", "LN_EQUAL", 500L,
              cell_class = "benchmark"),
    make_cell("M005_LN_NO_ROOT_n50", "LN_NO_ROOT", 50L,
              cell_class = "benchmark"),
    make_cell("M006_LN_NO_ROOT_n100", "LN_NO_ROOT", 100L,
              cell_class = "benchmark"),
    make_cell("M007_LN_NO_ROOT_n250", "LN_NO_ROOT", 250L,
              cell_class = "benchmark"),
    make_cell("M008_LN_NO_ROOT_n500", "LN_NO_ROOT", 500L,
              cell_class = "benchmark"),
    make_cell(
      "M009_MIXED_SIGN_n100", "MIXED_SIGN", 100L,
      cell_class = "benchmark", moment_types = c("absolute", "signed"),
      root_intervals = mixed_sign_roots,
      run_wald_refinement = TRUE,
      run_root_multiplier_refinement = TRUE
    ),
    make_cell(
      "M010_MIXED_SIGN_n250", "MIXED_SIGN", 250L,
      cell_class = "benchmark", moment_types = c("absolute", "signed"),
      root_intervals = mixed_sign_roots,
      run_wald_refinement = TRUE,
      run_root_multiplier_refinement = TRUE
    ),

    # Power cells: regular alternatives with prespecified simple roots.
    make_cell(
      "M011_LN_MODERATE_n100", "LN_MODERATE", 100L,
      cell_class = "power", root_intervals = q_lognormal,
      run_wald_refinement = TRUE, run_root_multiplier_refinement = TRUE
    ),
    make_cell(
      "M012_LN_MODERATE_n250", "LN_MODERATE", 250L,
      cell_class = "power", root_intervals = q_lognormal,
      run_wald_refinement = TRUE, run_root_multiplier_refinement = TRUE
    ),
    make_cell(
      "M013_LN_MODERATE_n500", "LN_MODERATE", 500L,
      cell_class = "power", root_intervals = q_lognormal,
      run_wald_refinement = TRUE, run_root_multiplier_refinement = TRUE
    ),
    make_cell(
      "M014_LN_MODERATE_n1000", "LN_MODERATE", 1000L,
      cell_class = "power", root_intervals = q_lognormal,
      run_wald_refinement = TRUE, run_root_multiplier_refinement = TRUE
    ),
    make_cell(
      "M015_LN_WEAK_n100", "LN_WEAK", 100L,
      cell_class = "power", root_intervals = q_lognormal,
      run_wald_refinement = TRUE, run_root_multiplier_refinement = TRUE
    ),
    make_cell(
      "M016_LN_WEAK_n250", "LN_WEAK", 250L,
      cell_class = "power", root_intervals = q_lognormal,
      run_wald_refinement = TRUE, run_root_multiplier_refinement = TRUE
    ),
    make_cell(
      "M017_LN_WEAK_n500", "LN_WEAK", 500L,
      cell_class = "power", root_intervals = q_lognormal,
      run_wald_refinement = TRUE, run_root_multiplier_refinement = TRUE
    ),
    make_cell(
      "M018_LN_WEAK_n1000", "LN_WEAK", 1000L,
      cell_class = "power", root_intervals = q_lognormal,
      run_wald_refinement = TRUE, run_root_multiplier_refinement = TRUE
    ),
    make_cell(
      "M019_RARE_CROSSING_P05_n100", "RARE_CROSSING_P05", 100L,
      cell_class = "power", root_intervals = q_p1,
      run_wald_refinement = TRUE, run_root_multiplier_refinement = TRUE
    ),
    make_cell(
      "M020_RARE_CROSSING_P05_n250", "RARE_CROSSING_P05", 250L,
      cell_class = "power", root_intervals = q_p1,
      run_wald_refinement = TRUE, run_root_multiplier_refinement = TRUE
    ),
    make_cell(
      "M021_RARE_CROSSING_P05_n500", "RARE_CROSSING_P05", 500L,
      cell_class = "power", root_intervals = q_p1,
      run_wald_refinement = TRUE, run_root_multiplier_refinement = TRUE
    ),
    make_cell(
      "M022_RARE_CROSSING_P05_n1000", "RARE_CROSSING_P05", 1000L,
      cell_class = "power", root_intervals = q_p1,
      run_wald_refinement = TRUE, run_root_multiplier_refinement = TRUE
    ),
    make_cell(
      "M023_TANGENCY_STRONG_n500", "TANGENCY_STRONG", 500L,
      cell_class = "stress", interval_ids = "K1",
      tangency_neighborhood = c(0.50, 1.45)
    ),
    make_cell(
      "M024_TANGENCY_STRONG_n1000", "TANGENCY_STRONG", 1000L,
      cell_class = "stress", interval_ids = "K1",
      tangency_neighborhood = c(0.50, 1.45)
    ),
    make_cell(
      "M025_TWO_ROOT_STRONG_n500", "TWO_ROOT_STRONG", 500L,
      cell_class = "stress", interval_ids = "K3",
      root_intervals = q_two_root, run_wald_refinement = TRUE,
      run_root_multiplier_refinement = TRUE
    ),
    make_cell(
      "M026_TWO_ROOT_STRONG_n1000", "TWO_ROOT_STRONG", 1000L,
      cell_class = "stress", interval_ids = "K3",
      root_intervals = q_two_root, run_wald_refinement = TRUE,
      run_root_multiplier_refinement = TRUE
    ),

    # Stress cells: reported in full, but 95% finite-sample coverage is not an
    # acceptance criterion for this class.
    make_cell(
      "M027_LN_HIGH_VARIANCE_n250", "LN_HIGH_VARIANCE", 250L,
      cell_class = "stress", root_intervals = q_lognormal,
      run_wald_refinement = TRUE, run_root_multiplier_refinement = TRUE
    ),
    make_cell(
      "M028_LN_HIGH_VARIANCE_n500", "LN_HIGH_VARIANCE", 500L,
      cell_class = "stress", root_intervals = q_lognormal,
      run_wald_refinement = TRUE, run_root_multiplier_refinement = TRUE
    ),
    make_cell(
      "M029_LN_HIGH_VARIANCE_n1000", "LN_HIGH_VARIANCE", 1000L,
      cell_class = "stress", root_intervals = q_lognormal,
      run_wald_refinement = TRUE, run_root_multiplier_refinement = TRUE
    ),
    make_cell("M030_RARE_EQUAL_P01_n100", "RARE_EQUAL_P01", 100L,
              cell_class = "stress"),
    make_cell("M031_RARE_EQUAL_P01_n250", "RARE_EQUAL_P01", 250L,
              cell_class = "stress"),
    make_cell("M032_RARE_EQUAL_P01_n500", "RARE_EQUAL_P01", 500L,
              cell_class = "stress"),
    make_cell("M033_RARE_NO_ROOT_P01_n100", "RARE_NO_ROOT_P01", 100L,
              cell_class = "stress"),
    make_cell("M034_RARE_NO_ROOT_P01_n250", "RARE_NO_ROOT_P01", 250L,
              cell_class = "stress"),
    make_cell("M035_RARE_NO_ROOT_P01_n500", "RARE_NO_ROOT_P01", 500L,
              cell_class = "stress"),
    make_cell(
      "M036_RARE_CROSSING_P01_n100", "RARE_CROSSING_P01", 100L,
      cell_class = "stress", root_intervals = q_p1,
      run_wald_refinement = TRUE, run_root_multiplier_refinement = TRUE
    ),
    make_cell(
      "M037_RARE_CROSSING_P01_n250", "RARE_CROSSING_P01", 250L,
      cell_class = "stress", root_intervals = q_p1,
      run_wald_refinement = TRUE, run_root_multiplier_refinement = TRUE
    ),
    make_cell(
      "M038_RARE_CROSSING_P01_n500", "RARE_CROSSING_P01", 500L,
      cell_class = "stress", root_intervals = q_p1,
      run_wald_refinement = TRUE, run_root_multiplier_refinement = TRUE
    ),
    make_cell(
      "M039_RARE_CROSSING_P01_n1000", "RARE_CROSSING_P01", 1000L,
      cell_class = "stress", root_intervals = q_p1,
      run_wald_refinement = TRUE, run_root_multiplier_refinement = TRUE
    ),
    make_cell("M040_TANGENCY_n500", "TANGENCY", 500L,
              cell_class = "stress", interval_ids = "K1",
              tangency_neighborhood = c(0.50, 1.45)),
    make_cell("M041_TANGENCY_n1000", "TANGENCY", 1000L,
              cell_class = "stress", interval_ids = "K1",
              tangency_neighborhood = c(0.50, 1.45)),
    make_cell(
      "M042_TWO_ROOT_n500", "TWO_ROOT", 500L,
      cell_class = "stress", interval_ids = "K3",
      root_intervals = q_two_root, run_wald_refinement = TRUE,
      run_root_multiplier_refinement = TRUE
    ),
    make_cell(
      "M043_TWO_ROOT_n1000", "TWO_ROOT", 1000L,
      cell_class = "stress", interval_ids = "K3",
      root_intervals = q_two_root, run_wald_refinement = TRUE,
      run_root_multiplier_refinement = TRUE
    ),
    make_cell(
      "M044_PARETO_n250", "PARETO", 250L,
      cell_class = "stress", root_intervals = q_p1,
      run_wald_refinement = TRUE, run_root_multiplier_refinement = TRUE
    ),
    make_cell(
      "M045_PARETO_n500", "PARETO", 500L,
      cell_class = "stress", root_intervals = q_p1,
      run_wald_refinement = TRUE, run_root_multiplier_refinement = TRUE
    ),
    make_cell(
      "M046_PARETO_n1000", "PARETO", 1000L,
      cell_class = "stress", root_intervals = q_p1,
      run_wald_refinement = TRUE, run_root_multiplier_refinement = TRUE
    ),
    make_cell(
      "M047_RARE_CROSSING_P01_nx100_ny400", "RARE_CROSSING_P01",
      100L, 400L, "stress", root_intervals = q_p1,
      run_wald_refinement = TRUE, run_root_multiplier_refinement = TRUE
    ),
    make_cell(
      "M048_RARE_CROSSING_P01_nx400_ny100", "RARE_CROSSING_P01",
      400L, 100L, "stress", root_intervals = q_p1,
      run_wald_refinement = TRUE, run_root_multiplier_refinement = TRUE
    ),

    # Confirmatory power cells selected analytically before Monte Carlo.
    make_cell(
      "M049_TANGENCY_SV_POWER_n2000", "TANGENCY_SV_POWER", 2000L,
      cell_class = "power", interval_ids = "K1",
      tangency_neighborhood = c(0.35, 0.85)
    ),
    make_cell(
      "M050_TANGENCY_SV_POWER_n4000", "TANGENCY_SV_POWER", 4000L,
      cell_class = "power", interval_ids = "K1",
      tangency_neighborhood = c(0.35, 0.85)
    ),
    make_cell(
      "M051_TANGENCY_SV_POWER_n8000", "TANGENCY_SV_POWER", 8000L,
      cell_class = "power", interval_ids = "K1",
      tangency_neighborhood = c(0.35, 0.85)
    ),
    make_cell(
      "M052_TWO_ROOT_SV_POWER_n2000", "TWO_ROOT_SV_POWER", 2000L,
      cell_class = "power", interval_ids = "K3",
      root_intervals = q_two_root_selected, run_wald_refinement = TRUE,
      run_root_multiplier_refinement = TRUE
    ),
    make_cell(
      "M053_TWO_ROOT_SV_POWER_n4000", "TWO_ROOT_SV_POWER", 4000L,
      cell_class = "power", interval_ids = "K3",
      root_intervals = q_two_root_selected, run_wald_refinement = TRUE,
      run_root_multiplier_refinement = TRUE
    ),
    make_cell(
      "M054_TWO_ROOT_SV_POWER_n8000", "TWO_ROOT_SV_POWER", 8000L,
      cell_class = "power", interval_ids = "K3",
      root_intervals = q_two_root_selected, run_wald_refinement = TRUE,
      run_root_multiplier_refinement = TRUE
    )
  )

  cell_ids <- vapply(cells, `[[`, character(1L), "cell_id")
  class_counts <- table(vapply(cells, `[[`, character(1L), "cell_class"))
  stopifnot(
    length(cells) == 54L,
    length(unique(cell_ids)) == 54L,
    identical(cell_ids, sprintf("M%03d_%s", seq_len(54L), sub(
      "^M[0-9]{3}_", "", cell_ids
    ))),
    identical(
      as.integer(class_counts[c("benchmark", "power", "stress")]),
      c(10L, 18L, 26L)
    ),
    all(vapply(cells, `[[`, character(1L), "variant_id") == "main_primary"),
    !any(vapply(
      cells, `[[`, logical(1L), "run_derivative_assisted_ablation"
    )),
    sum(vapply(cells, function(cell) {
      cell$scenario_id == "TANGENCY"
    }, logical(1L))) == 2L,
    sum(vapply(cells, function(cell) {
      cell$scenario_id == "TANGENCY_STRONG" && cell$cell_class == "stress"
    }, logical(1L))) == 2L,
    sum(vapply(cells, function(cell) {
      cell$scenario_id == "TWO_ROOT"
    }, logical(1L))) == 2L,
    sum(vapply(cells, function(cell) {
      cell$scenario_id == "TWO_ROOT_STRONG" && cell$cell_class == "stress"
    }, logical(1L))) == 2L
  )

  list(
    schema_version = "2.0.0",
    study_id = "main_final_v003",
    release_status = "ready",
    design_locked = TRUE,
    study_role = "independent_confirmatory_main",
    freeze_decision = paste(
      "frozen after the analytic signal-variance screen, independent power",
      "and null diagnostics, fully paired sensitivity study, and the paired",
      "depth-six versus depth-eight enclosure diagnostic"
    ),
    purpose = paste(
      "Independent main experiment with prespecified benchmark, power, and",
      "stress classes; selected tangency and two-root alternatives provide",
      "confirmatory power cells, while difficult mechanisms are retained and",
      "reported without imposing nominal-coverage success on stress cells"
    ),
    master_seed = 2026082008L,
    cells = cells,
    scenario_ids = unique(vapply(cells, `[[`, character(1L), "scenario_id")),
    sample_sizes = unique(lapply(cells, function(cell) {
      c(n_x = cell$n_x, n_y = cell$n_y)
    })),
    moment_types = c("absolute", "signed"),
    order_intervals = intervals,
    interval_calibration = "separate_restricted_suprema_from_common_process",
    cell_class_policy = list(
      benchmark = paste(
        "finite-sample calibration and implementation benchmarks; evaluate",
        "coverage and false-report control against prespecified Monte Carlo",
        "uncertainty"
      ),
      power = paste(
        "evaluate report, certification, and report-and-cover probabilities",
        "with no conditioning hidden in unconditional metrics"
      ),
      stress = paste(
        "retain and report all results and latent-tail strata; nominal 95%",
        "coverage is not an acceptance requirement"
      )
    ),
    root_interval_policy = list(
      lognormal_p4_over_3 = c(0.75, 1.90),
      simple_root_p1 = c(0.50, 1.45),
      two_root = list(c(0.50, 1.45), c(1.65, 2.35)),
      selected_tangency_root = 0.60,
      selected_two_root = list(c(0.45, 1.05), c(1.40, 2.10)),
      tangencies = "outer_set_only_no_regular_root_interval"
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
      root_bias_conditional_on_isolation = "conditional_on_root_isolation",
      root_rmse_conditional_on_isolation = "conditional_on_root_isolation",
      root_collection_bias_conditional_on_joint_isolation =
        "pooled_root_errors_only_after_joint_isolation_of_every_required_root",
      root_collection_rmse_conditional_on_joint_isolation =
        "sqrt_of_pooled_squared_root_errors_after_joint_isolation"
    ),
    expected_output_contract = list(
      cells = 54L,
      units = 109728L,
      tasks = 48L,
      independent_sample_streams = 109728L,
      independent_bootstrap_streams = 109728L,
      scientific_groups = 140L,
      replication_rows = 284480L,
      binary_metric_rows = 14980L,
      continuous_metric_rows = 17080L,
      paired_comparison_rows = 0L,
      enclosure_depth_audit_rows = 0L,
      across_interval_groups = 56L,
      across_interval_binary_metric_rows = 224L,
      root_metric_rows = 81L,
      root_collection_deduplicated_rows = 33L,
      root_metric_deduplicated_rows = 40L
    ),
    article_output_contract = list(
      design_rows = 17L,
      benchmark_editorial_rows = 12L,
      root_collection_rows = 33L,
      root_specific_rows = 40L,
      special_power_rows = 6L,
      stress_scientific_groups = 62L
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
    maximum_continuum_enclosure_bisections = 8L,
    maximum_enclosure_nodes = 57601L,
    maximum_truth_bisection_levels = 20L,
    maximum_truth_subintervals = 50000L,
    enclosure_absolute_tolerance = 1e-12,
    roundoff_inflation_required = TRUE,
    root_interval = c(0.75, 1.90),
    tangency_neighborhood = c(0.50, 1.45),
    two_root_intervals = list(c(0.50, 1.45), c(1.65, 2.35)),
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
    run_grid_only_ablation = TRUE,
    run_pointwise_ablation = TRUE,
    run_wald_refinement = FALSE,
    run_root_multiplier_refinement = FALSE,
    formal_wald_guard_set = "R",
    tail_diagnostic_orders = c(0.60, 0.75, 1.00, 1.50, 1.75, 2.00, 2.50),
    tail_diagnostic_probabilities = c(0.95, 0.99),
    tail_top_fraction = 0.01,
    mc_reps_per_cell = 2032L,
    task_packing = "global",
    total_units = 109728L,
    units_per_task = 2286L,
    units_per_worker_chunk = 2L,
    maximum_workers = 127L,
    checkpoint_reserve_seconds = 180L,
    internal_wall_budget_seconds = 13500L,
    article_output_profile = "main_final_v003"
  )
})
