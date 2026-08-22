list(
  application_id = "hillstrom_application_v001",
  application_version = "1.0.2",
  design_status = "locked_before_canonical_run",
  engine = list(
    version = "2.0.0-bridge.8.1",
    manifest_sha256 = "5cdf653d8db1ab2ed59f9f1a02fb83459bf75bd3d042939ddb76fcd755da8182",
    common_r_sha256 = "98f66209fa0e90834741261fe7d9597fe534335a055463e029e154381a074cf5",
    method_r_sha256 = "d12fae4433ba57efb3d1d49a86ea2abcdd35f971b01920bb8d057b41f0efde63"
  ),
  data = list(
    title = "MineThatData E-Mail Analytics and Data Mining Challenge",
    source_page = "https://blog.minethatdata.com/2008/03/minethatdata-e-mail-analytics-and-data.html",
    urls = c(
      "https://www.minethatdata.com/Kevin_Hillstrom_MineThatData_E-MailAnalytics_DataMiningChallenge_2008.03.20.csv"
    ),
    raw_filename = "hillstrom_2008_raw.csv",
    sha256 = "0e5893329d8b93cefecc571777672028290ab69865718020c78c7284f291aece",
    bytes = 3964977,
    rows = 64000L,
    expected_columns = c(
      "recency", "history_segment", "history", "mens", "womens",
      "zip_code", "newbie", "channel", "segment", "visit",
      "conversion", "spend"
    ),
    segment_counts = c(
      "Mens E-Mail" = 21307L,
      "Womens E-Mail" = 21387L,
      "No E-Mail" = 21306L
    ),
    positive_spend_counts = c(
      "Mens E-Mail" = 267L,
      "Womens E-Mail" = 189L,
      "No E-Mail" = 122L
    ),
    observed_spend_range = c(0, 499),
    recency_range = c(1L, 12L),
    category_levels = list(
      history_segment = c(
        "1) $0 - $100", "2) $100 - $200", "3) $200 - $350",
        "4) $350 - $500", "5) $500 - $750", "6) $750 - $1,000",
        "7) $1,000 +"
      ),
      zip_code = c("Rural", "Surburban", "Urban"),
      channel = c("Phone", "Web", "Multichannel"),
      segment = c("Mens E-Mail", "Womens E-Mail", "No E-Mail")
    ),
    redistribution = "not_included_no_formal_upstream_license"
  ),
  analysis = list(
    x_segment = "No E-Mail",
    y_segment = "Mens E-Mail",
    contrast = "Y_minus_X",
    outcome = "spend",
    moment_types = "absolute",
    order_intervals = list(K1 = c(0.25, 1.50)),
    alpha = 0.05,
    bootstrap_reps = 9999L,
    bootstrap_seed = 20260822L,
    rng_kind = c("Mersenne-Twister", "Inversion", "Rejection"),
    multiplier_distribution = "rademacher",
    bootstrap_batch_size = 50L,
    bootstrap_grid_spacing = 0.005,
    enclosure_grid_spacing = 0.0025,
    audit_grid_spacing = 0.00125,
    expected_bootstrap_grid_size = 251L,
    expected_enclosure_grid_size = 501L,
    expected_audit_grid_size = 1001L,
    variance_tolerance = 1e-12,
    guard_max_levels = 6L,
    enclosure_max_levels = 8L,
    enclosure_max_nodes = 57601L,
    enclosure_refine_unresolved_levels = 0L,
    enclosure_safety_margin = 1e-12,
    maximum_enclosure_half_width = 0.0025,
    numerical_tolerance = 1e-12,
    roundoff_inflation = TRUE,
    evaluation_chunk_size = 256L,
    anchor_orders = c(0.25, 0.50, 0.75, 1.00, 1.25, 1.50),
    tail_diagnostic_orders = c(0.25, 0.50, 0.75, 1.00, 1.25, 1.50),
    tail_top_fraction = 0.01,
    structural_crossing_budget = NA_integer_,
    run_wald_refinement = FALSE,
    run_derivative_assisted_inference = FALSE,
    retain_raw_data_in_outputs = FALSE
  ),
  compute = list(
    r_version = "4.4.2"
  )
)
