resolve_hillstrom_root <- function() {
  configured <- Sys.getenv("HILL_ROOT", unset = "")
  if (nzchar(configured)) {
    return(normalizePath(configured, winslash = "/", mustWork = TRUE))
  }
  file_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_argument)) {
    script <- normalizePath(sub("^--file=", "", file_argument[1L]),
                            winslash = "/", mustWork = TRUE)
    candidate <- normalizePath(file.path(dirname(script), ".."),
                               winslash = "/", mustWork = TRUE)
    if (file.exists(file.path(candidate, "VERSION")) &&
        dir.exists(file.path(candidate, "config"))) return(candidate)
  }
  candidate <- normalizePath(".", winslash = "/", mustWork = TRUE)
  if (file.exists(file.path(candidate, "VERSION")) &&
      dir.exists(file.path(candidate, "config"))) return(candidate)
  stop("Set HILL_ROOT to the hillstrom_application_v001 directory.", call. = FALSE)
}

resolve_cmo_root <- function(hill_root = resolve_hillstrom_root()) {
  configured <- Sys.getenv("CMO_ROOT", unset = "")
  candidates <- c(
    configured,
    file.path(hill_root, "..", "..", "reproducibility", "frozen_source")
  )
  candidates <- unique(candidates[nzchar(candidates)])
  for (candidate in candidates) {
    expanded <- tryCatch(
      normalizePath(candidate, winslash = "/", mustWork = TRUE),
      error = function(e) ""
    )
    if (nzchar(expanded) &&
        file.exists(file.path(expanded, "R", "common.R")) &&
        file.exists(file.path(expanded, "R", "method.R"))) return(expanded)
  }
  stop("Set CMO_ROOT to the frozen 2.0.0-bridge.8.1 engine directory.",
       call. = FALSE)
}

hill_root <- resolve_hillstrom_root()
cmo_engine_root <- resolve_cmo_root(hill_root)
Sys.setenv(CMO_ROOT = cmo_engine_root)
locked_common_r_sha256 <-
  "98f66209fa0e90834741261fe7d9597fe534335a055463e029e154381a074cf5"
common_r_path <- file.path(cmo_engine_root, "R", "common.R")
common_hash_output <- system2(
  "sha256sum", shQuote(common_r_path), stdout = TRUE, stderr = TRUE
)
common_hash_status <- attr(common_hash_output, "status")
if ((!is.null(common_hash_status) && as.integer(common_hash_status) != 0L) ||
    !length(common_hash_output)) {
  stop("Could not hash the frozen engine's R/common.R before sourcing it.",
       call. = FALSE)
}
common_hash_value <- strsplit(
  common_hash_output[1L], "[[:space:]]+"
)[[1L]][1L]
if (!identical(common_hash_value, locked_common_r_sha256)) {
  stop("R/common.R failed the pre-source SHA-256 gate.", call. = FALSE)
}
source(common_r_path, local = .GlobalEnv, chdir = FALSE)
rm(common_hash_output, common_hash_status, common_hash_value, common_r_path,
   locked_common_r_sha256)

sha256_file <- function(path) {
  path <- normalizePath(path, winslash = "/", mustWork = TRUE)
  output <- system2("sha256sum", shQuote(path), stdout = TRUE, stderr = TRUE)
  status <- attr(output, "status") %||% 0L
  assert_true(identical(as.integer(status), 0L),
              sprintf("sha256sum failed for %s", path))
  value <- strsplit(output[1L], "[[:space:]]+")[[1L]][1L]
  assert_true(grepl("^[0-9a-f]{64}$", value),
              sprintf("Invalid SHA-256 output for %s", path))
  value
}

with_working_directory <- function(path, expression) {
  previous <- getwd()
  on.exit(setwd(previous), add = TRUE)
  setwd(path)
  force(expression)
}

verify_engine_gate <- function(config) {
  assert_true(identical(trimws(readLines(file.path(cmo_engine_root, "VERSION"),
                                           warn = FALSE)[1L]),
                        config$engine$version),
              "The frozen engine version does not match the application lock.")
  manifest <- file.path(cmo_engine_root, "SOURCE_SHA256.txt")
  assert_true(file.exists(manifest),
              "The public engine source manifest is missing.")
  assert_true(identical(sha256_file(manifest), config$engine$manifest_sha256),
              "The public engine source manifest hash does not match the lock.")
  assert_true(identical(sha256_file(file.path(cmo_engine_root, "R", "common.R")),
                        config$engine$common_r_sha256),
              "R/common.R does not match the application lock.")
  assert_true(identical(sha256_file(file.path(cmo_engine_root, "R", "method.R")),
                        config$engine$method_r_sha256),
              "R/method.R does not match the application lock.")
  status <- with_working_directory(
    cmo_engine_root,
    system2("sha256sum", c("-c", "--status", "SOURCE_SHA256.txt"))
  )
  assert_true(identical(as.integer(status), 0L),
              "The public engine source verification failed.")
  invisible(TRUE)
}

source_hillstrom_engine <- function(config) {
  verify_engine_gate(config)
  source(file.path(cmo_engine_root, "R", "method.R"), local = .GlobalEnv,
         chdir = FALSE)
  invisible(TRUE)
}

load_hillstrom_config <- function(path = "config/application_v001.R") {
  candidate <- if (grepl("^/", path)) path else file.path(hill_root, path)
  full_path <- normalizePath(candidate, winslash = "/", mustWork = TRUE)
  canonical_path <- normalizePath(
    file.path(hill_root, "config", "application_v001.R"),
    winslash = "/", mustWork = TRUE
  )
  assert_true(identical(full_path, canonical_path),
              "Only the design-locked config/application_v001.R is permitted.")
  config <- dget(full_path)
  assert_true(is.list(config), "The application configuration must be a list.")
  assert_true(identical(config$application_id, "hillstrom_application_v001"),
              "Unexpected application_id.")
  assert_true(identical(config$application_version, "1.0.2"),
              "Unexpected application version.")
  interval <- config$analysis$order_intervals$K1
  assert_true(identical(as.numeric(interval), c(0.25, 1.50)),
              "The locked primary interval must be K1=[0.25,1.50].")
  assert_true(identical(config$analysis$moment_types, "absolute"),
              "Only the locked absolute-moment analysis is permitted.")
  assert_true(is.na(config$analysis$structural_crossing_budget),
              "No structural crossing budget may be imposed in this application.")
  config$config_path <- full_path
  config
}

normalise_column_names <- function(names_in) {
  names_in <- trimws(tolower(names_in))
  names_in <- gsub("[^a-z0-9]+", "_", names_in)
  gsub("^_+|_+$", "", names_in)
}

read_and_validate_hillstrom <- function(data_path, config, strict = TRUE) {
  data_path <- normalizePath(data_path, winslash = "/", mustWork = TRUE)
  if (isTRUE(strict)) {
    assert_true(identical(sha256_file(data_path), config$data$sha256),
                "The raw-data SHA-256 does not match the locked source file.")
    assert_true(identical(as.numeric(file.info(data_path)$size),
                          as.numeric(config$data$bytes)),
                "The raw-data byte count does not match the locked source file.")
  }
  data <- read.csv(data_path, stringsAsFactors = FALSE, check.names = FALSE,
                   na.strings = c("NA"), strip.white = FALSE)
  names(data) <- normalise_column_names(names(data))
  assert_true(!anyDuplicated(names(data)), "Column names are duplicated after normalization.")
  assert_true(identical(names(data), config$data$expected_columns),
              "The normalized CSV schema is not the locked 12-column schema.")
  if (isTRUE(strict)) {
    assert_true(nrow(data) == config$data$rows,
                "The CSV does not contain exactly 64,000 rows.")
  }
  assert_true(!anyNA(data), "The CSV contains missing values.")
  character_columns <- c("history_segment", "zip_code", "channel", "segment")
  assert_true(!any(vapply(data[character_columns], function(value) {
    any(!nzchar(trimws(as.character(value))))
  }, logical(1L))), "A categorical field is empty.")
  numeric_columns <- setdiff(names(data), character_columns)
  for (column in numeric_columns) {
    assert_true(is.numeric(data[[column]]) && all(is.finite(data[[column]])),
                sprintf("Column %s must be finite numeric data.", column))
  }
  for (column in c("mens", "womens", "newbie", "visit", "conversion")) {
    assert_true(all(data[[column]] %in% c(0, 1)),
                sprintf("Column %s is not binary.", column))
  }
  assert_true(all(data$recency == floor(data$recency) &
                    data$recency >= config$data$recency_range[1L] &
                    data$recency <= config$data$recency_range[2L]),
              "recency must contain integers from 1 through 12.")
  assert_true(all(data$history >= 0), "history contains a negative value.")
  assert_true(all(data$spend >= 0), "spend contains a negative value.")
  assert_true(all(data$conversion == as.integer(data$spend > 0)),
              "conversion and positive spend are inconsistent.")
  assert_true(all(data$visit >= data$conversion),
              "A conversion occurs without a recorded visit.")
  for (column in names(config$data$category_levels)) {
    assert_true(setequal(unique(data[[column]]),
                         config$data$category_levels[[column]]),
                sprintf("Categorical levels in %s do not match the locked source.",
                        column))
  }
  segment_order <- names(config$data$segment_counts)
  segment_counts <- table(factor(data$segment, levels = segment_order))
  names(segment_counts) <- segment_order
  positive_counts <- vapply(segment_order, function(segment) {
    sum(data$segment == segment & data$spend > 0)
  }, integer(1L))
  if (isTRUE(strict)) {
    assert_true(identical(as.integer(segment_counts),
                          as.integer(config$data$segment_counts)),
                "Treatment-arm counts do not match the locked data contract.")
    assert_true(identical(as.integer(positive_counts),
                          as.integer(config$data$positive_spend_counts)),
                "Positive-spend counts do not match the locked data contract.")
    assert_true(identical(range(data$spend), config$data$observed_spend_range),
                "The observed spend range does not match the locked data contract.")
  }
  x_rows <- data$segment == config$analysis$x_segment
  y_rows <- data$segment == config$analysis$y_segment
  assert_true(!any(x_rows & y_rows), "The two analysis samples overlap.")
  assert_true(any(x_rows) && any(y_rows), "One of the locked treatment arms is absent.")
  x <- as.numeric(data[[config$analysis$outcome]][x_rows])
  y <- as.numeric(data[[config$analysis$outcome]][y_rows])
  assert_true(all(is.finite(x)) && all(is.finite(y)) &&
                all(x >= 0) && all(y >= 0),
              "The two analysis vectors must be finite and nonnegative.")
  validation <- data.frame(
    check = c(
      "sha256", "bytes", "rows", "columns", "missing_values",
      "segment_counts", "positive_spend_counts", "spend_range",
      "binary_columns", "recency_range", "categorical_levels",
      "history_nonnegative", "conversion_spend_consistency",
      "visit_conversion_consistency", "samples_nonoverlapping"
    ),
    status = "PASSED",
    detail = c(
      if (isTRUE(strict)) config$data$sha256 else "synthetic_fixture",
      file.info(data_path)$size, nrow(data), ncol(data), sum(is.na(data)),
      paste(names(segment_counts), as.integer(segment_counts),
            sep = "=", collapse = ";"),
      paste(names(positive_counts), as.integer(positive_counts),
            sep = "=", collapse = ";"),
      paste(range(data$spend), collapse = ";"), "mens;womens;newbie;visit;conversion",
      paste(range(data$recency), collapse = ";"),
      paste(names(config$data$category_levels), collapse = ";"),
      min(data$history),
      sum(data$conversion != as.integer(data$spend > 0)),
      sum(data$visit < data$conversion),
      sum(x_rows & y_rows)
    ),
    stringsAsFactors = FALSE
  )
  list(data = data, x = x, y = y, validation = validation,
       data_path = data_path)
}

merge_adjacent_cells <- function(cells, selector, label) {
  selected <- cells[selector, c("left", "right"), drop = FALSE]
  if (!nrow(selected)) {
    return(data.frame(sign = character(0L), left = numeric(0L),
                      right = numeric(0L), left_anchor = numeric(0L),
                      right_anchor = numeric(0L), stringsAsFactors = FALSE))
  }
  tolerance <- 128 * .Machine$double.eps * max(1, abs(unlist(selected)))
  groups <- cumsum(c(TRUE, head(selected$right, -1L) <
                           tail(selected$left, -1L) - tolerance))
  output <- lapply(split(selected, groups), function(block) {
    data.frame(
      sign = label, left = min(block$left), right = max(block$right),
      left_anchor = mean(c(block$left[1L], block$right[1L])),
      right_anchor = mean(c(tail(block$left, 1L), tail(block$right, 1L))),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, output)
}

brackets_from_sign_intervals <- function(sign_intervals) {
  empty <- data.frame(
    lower = numeric(0L), upper = numeric(0L),
    left_sign = character(0L), right_sign = character(0L),
    endpoint_sign_sequence = character(0L), width = numeric(0L),
    stringsAsFactors = FALSE
  )
  if (nrow(sign_intervals) < 2L) return(empty)
  sign_intervals <- sign_intervals[
    order(sign_intervals$left, sign_intervals$right), , drop = FALSE
  ]
  changes <- which(head(sign_intervals$sign, -1L) !=
                     tail(sign_intervals$sign, -1L))
  if (!length(changes)) return(empty)
  output <- lapply(changes, function(index) {
    left <- sign_intervals[index, , drop = FALSE]
    right <- sign_intervals[index + 1L, , drop = FALSE]
    orientation <- if (left$sign == "negative" && right$sign == "positive") {
      "negative_to_positive"
    } else {
      "positive_to_negative"
    }
    data.frame(
      lower = left$right_anchor, upper = right$left_anchor,
      left_sign = left$sign, right_sign = right$sign,
      endpoint_sign_sequence = orientation,
      width = right$left_anchor - left$right_anchor,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, output)
}

assert_bracket_invariants <- function(brackets, sign_intervals, interval,
                                      tolerance = 1e-12) {
  expected_columns <- c(
    "lower", "upper", "left_sign", "right_sign",
    "endpoint_sign_sequence", "width"
  )
  assert_true(identical(names(brackets), expected_columns),
              "The certified-bracket table has unexpected columns.")
  if (!nrow(brackets)) return(invisible(TRUE))
  assert_true(
    all(is.finite(brackets$lower)) && all(is.finite(brackets$upper)) &&
      all(is.finite(brackets$width)),
    "A certified bracket has a nonfinite endpoint or width."
  )
  assert_true(
    all(brackets$lower < brackets$upper) &&
      all(brackets$lower >= interval[1L] - tolerance) &&
      all(brackets$upper <= interval[2L] + tolerance),
    "A certified bracket is not a strict subinterval of K1."
  )
  assert_true(
    all(brackets$left_sign %in% c("positive", "negative")) &&
      all(brackets$right_sign %in% c("positive", "negative")) &&
      all(brackets$left_sign != brackets$right_sign),
    "A certified bracket does not have opposite endpoint signs."
  )
  expected_sequences <- paste(
    brackets$left_sign, brackets$right_sign, sep = "_to_"
  )
  assert_true(
    identical(as.character(brackets$endpoint_sign_sequence),
              expected_sequences),
    "A certified bracket has an inconsistent endpoint-sign sequence."
  )
  scale <- pmax(1, abs(brackets$width),
                abs(brackets$upper - brackets$lower))
  assert_true(
    all(abs(brackets$width - (brackets$upper - brackets$lower)) <=
          tolerance * scale),
    "A certified bracket has an inconsistent width."
  )
  endpoint_is_certified <- function(value, sign) {
    any(sign_intervals$sign == sign &
          sign_intervals$left <= value + tolerance &
          sign_intervals$right >= value - tolerance)
  }
  assert_true(
    all(vapply(seq_len(nrow(brackets)), function(index) {
      endpoint_is_certified(brackets$lower[index], brackets$left_sign[index]) &&
        endpoint_is_certified(brackets$upper[index], brackets$right_sign[index])
    }, logical(1L))),
    "A bracket endpoint is outside its reported whole-cell sign region."
  )
  assert_true(all(diff(brackets$lower) >= -tolerance),
              "Certified brackets are not ordered by their lower endpoints.")
  if (nrow(brackets) > 1L) {
    assert_true(
      all(head(brackets$upper, -1L) <= tail(brackets$lower, -1L) + tolerance),
      "The interiors of certified brackets overlap."
    )
  }
  invisible(TRUE)
}

format_hillstrom_number <- function(value, digits = 6L) {
  formatC(value, digits = digits, format = "fg", flag = "#")
}

escape_hillstrom_latex <- function(value) {
  gsub("_", "\\_", as.character(value), fixed = TRUE)
}

render_hillstrom_table <- function(descriptive, anchor_moments,
                                   application_summary, outer_components,
                                   certified_brackets) {
  assert_true(nrow(application_summary) == 1L,
              "The article table requires a one-row application summary.")
  line_break <- "\\\\"
  format_order <- function(value) formatC(value, digits = 2L, format = "f")
  segment_for_role <- function(role) {
    index <- which(as.character(descriptive$role) == role)
    assert_true(length(index) == 1L,
                sprintf("The article table requires exactly one %s arm.", role))
    escape_hillstrom_latex(descriptive$segment[index])
  }
  interval_text <- sprintf(
    "$[%s,%s]$",
    format_order(application_summary$p_min),
    format_order(application_summary$p_max)
  )
  conclusion_text <- switch(
    as.character(application_summary$conclusion),
    UNIFORM_Y_DOMINANCE_NO_ROOT = sprintf(
      "%s moments uniformly exceed %s moments on %s; no critical order certified.",
      segment_for_role("Y"), segment_for_role("X"), interval_text
    ),
    UNIFORM_X_DOMINANCE_NO_ROOT = sprintf(
      "%s moments uniformly exceed %s moments on %s; no critical order certified.",
      segment_for_role("X"), segment_for_role("Y"), interval_text
    ),
    REVERSAL_AT_LEAST_ONE_ROOT = sprintf(
      "At least one ranking reversal is certified on %s.", interval_text
    ),
    INCONCLUSIVE = sprintf(
      "The simultaneous procedure is inconclusive on %s.", interval_text
    ),
    escape_hillstrom_latex(application_summary$conclusion)
  )
  component_word <- if (application_summary$outer_component_count == 1L) {
    "component"
  } else {
    "components"
  }
  arm_rows <- vapply(seq_len(nrow(descriptive)), function(index) {
    sprintf(
      "%s & %s & %d & %d & %s & %s & %s & %s %s",
      escape_hillstrom_latex(descriptive$role[index]),
      escape_hillstrom_latex(descriptive$segment[index]),
      descriptive$n[index], descriptive$n_positive[index],
      format_hillstrom_number(descriptive$positive_rate[index]),
      format_hillstrom_number(descriptive$mean_spend[index]),
      format_hillstrom_number(descriptive$sd_spend[index]),
      format_hillstrom_number(descriptive$conditional_mean_positive[index]),
      line_break
    )
  }, character(1L))
  anchor_rows <- vapply(seq_len(nrow(anchor_moments)), function(index) {
    sprintf(
      paste0(
        "%s & %s & %s & %s & $[%s,%s]$ & ",
        "$[%s,%s]$ & %s %s"
      ),
      format_order(anchor_moments$order[index]),
      format_hillstrom_number(anchor_moments$moment_x[index]),
      format_hillstrom_number(anchor_moments$moment_y[index]),
      format_hillstrom_number(anchor_moments$delta_hat[index]),
      format_hillstrom_number(anchor_moments$implemented_lower[index]),
      format_hillstrom_number(anchor_moments$implemented_upper[index]),
      format_hillstrom_number(anchor_moments$endpoint_lower[index]),
      format_hillstrom_number(anchor_moments$endpoint_upper[index]),
      escape_hillstrom_latex(anchor_moments$anchor_sign_status[index]),
      line_break
    )
  }, character(1L))
  outer_rows <- if (nrow(outer_components)) {
    vapply(seq_len(nrow(outer_components)), function(index) {
      sprintf(
        "Outer component %d & $[%s,%s]$; width %s %s",
        outer_components$component[index],
        format_hillstrom_number(outer_components$left[index]),
        format_hillstrom_number(outer_components$right[index]),
        format_hillstrom_number(outer_components$width[index]), line_break
      )
    }, character(1L))
  } else {
    paste("Outer components & none", line_break)
  }
  bracket_rows <- if (nrow(certified_brackets)) {
    vapply(seq_len(nrow(certified_brackets)), function(index) {
      sprintf(
        paste0(
          "Certified bracket %d & $[%s,%s]$; endpoint signs %s, %s %s"
        ),
        index,
        format_hillstrom_number(certified_brackets$lower[index]),
        format_hillstrom_number(certified_brackets$upper[index]),
        certified_brackets$left_sign[index],
        certified_brackets$right_sign[index], line_break
      )
    }, character(1L))
  } else {
    paste("Certified brackets & none", line_break)
  }
  c(
    "\\begingroup",
    "\\small",
    "\\setlength{\\tabcolsep}{4pt}",
    "\\renewcommand{\\arraystretch}{1.08}",
    "\\textbf{Panel A. Arm-level descriptive statistics}\\par\\smallskip",
    paste0(
      "\\noindent\\begin{tabularx}{\\linewidth}{@{}l",
      ">{\\raggedright\\arraybackslash}Xrrrrrr@{}}"
    ),
    "\\toprule",
    paste0(
      "Role & Arm & $n$ & $n_{>0}$ & $n_{>0}/n$ & $\\overline S$ & ",
      "$\\operatorname{SD}(S)$ & $\\overline S\\mid S>0$ ", line_break
    ),
    "\\midrule",
    arm_rows,
    "\\bottomrule",
    "\\end{tabularx}",
    "\\medskip",
    paste0(
      "\\textbf{Panel B. Moment contrast and simultaneous endpoints}",
      "\\par\\smallskip"
    ),
    paste0(
      "\\noindent\\begin{tabularx}{\\linewidth}{@{}rrrr",
      ">{\\centering\\arraybackslash}X",
      ">{\\centering\\arraybackslash}Xl@{}}"
    ),
    "\\toprule",
    paste0(
      "$p$ & $\\widehat M_X$ & $\\widehat M_Y$ & $\\widehat\\Delta$ & ",
      "Nominal 95\\% band & Cell enclosure & Sign ", line_break
    ),
    "\\midrule",
    anchor_rows,
    "\\bottomrule",
    "\\end{tabularx}",
    "\\medskip",
    "\\textbf{Panel C. Continuum certification summary}\\par\\smallskip",
    paste0(
      "\\noindent\\begin{tabularx}{\\linewidth}{@{}",
      ">{\\raggedright\\arraybackslash}p{0.27\\linewidth}",
      ">{\\raggedright\\arraybackslash}X@{}}"
    ),
    "\\toprule",
    sprintf(
      "Conclusion & %s %s", conclusion_text, line_break
    ),
    sprintf(
      paste0(
        "Outer-set summary & %d %s; total length %s %s"
      ),
      application_summary$outer_component_count,
      component_word,
      format_hillstrom_number(application_summary$outer_total_length),
      line_break
    ),
    outer_rows,
    bracket_rows,
    "\\bottomrule",
    "\\end{tabularx}",
    "\\par\\smallskip",
    paste0(
      "{\\footnotesize\\emph{Note:} Endpoint signs label certified anchors; ",
      "they do not assert uniqueness or crossing direction.}\\par"
    ),
    "\\endgroup"
  )
}

scientific_output_files <- function() {
  c(
    "application_summary.csv", "data_validation.csv",
    "descriptive_statistics.csv", "moment_band.csv",
    "moment_band_audit.csv", "anchor_moments.csv",
    "table_01_hillstrom_application.csv",
    "enclosure_cells.csv", "outer_set_components.csv",
    "certified_sign_intervals.csv", "certified_brackets.csv",
    "tail_diagnostics.csv", "bootstrap_suprema.csv", "audit_checks.csv",
    "rng_state_hashes.csv"
  )
}

verify_design_lock <- function() {
  lock <- file.path(hill_root, "DESIGN_LOCK_SHA256.txt")
  assert_true(file.exists(lock), "DESIGN_LOCK_SHA256.txt is missing.")
  status <- with_working_directory(
    hill_root,
    system2("sha256sum", c("-c", "--status", "DESIGN_LOCK_SHA256.txt"))
  )
  assert_true(identical(as.integer(status), 0L),
              "The Hillstrom design lock verification failed.")
  invisible(TRUE)
}

verify_source_manifest <- function() {
  manifest <- file.path(hill_root, "SOURCE_MANIFEST_SHA256.txt")
  assert_true(file.exists(manifest), "SOURCE_MANIFEST_SHA256.txt is missing.")
  status <- with_working_directory(
    hill_root,
    system2("sha256sum", c("-c", "--status", "SOURCE_MANIFEST_SHA256.txt"))
  )
  assert_true(identical(as.integer(status), 0L),
              "The Hillstrom source manifest verification failed.")
  invisible(TRUE)
}

write_result_manifest <- function(directory) {
  files <- sort(list.files(directory, recursive = TRUE, full.names = TRUE))
  files <- files[file.info(files)$isdir %in% FALSE]
  excluded <- "MANIFEST_SHA256.txt"
  files <- files[!basename(files) %in% excluded]
  relative <- substring(files, nchar(normalizePath(directory, winslash = "/",
                                                   mustWork = TRUE)) + 2L)
  hashes <- vapply(files, sha256_file, character(1L))
  lines <- sprintf("%s  ./%s", hashes, relative)
  atomic_write_lines(lines, file.path(directory, "MANIFEST_SHA256.txt"),
                     overwrite = TRUE)
  invisible(lines)
}

verify_result_manifest <- function(directory) {
  status <- with_working_directory(
    directory,
    system2("sha256sum", c("-c", "--status", "MANIFEST_SHA256.txt"))
  )
  identical(as.integer(status), 0L)
}
