#!/usr/bin/env Rscript

arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 2L) {
  stop(
    paste(
      "Usage: Rscript scripts/materialize_enclosure_depth_asset.R",
      "ASSET PROJECT"
    ),
    call. = FALSE
  )
}

asset <- normalizePath(arguments[[1L]], mustWork = TRUE)
project <- normalizePath(arguments[[2L]], mustWork = TRUE)
replications <- readRDS(asset)

required <- c(
  "study_id", "unit_id", "cell_id", "cell_class", "rep_id",
  "pairing_key", "variant_id", "scenario_id", "n_x", "n_y",
  "seed_id", "bootstrap_seed_id", "moment_type", "interval_id"
)
key <- c("unit_id", "moment_type", "interval_id")
stopifnot(
  is.data.frame(replications),
  nrow(replications) == 3048L,
  all(required %in% names(replications)),
  all(replications$study_id == "enclosure_depth_v001"),
  length(unique(replications$cell_id)) == 12L,
  length(unique(replications$unit_id)) == 3048L,
  length(unique(replications$seed_id)) == 1524L,
  length(unique(replications$bootstrap_seed_id)) == 1524L,
  length(unique(replications$pairing_key)) == 6L,
  setequal(
    unique(as.character(replications$variant_id)),
    c("depth6_reference", "depth8_refined")
  ),
  !anyDuplicated(replications[, key, drop = FALSE])
)

pair_index <- interaction(
  as.character(replications$pairing_key),
  replications$rep_id,
  drop = TRUE,
  sep = "::"
)
pair_groups <- split(seq_len(nrow(replications)), pair_index)
valid_pairs <- vapply(pair_groups, function(index) {
  length(index) == 2L &&
    setequal(
      as.character(replications$variant_id[index]),
      c("depth6_reference", "depth8_refined")
    ) &&
    length(unique(replications$seed_id[index])) == 1L &&
    length(unique(replications$bootstrap_seed_id[index])) == 1L &&
    length(unique(replications$scenario_id[index])) == 1L &&
    length(unique(replications$n_x[index])) == 1L &&
    length(unique(replications$n_y[index])) == 1L &&
    length(unique(replications$cell_class[index])) == 1L &&
    length(unique(replications$moment_type[index])) == 1L &&
    length(unique(replications$interval_id[index])) == 1L
}, logical(1L))
stopifnot(length(pair_groups) == 1524L, all(valid_pairs))

raw_directory <- file.path(
  project, "results", "raw", "enclosure_depth_v001", "task_00001"
)
status_directory <- file.path(project, "status", "enclosure_depth_v001")
dir.create(raw_directory, recursive = TRUE, showWarnings = FALSE)
dir.create(status_directory, recursive = TRUE, showWarnings = FALSE)

saveRDS(
  list(results = replications, errors = data.frame()),
  file.path(raw_directory, "wave_00001.rds"),
  version = 3L,
  compress = FALSE
)
validation <- list(
  complete = TRUE,
  numerically_clean = TRUE,
  expected_units = 3048L,
  expected_tasks = 12L,
  independent_sample_streams = 1524L,
  independent_bootstrap_streams = 1524L
)
saveRDS(
  validation,
  file.path(status_directory, "validation.rds"),
  version = 3L
)

cat("ENCLOSURE-DEPTH ASSET MATERIALIZATION: PASSED\n")
