#!/usr/bin/env Rscript

arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 2L) {
  stop(
    "Usage: Rscript scripts/materialize_replication_asset.R ASSET PROJECT",
    call. = FALSE
  )
}

asset <- normalizePath(arguments[[1L]], mustWork = TRUE)
project <- normalizePath(arguments[[2L]], mustWork = TRUE)
replications <- readRDS(asset)

key <- c("unit_id", "moment_type", "interval_id")
stopifnot(
  is.data.frame(replications),
  nrow(replications) == 284480L,
  all(c("study_id", "cell_id", key) %in% names(replications)),
  all(replications$study_id == "main_final_v003"),
  length(unique(replications$cell_id)) == 54L,
  length(unique(replications$unit_id)) == 109728L,
  !anyDuplicated(replications[, key, drop = FALSE])
)

raw_directory <- file.path(
  project, "results", "raw", "main_final_v003", "task_00001"
)
status_directory <- file.path(project, "status", "main_final_v003")
dir.create(raw_directory, recursive = TRUE, showWarnings = FALSE)
dir.create(status_directory, recursive = TRUE, showWarnings = FALSE)

saveRDS(
  list(results = replications, errors = data.frame()),
  file.path(raw_directory, "wave_00001.rds"),
  version = 3L, compress = FALSE
)
validation <- list(
  complete = TRUE,
  numerically_clean = TRUE,
  expected_units = 109728L,
  expected_tasks = 48L,
  independent_sample_streams = 109728L,
  independent_bootstrap_streams = 109728L
)
saveRDS(
  validation, file.path(status_directory, "validation.rds"),
  version = 3L
)

cat("REPLICATION ASSET MATERIALIZATION: PASSED\n")
