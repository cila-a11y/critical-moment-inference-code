#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
EXPECTED_R="4.4.2"
CONFIG="config/main_final_v003.R"
STUDY="main_final_v003"
OUTPUT="$ROOT/generated/rebuild"

command -v Rscript >/dev/null 2>&1 || {
  echo "Rscript is required." >&2
  exit 69
}
ACTUAL_R="$(Rscript --vanilla -e 'cat(paste(R.version$major, R.version$minor, sep="."))')"
test "$ACTUAL_R" = "$EXPECTED_R" || {
  echo "R $EXPECTED_R is required; found R $ACTUAL_R." >&2
  exit 69
}

WORK="$(mktemp -d "${TMPDIR:-/tmp}/cmo-rebuild.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT
cp -R "$ROOT/simulation" "$WORK/simulation"
mkdir -p \
  "$WORK/simulation/results/summary/$STUDY" \
  "$WORK/simulation/status/$STUDY"

for file in \
  summary.csv \
  binary_metrics.csv \
  continuous_metrics.csv \
  across_interval_summary.csv \
  across_interval_binary_metrics.csv \
  root_collection_summary_deduplicated.csv \
  root_metrics_deduplicated_by_root.csv \
  numerical_failure_summary_by_cell.csv \
  tail_stratified.csv
do
  cp "$ROOT/inputs/simulation/$file" \
    "$WORK/simulation/results/summary/$STUDY/$file"
done

Rscript --vanilla -e \
  'saveRDS(list(complete=TRUE), commandArgs(TRUE)[1], version=3L)' \
  "$WORK/simulation/status/$STUDY/validation.rds"

CMO_ROOT="$WORK/simulation" Rscript --vanilla \
  "$WORK/simulation/R/build_article_outputs.R" --config="$CONFIG"
SUMMARY="$WORK/simulation/results/summary/$STUDY"
Rscript --vanilla "$ROOT/scripts/build_final_tables.R" \
  "$SUMMARY" "$SUMMARY/article" \
  "$WORK/simulation/config/main_final_v003.R"

STAGE="$WORK/stage"
mkdir -p "$STAGE"
ARTICLE="$WORK/simulation/results/summary/$STUDY/article"
for file in \
  table_01_design.csv \
  table_02_benchmark.csv \
  table_03_root_inference.csv \
  table_04_special_power.csv \
  table_05_tail_strata.csv \
  table_07_ablation.csv \
  figure_01_tangency_contraction.png \
  figure_02_two_root_power.png
do
  cp "$ARTICLE/$file" "$STAGE/$file"
done

rm -rf "$OUTPUT"
mkdir -p "$(dirname "$OUTPUT")"
mv "$STAGE" "$OUTPUT"
python3 "$ROOT/scripts/verify.py" --compare-simulation "$OUTPUT"
echo "SIMULATION OUTPUT REBUILD: PASSED"
