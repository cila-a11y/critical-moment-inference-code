#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
EXPECTED_R="4.4.2"
EXPECTED_BYTES="405328"
EXPECTED_SHA256="f66e1ba64b37c8f76b2f30ebef554187e68c8cfcf555d8bbe04d4e2f7ee27ecd"
ASSET="${1:-$ROOT/inputs/simulation/enclosure_depth_v001/replication_results.rds}"
CONFIG="config/enclosure_depth_v001.R"
STUDY="enclosure_depth_v001"
OUTPUT="$ROOT/generated/enclosure-depth-reproduction"

command -v Rscript >/dev/null 2>&1 || {
  echo "Rscript is required." >&2
  exit 69
}
command -v sha256sum >/dev/null 2>&1 || {
  echo "GNU sha256sum is required." >&2
  exit 69
}
ACTUAL_R="$(Rscript --vanilla -e 'cat(paste(R.version$major, R.version$minor, sep="."))')"
test "$ACTUAL_R" = "$EXPECTED_R" || {
  echo "R $EXPECTED_R is required; found R $ACTUAL_R." >&2
  exit 69
}
test -f "$ASSET" || {
  echo "Missing replication asset: $ASSET" >&2
  exit 66
}
test "$(wc -c < "$ASSET" | tr -d '[:space:]')" = "$EXPECTED_BYTES" || {
  echo "The enclosure-depth asset has the wrong byte count." >&2
  exit 65
}
test "$(sha256sum "$ASSET" | awk '{print $1}')" = "$EXPECTED_SHA256" || {
  echo "The enclosure-depth asset failed its SHA-256 check." >&2
  exit 65
}

WORK="$(mktemp -d "${TMPDIR:-/tmp}/cmo-depth-reproduction.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT
cp -R "$ROOT/simulation" "$WORK/simulation"
(
  cd "$WORK/simulation"
  sha256sum -c --status SOURCE_SHA256.txt
) || {
  echo "The copied simulation source failed its SHA-256 check." >&2
  exit 65
}

CMO_ROOT="$WORK/simulation" Rscript --vanilla \
  "$WORK/simulation/R/prepare_design.R" --config="$CONFIG"
Rscript --vanilla "$ROOT/scripts/materialize_enclosure_depth_asset.R" \
  "$ASSET" "$WORK/simulation"
CMO_ROOT="$WORK/simulation" Rscript --vanilla \
  "$WORK/simulation/R/summarize_results.R" --config="$CONFIG"

SUMMARY="$WORK/simulation/results/summary/$STUDY"
STAGE="$WORK/stage"
mkdir -p "$STAGE"
for file in \
  summary.csv \
  binary_metrics.csv \
  continuous_metrics.csv \
  across_interval_summary.csv \
  across_interval_binary_metrics.csv \
  root_metrics_conditional_on_joint_isolation.csv \
  tail_stratified.csv \
  paired_comparisons.csv \
  enclosure_depth_audit.csv \
  numerical_failures.csv \
  numerical_failure_summary_by_cell.csv
do
  cp "$SUMMARY/$file" "$STAGE/$file"
done

python3 "$ROOT/scripts/verify_enclosure_depth.py" --compare "$STAGE"
rm -rf "$OUTPUT"
mkdir -p "$(dirname "$OUTPUT")"
mv "$STAGE" "$OUTPUT"
echo "ENCLOSURE-DEPTH REAGGREGATION: PASSED"
