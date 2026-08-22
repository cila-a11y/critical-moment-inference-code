#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
EXPECTED_R="4.4.2"
EXPECTED_BYTES="36063100"
EXPECTED_SHA256="eabac798fb2da2a126203d965e5de41012af7aad3b74358ee70302b8dd35606b"
ASSET="${1:-$ROOT/main_final_v003_replication_results.rds}"
CONFIG="config/main_final_v003.R"
STUDY="main_final_v003"
OUTPUT="$ROOT/generated/reproduction"

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
  echo "Missing release asset: $ASSET" >&2
  exit 66
}
test "$(wc -c < "$ASSET" | tr -d '[:space:]')" = "$EXPECTED_BYTES" || {
  echo "The replication asset has the wrong byte count." >&2
  exit 65
}
test "$(sha256sum "$ASSET" | awk '{print $1}')" = "$EXPECTED_SHA256" || {
  echo "The replication asset failed its SHA-256 check." >&2
  exit 65
}

WORK="$(mktemp -d "${TMPDIR:-/tmp}/cmo-reproduction.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT
cp -R "$ROOT/simulation" "$WORK/simulation"

CMO_ROOT="$WORK/simulation" Rscript --vanilla \
  "$WORK/simulation/R/prepare_design.R" --config="$CONFIG"
Rscript --vanilla "$ROOT/scripts/materialize_replication_asset.R" \
  "$ASSET" "$WORK/simulation"
CMO_ROOT="$WORK/simulation" Rscript --vanilla \
  "$WORK/simulation/R/summarize_results.R" --config="$CONFIG"
CMO_ROOT="$WORK/simulation" Rscript --vanilla \
  "$WORK/simulation/R/build_article_outputs.R" --config="$CONFIG"

SUMMARY="$WORK/simulation/results/summary/$STUDY"
Rscript --vanilla "$ROOT/scripts/build_final_tables.R" \
  "$SUMMARY" "$SUMMARY/article" \
  "$WORK/simulation/config/main_final_v003.R"

STAGE="$WORK/stage"
ARTICLE="$SUMMARY/article"
mkdir -p "$STAGE"
for file in \
  table_01_design.csv \
  table_02_benchmark.csv \
  table_03_root_inference.csv \
  table_04_special_power.csv \
  table_05_tail_strata.csv \
  figure_01_tangency_contraction.png \
  figure_02_two_root_power.png
do
  cp "$ARTICLE/$file" "$STAGE/$file"
done

rm -rf "$OUTPUT"
mkdir -p "$(dirname "$OUTPUT")"
mv "$STAGE" "$OUTPUT"
python3 "$ROOT/scripts/verify.py" --compare-simulation "$OUTPUT"
echo "SIMULATION REAGGREGATION FROM REPLICATIONS: PASSED"
