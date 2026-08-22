#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DATA="${1:-$ROOT/data/hillstrom_2008_raw.csv}"
EXPECTED_R="4.4.2"
EXPECTED_BYTES="3964977"
EXPECTED_SHA256="0e5893329d8b93cefecc571777672028290ab69865718020c78c7284f291aece"
OUTPUT="$ROOT/generated/hillstrom"

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
test -f "$DATA" || {
  echo "Missing Hillstrom data. Run scripts/acquire_hillstrom.sh first." >&2
  exit 66
}
test "$(wc -c < "$DATA" | tr -d '[:space:]')" = "$EXPECTED_BYTES" || {
  echo "The Hillstrom file has the wrong byte count." >&2
  exit 65
}
test "$(sha256sum "$DATA" | awk '{print $1}')" = "$EXPECTED_SHA256" || {
  echo "The Hillstrom file failed its SHA-256 check." >&2
  exit 65
}

WORK="$(mktemp -d "${TMPDIR:-/tmp}/cmo-hillstrom.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT
cp -R "$ROOT/application" "$WORK/application"
mkdir -p "$WORK/application/status/hillstrom_application_v001"

PROVENANCE="$WORK/application/status/hillstrom_application_v001/data_provenance.csv"
Rscript --vanilla - "$PROVENANCE" <<'RSCRIPT'
arguments <- commandArgs(trailingOnly = TRUE)
write.csv(
  data.frame(
    dataset = "MineThatData E-Mail Analytics and Data Mining Challenge",
    source_mode = "locked_download_or_local_copy",
    source_used = "upstream file matching the locked SHA-256",
    source_page = paste0(
      "https://blog.minethatdata.com/2008/03/",
      "minethatdata-e-mail-analytics-and-data.html"
    ),
    accessed_utc = format(Sys.time(), tz = "UTC", format = "%Y-%m-%dT%H:%M:%SZ"),
    raw_filename = "hillstrom_2008_raw.csv",
    bytes = 3964977,
    sha256 = "0e5893329d8b93cefecc571777672028290ab69865718020c78c7284f291aece",
    license_status = "no formal upstream redistribution license located",
    redistributed = "FALSE",
    stringsAsFactors = FALSE
  ),
  arguments[[1]], row.names = FALSE, na = ""
)
RSCRIPT

HILL_ROOT="$WORK/application" CMO_ROOT="$ROOT/simulation" HILL_DATA="$DATA" \
  Rscript --vanilla "$WORK/application/R/run_application.R" \
  --config=config/application_v001.R --data="$DATA" --run-id=primary
HILL_ROOT="$WORK/application" CMO_ROOT="$ROOT/simulation" HILL_DATA="$DATA" \
  Rscript --vanilla "$WORK/application/R/verify_application.R" \
  --config=config/application_v001.R --data="$DATA" --run-id=primary

RUN="$WORK/application/results/hillstrom_application_v001/primary"
STAGE="$WORK/stage"
mkdir -p "$STAGE"
for file in \
  application_summary.csv \
  data_validation.csv \
  descriptive_statistics.csv \
  moment_band.csv \
  moment_band_audit.csv \
  anchor_moments.csv \
  enclosure_cells.csv \
  outer_set_components.csv \
  certified_sign_intervals.csv \
  certified_brackets.csv \
  tail_diagnostics.csv \
  bootstrap_suprema.csv \
  audit_checks.csv \
  rng_state_hashes.csv
do
  cp "$RUN/$file" "$STAGE/$file"
done
cp "$RUN/table_01_hillstrom_application.csv" \
  "$STAGE/table_06_hillstrom_application.csv"
cp "$RUN/figure_01_hillstrom_moment_contrast.pdf" \
  "$STAGE/figure_03_hillstrom_moment_contrast.pdf"

rm -rf "$OUTPUT"
mkdir -p "$(dirname "$OUTPUT")"
mv "$STAGE" "$OUTPUT"
python3 "$ROOT/scripts/verify.py" --compare-hillstrom "$OUTPUT"
echo "HILLSTROM APPLICATION REPRODUCTION: PASSED"
