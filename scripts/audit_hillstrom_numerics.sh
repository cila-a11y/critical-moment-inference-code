#!/usr/bin/env bash
set -euo pipefail
export LC_ALL=C

repo=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
data=${1:-"$repo/data/hillstrom_2008_raw.csv"}
output=${2:-"$repo/generated/verified_hillstrom_audit"}
expected_sha=0e5893329d8b93cefecc571777672028290ab69865718020c78c7284f291aece
expected_bytes=3964977

for command_name in python3 sha256sum cc; do
  command -v "$command_name" >/dev/null 2>&1 || {
    echo "Required command not found: $command_name" >&2
    exit 2
  }
done

if [[ ! -f "$data" ]]; then
  echo "Hillstrom data file not found: $data" >&2
  exit 2
fi

actual_bytes=$(wc -c < "$data" | tr -d '[:space:]')
actual_sha=$(sha256sum "$data" | awk '{print $1}')
[[ "$actual_bytes" == "$expected_bytes" ]] || {
  echo "Hillstrom byte count changed: $actual_bytes" >&2
  exit 2
}
[[ "$actual_sha" == "$expected_sha" ]] || {
  echo "Hillstrom SHA-256 changed: $actual_sha" >&2
  exit 2
}

work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT
mkdir -p "$output"

python3 - "$data" "$repo/results/hillstrom/bootstrap_suprema.csv" \
  "$work/counts.csv" <<'PY'
import csv
import sys
from collections import Counter
from decimal import Decimal
from pathlib import Path

data_path, suprema_path, counts_path = map(Path, sys.argv[1:])
groups = {"No E-Mail": 0, "Mens E-Mail": 1}
counts = {0: Counter(), 1: Counter()}
with data_path.open(newline="", encoding="utf-8-sig") as handle:
    for row in csv.DictReader(handle):
        group = groups.get(row["segment"])
        if group is not None:
            counts[group][row["spend"]] += 1

if sum(counts[0].values()) != 21306 or sum(counts[1].values()) != 21307:
    raise SystemExit("Unexpected Hillstrom arm counts")

with counts_path.open("w", newline="", encoding="utf-8") as handle:
    writer = csv.writer(handle, lineterminator="\n")
    for group in (0, 1):
        for value, count in sorted(
            counts[group].items(), key=lambda item: Decimal(item[0])
        ):
            writer.writerow((group, value, count))

with suprema_path.open(newline="", encoding="utf-8") as handle:
    values = [Decimal(row["supremum"]) for row in csv.DictReader(handle)]
if len(values) != 9999:
    raise SystemExit("Expected 9999 frozen multiplier suprema")
critical_value = sorted(values)[9499]
if str(critical_value) != "2.22823852982312":
    raise SystemExit(f"Unexpected multiplier critical value: {critical_value}")
PY

if [[ -n ${CMO_MPFR_INCLUDEDIR:-} || -n ${CMO_MPFR_LIBDIR:-} ]]; then
  [[ -n ${CMO_MPFR_INCLUDEDIR:-} && -n ${CMO_MPFR_LIBDIR:-} ]] || {
    echo "Set both CMO_MPFR_INCLUDEDIR and CMO_MPFR_LIBDIR" >&2
    exit 2
  }
  [[ -f "$CMO_MPFR_INCLUDEDIR/mpfr.h" && -f "$CMO_MPFR_INCLUDEDIR/gmp.h" ]] || {
    echo "MPFR/GMP development headers were not found in $CMO_MPFR_INCLUDEDIR" >&2
    exit 2
  }
  cc -O2 -std=c11 "$repo/scripts/audit_hillstrom_numerics.c" \
    -I"$CMO_MPFR_INCLUDEDIR" \
    -L"$CMO_MPFR_LIBDIR" -Wl,-rpath,"$CMO_MPFR_LIBDIR" \
    -lmpfr -lgmp -lm -o "$work/audit_hillstrom_numerics"
elif command -v pkg-config >/dev/null 2>&1 && pkg-config --exists mpfr; then
  # shellcheck disable=SC2046
  cc -O2 -std=c11 "$repo/scripts/audit_hillstrom_numerics.c" \
    $(pkg-config --cflags --libs mpfr) -lgmp -lm \
    -o "$work/audit_hillstrom_numerics"
elif printf '#include <mpfr.h>\n#include <gmp.h>\n' | cc -E - >/dev/null 2>&1; then
  cc -O2 -std=c11 "$repo/scripts/audit_hillstrom_numerics.c" \
    -lmpfr -lgmp -lm -o "$work/audit_hillstrom_numerics"
else
  echo "MPFR/GMP development headers and libraries were not found" >&2
  echo "Install libmpfr-dev and libgmp-dev, or set CMO_MPFR_INCLUDEDIR and CMO_MPFR_LIBDIR" >&2
  exit 2
fi
"$work/audit_hillstrom_numerics" "$work/counts.csv" > "$work/result.txt"

python3 - "$work/result.txt" "$output" "$expected_sha" <<'PY'
import csv
import sys
from decimal import Decimal
from pathlib import Path

result_path = Path(sys.argv[1])
output = Path(sys.argv[2])
source_sha = sys.argv[3]
values = {}
for line in result_path.read_text(encoding="utf-8").splitlines():
    key, value = line.split("=", 1)
    values[key] = value

required = {
    "n_x": "21306",
    "n_y": "21307",
    "precision_bits": "256",
    "cells": "4000",
    "accuracy_target_met": "TRUE",
    "positive_cells": "4000",
    "outer_cells": "0",
}
for key, expected in required.items():
    if values.get(key) != expected:
        raise SystemExit(f"Audit check failed for {key}: {values.get(key)}")

def checked_decimal(key):
    value = Decimal(values[key])
    if not value.is_finite():
        raise SystemExit(f"Audit produced a non-finite value for {key}")
    return value

guard_ratio = checked_decimal("certified_guard_ratio_lower")
variance_lower = checked_decimal("certified_variance_lower")
enclosure_excess = checked_decimal("enclosure_excess_upper")
band_lower = checked_decimal("minimum_whole_cell_lower_bound")

if guard_ratio <= Decimal("1e-12"):
    raise SystemExit("The certified relative-variance guard is not positive")
if variance_lower <= 0:
    raise SystemExit("The certified variance lower bound is not positive")
if enclosure_excess > Decimal("0.1"):
    raise SystemExit("The whole-cell enclosure exceeds the accuracy target")
if band_lower <= 0:
    raise SystemExit("At least one whole-cell lower bound is not positive")

summary_fields = [
    "audit_id", "source_sha256", "critical_value", "precision_bits", "cells",
    "n_x", "n_y", "certified_guard_ratio_lower",
    "certified_variance_lower", "accuracy_target",
    "enclosure_excess_upper", "accuracy_target_met",
    "minimum_whole_cell_lower_bound", "positive_cells", "outer_cells",
]
summary = {
    "audit_id": "hillstrom_verified_numerics_v001",
    "source_sha256": source_sha,
    "critical_value": "2.22823852982312",
    "precision_bits": values["precision_bits"],
    "cells": values["cells"],
    "n_x": values["n_x"],
    "n_y": values["n_y"],
    "certified_guard_ratio_lower": values["certified_guard_ratio_lower"],
    "certified_variance_lower": values["certified_variance_lower"],
    "accuracy_target": "0.1",
    "enclosure_excess_upper": values["enclosure_excess_upper"],
    "accuracy_target_met": values["accuracy_target_met"],
    "minimum_whole_cell_lower_bound": values["minimum_whole_cell_lower_bound"],
    "positive_cells": values["positive_cells"],
    "outer_cells": values["outer_cells"],
}
output.mkdir(parents=True, exist_ok=True)
with (output / "verified_numerical_audit_summary.csv").open(
    "w", newline="", encoding="utf-8"
) as handle:
    writer = csv.DictWriter(handle, fieldnames=summary_fields, lineterminator="\n")
    writer.writeheader()
    writer.writerow(summary)

with (output / "verified_numerical_audit_anchors.csv").open(
    "w", newline="", encoding="utf-8"
) as handle:
    writer = csv.DictWriter(
        handle, fieldnames=("order", "lower", "upper", "sign"),
        lineterminator="\n"
    )
    writer.writeheader()
    for index in range(1, 7):
        writer.writerow({
            "order": values[f"anchor_{index}_order"],
            "lower": values[f"anchor_{index}_lower"],
            "upper": values[f"anchor_{index}_upper"],
            "sign": "positive",
        })
PY

if [[ ${CMO_AUDIT_COMPARE:-1} == 1 ]]; then
  cmp "$output/verified_numerical_audit_summary.csv" \
    "$repo/results/hillstrom/verified_numerical_audit_summary.csv"
  cmp "$output/verified_numerical_audit_anchors.csv" \
    "$repo/results/hillstrom/verified_numerical_audit_anchors.csv"
fi

echo "HILLSTROM VERIFIED NUMERICAL AUDIT: PASSED"
