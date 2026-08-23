#!/usr/bin/env python3
"""Verify the frozen enclosure-depth audit using the Python standard library."""

from __future__ import annotations

import argparse
import csv
import hashlib
import math
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results" / "simulation" / "enclosure_depth_v001"
ASSET = (
    ROOT
    / "inputs"
    / "simulation"
    / "enclosure_depth_v001"
    / "replication_results.rds"
)
CONFIG = ROOT / "simulation" / "config" / "enclosure_depth_v001.R"

EXPECTED_CONFIG_SHA256 = (
    "1ee6977fd82e0852f60e25b9afdc1c094d9461c8e7fe2bed8179bbfa15bfcd7e"
)
EXPECTED_ASSET_SHA256 = (
    "f66e1ba64b37c8f76b2f30ebef554187e68c8cfcf555d8bbe04d4e2f7ee27ecd"
)
EXPECTED_ASSET_BYTES = 405_328

EXPECTED_RESULTS = {
    "summary.csv": (
        12,
        "056541351fed22c47dd88e275835bed8e5b6c0722a391642473d79383bbf98aa",
    ),
    "binary_metrics.csv": (
        1284,
        "b27aacfd9365a837b8c82413836bd175c676e2715e1d0e3a01362af702eb1892",
    ),
    "continuous_metrics.csv": (
        1464,
        "43f35db0f64bba274ffe6faa81a625644667b3f6e36a978cea044ae590c51fd0",
    ),
    "paired_comparisons.csv": (
        1374,
        "9756ac144d29d37d34f6dc746b01c06864c7aae737b29c8ef9a38085c2e4eb0a",
    ),
    "enclosure_depth_audit.csv": (
        6,
        "e8289985ee14f915cf4603356324d6ebb2b962cc22d19c34aa1769e78eee2db5",
    ),
    "across_interval_summary.csv": (
        12,
        "37e8678e2e2cb06b52a03671b8cdf35791277bbdea3562ebd3702e28947352d7",
    ),
    "across_interval_binary_metrics.csv": (
        48,
        "ea4362ed479453a39bc594dc89a34fb028959fac279df727715ad7d299d7a19a",
    ),
    "root_metrics_conditional_on_joint_isolation.csv": (
        12,
        "5288def60cd97837bba0d7ae22e9d66221c08a1e5ee5c67a30372299477553f8",
    ),
    "numerical_failure_summary_by_cell.csv": (
        12,
        "da7c7446e9114cd4df0d273aec4f6bf31e5aaa8c90005599a1d4d956347c3195",
    ),
    "numerical_failures.csv": (
        0,
        "bd85bcdb8d4e613a79cb62d0903946ad10c83e63dc75f67614c159c0dbf4d184",
    ),
    "tail_stratified.csv": (
        12,
        "e7ff304c8fa8a09db3b8ae117110c43d5fa209864b93a3c31da94dcc1e3a8c34",
    ),
}
EXPECTED_SUMMARY_MANIFEST_SHA256 = (
    "e8bfabe97fb6342df63ff837d5e5ed154212bfa55b714fdc68d1f68106745c8c"
)

# These bridge.7 computational files are byte-identical to the public engine.
# The later summarizer adds reporting fields and is checked by reaggregation.
EXPECTED_CORE_SOURCE_SHA256 = {
    "augmented.R": "4e44564217ff0f84ce477e9d3808ec11f3035bb44fec064a899d1a95ef01c52d",
    "common.R": "98f66209fa0e90834741261fe7d9597fe534335a055463e029e154381a074cf5",
    "dgp.R": "244d2b2874c779b4da67b2a5e4e9408943a5503f9f1165c8d5f795ebcf11f4f9",
    "method.R": "d12fae4433ba57efb3d1d49a86ea2abcdd35f971b01920bb8d057b41f0efde63",
    "prepare_design.R": "71e3379cc11779b0d2872eabb8cdc0bbaee0479f805d990c64f720694da8900f",
    "reporting.R": "202149691dcc9b448d6c93faed204fddef37e7f2bf4e2e16f803a1fa1397508f",
    "run_task.R": "08143087bbf90cb7a591678d4213a4c0571651ec1452085e866ac7791dbebfd4",
    "validate_results.R": "2b626f1ca315c4a0a48d5e12fa98d5cb40d97b699e4693d49248b3b93a40631b",
}


class Audit:
    def __init__(self) -> None:
        self.checks = 0
        self.errors: list[str] = []

    def require(self, condition: bool, message: str) -> None:
        self.checks += 1
        if not condition:
            self.errors.append(message)


def digest(path: Path) -> str:
    hasher = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            hasher.update(block)
    return hasher.hexdigest()


def read_csv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        return list(reader.fieldnames or []), list(reader)


def integer(row: dict[str, str], field: str) -> int:
    return int(row[field])


def number(row: dict[str, str], field: str) -> float:
    return float(row[field])


def verify_frozen_files(audit: Audit) -> None:
    audit.require(CONFIG.is_file(), "The enclosure-depth configuration is missing")
    if CONFIG.is_file():
        audit.require(
            digest(CONFIG) == EXPECTED_CONFIG_SHA256,
            "The enclosure-depth configuration changed",
        )

    audit.require(ASSET.is_file(), "The replication-level depth asset is missing")
    if ASSET.is_file():
        audit.require(
            ASSET.stat().st_size == EXPECTED_ASSET_BYTES,
            "The replication-level depth asset has the wrong byte count",
        )
        audit.require(
            digest(ASSET) == EXPECTED_ASSET_SHA256,
            "The replication-level depth asset changed",
        )

    for name, (expected_rows, expected_hash) in EXPECTED_RESULTS.items():
        path = RESULTS / name
        audit.require(path.is_file(), f"Missing frozen depth result: {name}")
        if not path.is_file():
            continue
        audit.require(digest(path) == expected_hash, f"Frozen result changed: {name}")
        _, rows = read_csv(path)
        audit.require(len(rows) == expected_rows, f"Wrong frozen row count: {name}")

    manifest = RESULTS / "summary_manifest.txt"
    audit.require(manifest.is_file(), "The depth summary manifest is missing")
    if manifest.is_file():
        audit.require(
            digest(manifest) == EXPECTED_SUMMARY_MANIFEST_SHA256,
            "The depth summary manifest changed",
        )

    source_directory = ROOT / "simulation" / "R"
    for name, expected_hash in EXPECTED_CORE_SOURCE_SHA256.items():
        path = source_directory / name
        audit.require(path.is_file(), f"Missing core simulation source: {name}")
        if path.is_file():
            audit.require(
                digest(path) == expected_hash,
                f"Core simulation source no longer matches bridge.7: {name}",
            )


def verify_scientific_contract(audit: Audit) -> None:
    _, summary = read_csv(RESULTS / "summary.csv")
    audit.require(len(summary) == 12, "The depth summary must contain 12 cells")
    audit.require(
        all(row["study_id"] == "enclosure_depth_v001" for row in summary),
        "The depth study identifier changed",
    )
    audit.require(
        {row["variant_id"] for row in summary}
        == {"depth6_reference", "depth8_refined"},
        "The depth variants changed",
    )
    audit.require(
        all(row["expected_replications"] == "254" for row in summary),
        "The expected replication count per cell changed",
    )
    audit.require(
        all(row["successful_replications"] == "254" for row in summary),
        "A depth cell is incomplete",
    )
    audit.require(
        all(row["replication_failure_rate"] == "0" for row in summary),
        "A depth cell has replication failures",
    )
    audit.require(
        all(row["bootstrap_reps"] == "999" for row in summary),
        "The bootstrap replication count changed",
    )

    _, depth_rows = read_csv(RESULTS / "enclosure_depth_audit.csv")
    audit.require(len(depth_rows) == 6, "The paired audit must contain six cells")
    audit.require(
        all(row["reference_variant_id"] == "depth6_reference" for row in depth_rows),
        "The paired reference variant changed",
    )
    audit.require(
        all(row["comparison_variant_id"] == "depth8_refined" for row in depth_rows),
        "The paired comparison variant changed",
    )
    audit.require(
        all(row["reference_level"] == "6" and row["comparison_level"] == "8"
            for row in depth_rows),
        "The paired enclosure levels changed",
    )
    audit.require(
        sum(integer(row, "paired_replications") for row in depth_rows) == 1524,
        "The paired-replication total changed",
    )
    audit.require(
        all(integer(row, "subset_checks") == 254 for row in depth_rows),
        "The subset check denominator changed",
    )
    for field in (
        "subset_failures",
        "root_losses_under_continuum_coverage",
        "length_bound_failures",
        "decision_disagreements",
        "decision_disagreement_replications",
    ):
        audit.require(
            sum(integer(row, field) for row in depth_rows) == 0,
            f"The paired audit has nonzero {field}",
        )
    audit.require(
        sum(integer(row, "root_retention_checks") for row in depth_rows) == 1454,
        "The continuum-covered root-retention denominator changed",
    )
    audit.require(
        sum(integer(row, "length_bound_checks") for row in depth_rows) == 1524,
        "The length-bound check denominator changed",
    )
    audit.require(
        max(number(row, "maximum_length_bound_excess") for row in depth_rows)
        <= 1e-10,
        "A numerical length bound was exceeded",
    )
    audit.require(
        min(number(row, "minimum_length_reduction") for row in depth_rows)
        >= -1e-10,
        "The deeper outer set increased length beyond tolerance",
    )

    _, failures = read_csv(RESULTS / "numerical_failure_summary_by_cell.csv")
    audit.require(len(failures) == 12, "The numerical audit must contain 12 cells")
    audit.require(
        sum(integer(row, "expected_units") for row in failures) == 3048,
        "The numerical-audit unit total changed",
    )
    audit.require(
        all(row["successful_units"] == "254" for row in failures),
        "A numerical-audit cell is incomplete",
    )
    audit.require(
        sum(integer(row, "numerical_failure_units") for row in failures) == 0,
        "The frozen depth study has numerical failures",
    )

    _, binary = read_csv(RESULTS / "binary_metrics.csv")
    for metric in ("enclosure_node_limit_hit", "enclosure_variance_limit_hit"):
        selected = [row for row in binary if row["metric"] == metric]
        audit.require(len(selected) == 12, f"Missing binary metric: {metric}")
        audit.require(
            sum(integer(row, "denominator") for row in selected) == 3048,
            f"Wrong denominator for {metric}",
        )
        audit.require(
            sum(integer(row, "successes") for row in selected) == 0,
            f"The frozen depth study has nonzero {metric}",
        )
    depth_limits = [
        row for row in binary if row["metric"] == "enclosure_depth_limit_hit"
    ]
    audit.require(
        sum(integer(row, "successes") for row in depth_limits) == 1010,
        "The recorded depth-limit incidence changed",
    )

    _, paired = read_csv(RESULTS / "paired_comparisons.csv")
    audit.require(
        all(row["reference_variant_id"] == "depth6_reference" for row in paired),
        "A paired comparison has the wrong reference",
    )
    audit.require(
        all(row["comparison_variant_id"] == "depth8_refined" for row in paired),
        "A paired comparison has the wrong direction",
    )
    selected_metrics: dict[str, list[dict[str, str]]] = {}
    for metric in ("elapsed_seconds", "enclosure_nodes_used", "peak_r_heap_mb"):
        selected_metrics[metric] = [row for row in paired if row["metric"] == metric]
        audit.require(
            len(selected_metrics[metric]) == 6,
            f"The paired cost metric {metric} must have six rows",
        )
        audit.require(
            all(row["metric_pair_denominator"] == "254"
                for row in selected_metrics[metric]),
            f"The paired cost denominator changed for {metric}",
        )

    elapsed = [number(row, "mean_difference") for row in selected_metrics["elapsed_seconds"]]
    nodes = [number(row, "mean_difference") for row in selected_metrics["enclosure_nodes_used"]]
    audit.require(
        min(elapsed) >= 0.012 and max(elapsed) <= 0.058,
        "The elapsed-time differences no longer lie between 0.012 and 0.058 seconds",
    )
    audit.require(
        min(nodes) >= 8.0 and max(nodes) <= 20.0,
        "The enclosure-node differences no longer lie between 8 and 20 nodes",
    )
    heap_relative_changes = [
        abs(number(row, "mean_difference")) / number(row, "reference_mean")
        for row in selected_metrics["peak_r_heap_mb"]
    ]
    audit.require(
        max(heap_relative_changes) < 0.005,
        "The recorded R-heap change is no longer below 0.5 percent",
    )

    manifest_values: dict[str, str] = {}
    for line in (RESULTS / "summary_manifest.txt").read_text(encoding="utf-8").splitlines():
        if "=" in line:
            key, value = line.split("=", 1)
            manifest_values[key] = value
    expected_manifest_values = {
        "study_id": "enclosure_depth_v001",
        "schema_version": "2.0.0",
        "replication_rows": "3048",
        "paired_comparison_rows": "1374",
        "enclosure_depth_audit_rows": "6",
        "expected_units": "3048",
        "numerical_failures": "0",
        "numerically_clean": "TRUE",
        "paired_difference_direction": (
            "comparison_variant_minus_reference_variant"
        ),
    }
    for key, value in expected_manifest_values.items():
        audit.require(
            manifest_values.get(key) == value,
            f"The summary-manifest value changed: {key}",
        )


def values_agree(expected: str, actual: str) -> bool:
    if expected == actual:
        return True
    try:
        left = float(expected)
        right = float(actual)
    except ValueError:
        return False
    if math.isnan(left) or math.isnan(right):
        return math.isnan(left) and math.isnan(right)
    return math.isclose(left, right, rel_tol=1e-12, abs_tol=1e-12)


def compare_reaggregation(audit: Audit, generated: Path) -> None:
    audit.require(generated.is_dir(), "The generated depth directory is missing")
    if not generated.is_dir():
        return
    for name in EXPECTED_RESULTS:
        expected_path = RESULTS / name
        actual_path = generated / name
        audit.require(actual_path.is_file(), f"Generated depth result is missing: {name}")
        if not actual_path.is_file():
            continue
        expected_fields, expected_rows = read_csv(expected_path)
        actual_fields, actual_rows = read_csv(actual_path)
        missing_fields = [field for field in expected_fields if field not in actual_fields]
        audit.require(
            not missing_fields,
            f"Generated {name} lacks frozen fields: {', '.join(missing_fields)}",
        )
        audit.require(
            len(actual_rows) == len(expected_rows),
            f"Generated row count differs: {name}",
        )
        if missing_fields or len(actual_rows) != len(expected_rows):
            continue
        for row_number, (expected_row, actual_row) in enumerate(
            zip(expected_rows, actual_rows), start=1
        ):
            for field in expected_fields:
                audit.require(
                    values_agree(expected_row[field], actual_row[field]),
                    f"Generated value differs in {name}, row {row_number}, field {field}",
                )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--compare", type=Path)
    args = parser.parse_args()

    audit = Audit()
    verify_frozen_files(audit)
    verify_scientific_contract(audit)
    if args.compare:
        compare_reaggregation(audit, args.compare.resolve())

    if audit.errors:
        print("ENCLOSURE-DEPTH AUDIT: FAILED", file=sys.stderr)
        for error in audit.errors:
            print(f"- {error}", file=sys.stderr)
        return 1
    print(f"ENCLOSURE-DEPTH AUDIT: PASSED ({audit.checks} checks)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
