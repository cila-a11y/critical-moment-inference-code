#!/usr/bin/env python3
"""Verify the public reproducibility package using the Python standard library."""

from __future__ import annotations

import argparse
import csv
import hashlib
import re
import struct
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "SHA256SUMS.txt"

SIMULATION_FILES = (
    "table_01_design.csv",
    "table_02_benchmark.csv",
    "table_03_root_inference.csv",
    "table_04_special_power.csv",
    "table_05_tail_strata.csv",
)
SIMULATION_FIGURES = (
    "figure_01_tangency_contraction.png",
    "figure_02_two_root_power.png",
)
HILLSTROM_FILES = (
    "application_summary.csv",
    "data_validation.csv",
    "descriptive_statistics.csv",
    "moment_band.csv",
    "moment_band_audit.csv",
    "anchor_moments.csv",
    "enclosure_cells.csv",
    "outer_set_components.csv",
    "certified_sign_intervals.csv",
    "certified_brackets.csv",
    "tail_diagnostics.csv",
    "bootstrap_suprema.csv",
    "audit_checks.csv",
    "rng_state_hashes.csv",
    "table_06_hillstrom_application.csv",
)

BANNED_PATTERNS = (
    (re.compile(rb"/(?:projects|home)/", re.IGNORECASE), "absolute cluster path"),
    (re.compile(rb"[a-z]:\\\\users\\\\", re.IGNORECASE), "Windows user path"),
    (re.compile(rb"id_(?:rsa|ed25519)", re.IGNORECASE), "private-key filename"),
    (re.compile(rb"f[0-9]{9}hpc[a-z0-9]+", re.IGNORECASE), "compute account"),
    (re.compile(rb"[a-z0-9.-]+\\.macc\\.fccn\\.pt", re.IGNORECASE), "cluster host"),
    (
        re.compile(rb"(?:job|jobid|slurm_job_id)[=,: ]+[0-9]{6,}", re.IGNORECASE),
        "scheduler job identifier",
    ),
)
BANNED_SUFFIXES = (".tar.gz", ".zip", ".bundle", ".log", ".out", ".err")
IGNORED_TOP_LEVEL = {".git", "data", "generated"}
IGNORED_PARTS = {"__pycache__"}
IGNORED_ROOT_FILES = {"main_final_v003_replication_results.rds"}


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


def csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def manifest_entries(audit: Audit) -> dict[str, str]:
    audit.require(MANIFEST.is_file(), "SHA256SUMS.txt is missing")
    if not MANIFEST.is_file():
        return {}
    entries: dict[str, str] = {}
    pattern = re.compile(r"^([0-9a-f]{64})  (.+)$")
    for line in MANIFEST.read_text(encoding="utf-8").splitlines():
        match = pattern.fullmatch(line)
        audit.require(match is not None, f"Malformed manifest line: {line}")
        if match:
            relative = match.group(2)
            audit.require(
                not relative.startswith("/") and ".." not in Path(relative).parts,
                f"Unsafe manifest path: {relative}",
            )
            audit.require(relative not in entries, f"Duplicate manifest path: {relative}")
            entries[relative] = match.group(1)
    return entries


def public_files() -> set[str]:
    files: set[str] = set()
    for path in ROOT.rglob("*"):
        if not path.is_file():
            continue
        relative = path.relative_to(ROOT)
        if relative.parts[0] in IGNORED_TOP_LEVEL:
            continue
        if len(relative.parts) == 1 and relative.name in IGNORED_ROOT_FILES:
            continue
        if any(part in IGNORED_PARTS for part in relative.parts) or path.suffix == ".pyc":
            continue
        if relative.as_posix() == "SHA256SUMS.txt":
            continue
        files.add(relative.as_posix())
    return files


def verify_manifest_and_privacy(audit: Audit) -> None:
    entries = manifest_entries(audit)
    audit.require(set(entries) == public_files(), "Manifest file list changed")
    for relative, expected in entries.items():
        path = ROOT / relative
        audit.require(path.is_file(), f"Missing file: {relative}")
        if not path.is_file():
            continue
        audit.require(digest(path) == expected, f"SHA-256 changed: {relative}")
        lower_name = relative.lower()
        audit.require(
            not lower_name.endswith(BANNED_SUFFIXES),
            f"Forbidden archive/log file: {relative}",
        )
        content = path.read_bytes()
        for pattern, label in BANNED_PATTERNS:
            audit.require(
                pattern.search(content) is None,
                f"Forbidden {label} in {relative}",
            )


def verify_simulation(audit: Audit) -> None:
    input_directory = ROOT / "inputs" / "simulation"
    expected_rows = {
        "summary.csv": 140,
        "binary_metrics.csv": 14980,
        "continuous_metrics.csv": 17080,
        "across_interval_summary.csv": 56,
        "across_interval_binary_metrics.csv": 224,
        "root_collection_summary_deduplicated.csv": 33,
        "root_metrics_deduplicated_by_root.csv": 40,
        "numerical_failure_summary_by_cell.csv": 54,
        "tail_stratified.csv": 227,
    }
    for name, count in expected_rows.items():
        path = input_directory / name
        audit.require(path.is_file(), f"Missing simulation input: {name}")
        if path.is_file():
            audit.require(len(csv_rows(path)) == count, f"Wrong row count: {name}")

    failures = csv_rows(input_directory / "numerical_failure_summary_by_cell.csv")
    audit.require(
        all(row["numerical_failure_units"] == "0" for row in failures),
        "Simulation numerical failures are not zero",
    )
    audit.require(
        sum(int(row["expected_units"]) for row in failures) == 109728,
        "Simulation unit total changed",
    )
    audit.require(
        all(row["expected_units"] == "2032" for row in failures),
        "Simulation replications per cell changed",
    )
    summary = csv_rows(input_directory / "summary.csv")
    audit.require(
        all(row["expected_replications"] == "2032" for row in summary),
        "Simulation summary replication count changed",
    )
    audit.require(
        all(row["bootstrap_reps"] == "999" for row in summary),
        "Simulation bootstrap count changed",
    )

    output_directory = ROOT / "results" / "simulation"
    for name in SIMULATION_FILES + SIMULATION_FIGURES:
        audit.require((output_directory / name).is_file(), f"Missing result: {name}")
    expected_output_rows = {
        "table_01_design.csv": 8,
        "table_02_benchmark.csv": 10,
        "table_03_root_inference.csv": 12,
        "table_04_special_power.csv": 6,
        "table_05_tail_strata.csv": 5,
    }
    for name, count in expected_output_rows.items():
        path = output_directory / name
        if path.is_file():
            audit.require(len(csv_rows(path)) == count, f"Wrong result row count: {name}")
    for name in SIMULATION_FIGURES:
        path = output_directory / name
        if path.is_file():
            audit.require(
                png_dimensions(path) == (1800, 1400),
                f"Reference figure has wrong dimensions: {name}",
            )

    tail = csv_rows(output_directory / "table_05_tail_strata.csv")
    audit.require(
        [row["observed_rare_draws"] for row in tail]
        == ["0", "1", "2", "3--5", "6 or more"],
        "Tail-stratum order changed",
    )
    audit.require(
        sum(int(row["replications"]) for row in tail) == 2032,
        "Tail-stratum replication total changed",
    )


def verify_hillstrom(audit: Audit) -> None:
    directory = ROOT / "results" / "hillstrom"
    for name in HILLSTROM_FILES:
        audit.require((directory / name).is_file(), f"Missing Hillstrom result: {name}")
    summary = csv_rows(directory / "application_summary.csv")
    audit.require(len(summary) == 1, "Hillstrom summary must have one row")
    if len(summary) == 1:
        row = summary[0]
        audit.require(row["n_x"] == "21306", "Hillstrom n_x changed")
        audit.require(row["n_y"] == "21307", "Hillstrom n_y changed")
        audit.require(row["bootstrap_reps"] == "9999", "Hillstrom B changed")
        audit.require(row["bootstrap_grid_size"] == "251", "Hillstrom grid changed")
        audit.require(row["audit_grid_size"] == "1001", "Hillstrom audit grid changed")
        audit.require(
            row["conclusion"] == "UNIFORM_Y_DOMINANCE_NO_ROOT",
            "Hillstrom conclusion changed",
        )
        audit.require(row["no_root_on_K_certified"] == "TRUE", "No-root certificate changed")
        audit.require(
            row["uniform_y_dominance_certified"] == "TRUE",
            "Uniform-dominance certificate changed",
        )
        audit.require(row["outer_component_count"] == "0", "Outer set changed")

    expected_rows = {
        "data_validation.csv": 15,
        "descriptive_statistics.csv": 2,
        "moment_band.csv": 251,
        "moment_band_audit.csv": 1001,
        "anchor_moments.csv": 6,
        "enclosure_cells.csv": 500,
        "outer_set_components.csv": 0,
        "certified_sign_intervals.csv": 1,
        "certified_brackets.csv": 0,
        "tail_diagnostics.csv": 24,
        "bootstrap_suprema.csv": 9999,
        "audit_checks.csv": 3,
        "rng_state_hashes.csv": 2,
        "table_06_hillstrom_application.csv": 8,
    }
    for name, count in expected_rows.items():
        path = directory / name
        if path.is_file():
            audit.require(len(csv_rows(path)) == count, f"Wrong Hillstrom row count: {name}")
    validation = csv_rows(directory / "data_validation.csv")
    audit.require(
        all(row["status"] == "PASSED" for row in validation),
        "A Hillstrom data-validation check failed",
    )
    checks = csv_rows(directory / "audit_checks.csv")
    audit.require(
        all(row["status"] == "PASSED" for row in checks),
        "A Hillstrom continuum audit failed",
    )
    descriptive = csv_rows(directory / "descriptive_statistics.csv")
    audit.require(
        [(row["role"], row["n"], row["n_positive"]) for row in descriptive]
        == [("X", "21306", "122"), ("Y", "21307", "267")],
        "Hillstrom descriptive counts changed",
    )
    bootstrap = csv_rows(directory / "bootstrap_suprema.csv")
    audit.require(
        [int(row["replication"]) for row in bootstrap] == list(range(1, 10000)),
        "Hillstrom bootstrap indices changed",
    )
    figure = directory / "figure_03_hillstrom_moment_contrast.pdf"
    audit.require(figure.is_file(), "Hillstrom figure is missing")
    if figure.is_file():
        audit.require(figure.read_bytes()[:4] == b"%PDF", "Hillstrom figure is not PDF")


def png_dimensions(path: Path) -> tuple[int, int] | None:
    data = path.read_bytes()[:24]
    if len(data) != 24 or data[:8] != b"\x89PNG\r\n\x1a\n":
        return None
    return struct.unpack(">II", data[16:24])


def compare_simulation(audit: Audit, generated: Path) -> None:
    expected = ROOT / "results" / "simulation"
    for name in SIMULATION_FILES:
        actual = generated / name
        audit.require(actual.is_file(), f"Generated simulation file is missing: {name}")
        if actual.is_file():
            audit.require(
                digest(actual) == digest(expected / name),
                f"Generated simulation values differ: {name}",
            )
    for name in SIMULATION_FIGURES:
        actual = generated / name
        audit.require(actual.is_file(), f"Generated figure is missing: {name}")
        if actual.is_file():
            audit.require(
                png_dimensions(actual) == (1800, 1400),
                f"Generated figure has wrong dimensions: {name}",
            )


def compare_hillstrom(audit: Audit, generated: Path) -> None:
    expected = ROOT / "results" / "hillstrom"
    for name in HILLSTROM_FILES:
        actual = generated / name
        audit.require(actual.is_file(), f"Generated Hillstrom file is missing: {name}")
        if actual.is_file():
            audit.require(
                digest(actual) == digest(expected / name),
                f"Generated Hillstrom values differ: {name}",
            )
    figure = generated / "figure_03_hillstrom_moment_contrast.pdf"
    audit.require(figure.is_file(), "Generated Hillstrom figure is missing")
    if figure.is_file():
        audit.require(figure.read_bytes()[:4] == b"%PDF", "Generated figure is not PDF")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--compare-simulation", type=Path)
    parser.add_argument("--compare-hillstrom", type=Path)
    args = parser.parse_args()

    audit = Audit()
    verify_manifest_and_privacy(audit)
    verify_simulation(audit)
    verify_hillstrom(audit)
    if args.compare_simulation:
        compare_simulation(audit, args.compare_simulation.resolve())
    if args.compare_hillstrom:
        compare_hillstrom(audit, args.compare_hillstrom.resolve())

    if audit.errors:
        print("PUBLIC REPRODUCIBILITY CHECK: FAILED", file=sys.stderr)
        for error in audit.errors:
            print(f"- {error}", file=sys.stderr)
        return 1
    print(f"PUBLIC REPRODUCIBILITY CHECK: PASSED ({audit.checks} checks)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
