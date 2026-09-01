# Critical moment inference: reproducibility code

This repository contains the code, locked configurations, aggregate reference
inputs, and reference outputs needed to reproduce the numerical tables and
figures in *Simultaneous nonparametric inference and structural root-count
certification for critical moment orders*.

## Requirements

- R 4.4.2 with Cairo graphics support;
- Python 3;
- Bash, GNU `sha256sum`, `flock`, and `make`;
- `curl` only for downloading the Hillstrom data.
- GitHub CLI only for the one-line release-asset download shown below.

The verified Hillstrom numerical audit additionally requires a C11 compiler
and the MPFR and GMP development headers and libraries. Its supported target
is Linux, WSL, or a Linux container; native Windows is not supported.

The analysis uses base and recommended R packages only.

## Verify the public package

From the repository root:

```sh
make verify
```

This checks every committed file against `SHA256SUMS.txt`, rejects private or
unexpected artifacts, and verifies the locked simulation and application
contracts.

## Simulation materials for Tables 1--7 and Figures 1--2

Tables 1 and 2 describe the data-generating mechanisms and the locked Monte
Carlo cell catalogue specified in `simulation/R/dgp.R` and
`simulation/config/main_final_v003.R`. The fast rebuild uses the committed
aggregate inputs to regenerate the numerical source files for Tables 3--7
and Figures 1--2:

```sh
make rebuild
```

Outputs are written to `generated/rebuild/`. Table values are checked against
the committed references. The stable pipeline filenames map to the manuscript
as follows:

- `table_02_benchmark.csv`: Table 3;
- `table_03_root_inference.csv`: Table 4;
- `table_04_special_power.csv`: Table 5;
- `table_05_tail_strata.csv`: Table 6;
- `table_07_ablation.csv`: Table 7.

`table_01_design.csv` is a compact design summary. PNG rendering is checked
at the locked dimensions of 1800 by 1400 pixels.

For an independent reaggregation, download
`main_final_v003_replication_results.rds`, which is intentionally stored in
release `v1.0.0`, into the repository root, then run:

```sh
gh release download v1.0.0 \
  --repo cila-a11y/critical-moment-inference-code \
  --pattern main_final_v003_replication_results.rds
make reproduce
```

The asset contains 284,480 replication-level result rows, is 36,063,100 bytes,
and has SHA-256
`eabac798fb2da2a126203d965e5de41012af7aad3b74358ee70302b8dd35606b`.
Later source releases do not duplicate this asset. The script verifies its
row count, byte count, and digest before reading it.

To rerun all 109,728 Monte Carlo units under the locked design:

```sh
make simulate CPUS=2
```

The full run comprises 48 resumable tasks and writes checkpoints only under
`generated/`. The stages may also be invoked separately: first
`bash scripts/run_full_simulation.sh prepare`, then the `task` command once for each
integer from 1 through 48, and finally
`bash scripts/run_full_simulation.sh finalize`.

## Verify and reaggregate the enclosure-depth audit

The frozen fully paired audit compares maximum continuum-enclosure depths six
and eight while keeping the relative-variance guard at depth six. It contains
1,524 independent pairs, corresponding to 3,048 variant evaluations. Its
scientific and numerical invariants can be checked without R:

```sh
make verify-depth
```

The committed replication-level asset is 405,328 bytes and has SHA-256
`f66e1ba64b37c8f76b2f30ebef554187e68c8cfcf555d8bbe04d4e2f7ee27ecd`.
With R 4.4.2, reaggregate the frozen rows and compare every historical output
field against the committed results with:

```sh
make reproduce-depth
```

The exact configuration and a concise account of the audit are stored with
the frozen results under `results/simulation/enclosure_depth_v001/`.

## Reproduce Table 8 and Figure 3

The Hillstrom CSV was released with the
[2008 MineThatData E-Mail Analytics Challenge](https://blog.minethatdata.com/2008/03/minethatdata-e-mail-analytics-and-data.html).
It is not redistributed here because no formal upstream redistribution licence
was located. Because the original download URL is no longer available, the
acquisition script uses a
[commit-pinned public archival copy](https://github.com/AdityaDabrase/ab-testing-email-marketing/blob/5866d7e8b50f46239c80d0ffa543fda501939ecc/data/raw/hillstrom.csv).
It downloads or accepts only the exact locked original file: 3,964,977 bytes
with SHA-256
`0e5893329d8b93cefecc571777672028290ab69865718020c78c7284f291aece`.

```sh
make data
make hillstrom
```

Alternatively, place the original file at `data/hillstrom_2008_raw.csv`;
`make hillstrom` still requires the locked byte count and SHA-256 before any
computation.

The application compares the `No E-Mail` and `Mens E-Mail` groups, runs the
locked 9,999-repetition simultaneous multiplier procedure, reconstructs the
continuum enclosure, and checks the generated scientific files against the
committed references in `results/hillstrom/`.
Numeric CSV fields are compared at the same relative tolerance of `1e-11`
used by the independent R post-run contract; schemas, row order, replication
indices, and categorical conclusions must agree exactly. RNG-state hashes are
checked within each run because R's version-3 serialisation records
locale-dependent encoding metadata; cross-system reproduction is checked
through all 9,999 indexed bootstrap suprema and the derived outputs.
The stable pipeline file `table_06_hillstrom_application.csv` contains the
descriptive and implemented-band values for manuscript Table 8. The verified
whole-cell column is supplied by
`verified_numerical_audit_anchors.csv`.

### Run the verified 4,000-cell numerical audit

After placing the locked Hillstrom CSV at `data/hillstrom_2008_raw.csv`, run:

```sh
make audit-hillstrom
```

This is separate from `make hillstrom`: it does not rerun the multiplier
bootstrap or require R. It derives the critical value from the 9,999 committed
multiplier suprema, compiles the interval-arithmetic audit, and checks all
4,000 cells at 256-bit MPFR precision with directed rounding. The generated
summary and six Table 8 anchors are compared byte for byte with
`results/hillstrom/verified_numerical_audit_summary.csv` and
`results/hillstrom/verified_numerical_audit_anchors.csv`.
The recorded decimal endpoints are rounded outwards from the directed MPFR
enclosures.

On Debian or Ubuntu, the additional packages are `build-essential`,
`libmpfr-dev`, `libgmp-dev`, and `pkg-config`. For a non-system installation,
set both `CMO_MPFR_INCLUDEDIR` and `CMO_MPFR_LIBDIR` to its include and library
directories. The audit first verifies the locked Hillstrom byte count and
SHA-256; the raw data are not redistributed.

### Rebuild the manuscript-specific Figure 3

The JNS manuscript uses a presentation-only vector redraw of Figure 3. It is
built from the already committed `results/hillstrom/moment_band.csv`; no data
or numerical result is recomputed. From the repository root, run:

```sh
make figure3-jns
```

The rebuilt file is written to
`generated/jns_figure03/figures/figure_03_hillstrom_moment_contrast_jns.pdf`.
The exact vector file embedded in the submitted manuscript is frozen at
`results/hillstrom/figure_03_hillstrom_moment_contrast_jns.pdf`, and the exact
historical plotting script is `scripts/build_figure03_jns.R`. The frozen PDF
has SHA-256
`9ba4216292a2b086974526728052b7e8fa205197d526b2769988b2cdcf55cc73`.

The rebuilt PDF can differ byte for byte because Cairo records its version and
creation time in PDF metadata. The plotted values and presentation are fixed
by the committed CSV and script; `make verify` checks the frozen manuscript
artifact byte for byte through `SHA256SUMS.txt`.

## Contents

- `simulation/`: simulation engine and locked main and enclosure-depth designs;
- `application/`: locked Hillstrom application;
- `inputs/simulation/`: aggregate inputs and frozen replication-level depth audit;
- `results/`: reference tables and figures;
- `scripts/`: verification and reproduction entry points.

The manuscript, raw data, machine-specific job files, logs, and administrative
archives are intentionally excluded.
