# Critical moment inference: reproducibility code

This repository contains the code, locked configurations, aggregate reference
inputs, and reference outputs needed to reproduce the numerical tables and
figures in *Simultaneous nonparametric inference for critical moment orders
under prespecified tail-crossing bounds*.

## Requirements

- R 4.4.2 with Cairo graphics support;
- Python 3;
- Bash, GNU `sha256sum`, `flock`, and `make`;
- `curl` only for downloading the Hillstrom data.
- GitHub CLI only for the one-line release-asset download shown below.

The analysis uses base and recommended R packages only.

## Verify the public package

From the repository root:

```sh
make verify
```

This checks every committed file against `SHA256SUMS.txt`, rejects private or
unexpected artifacts, and verifies the locked simulation and application
contracts.

## Rebuild Tables 1--5 and Figures 1--2

The fast rebuild uses the committed aggregate simulation inputs:

```sh
make rebuild
```

Outputs are written to `generated/rebuild/` and checked against the committed
table data. PNG rendering is checked at the locked dimensions of
1800 by 1400 pixels.

For an independent reaggregation, download the release asset
`main_final_v003_replication_results.rds` from release `v1.0.0` into the
repository root, then run:

```sh
gh release download v1.0.0 --pattern main_final_v003_replication_results.rds
make reproduce
```

The asset contains 284,480 replication-level result rows, is 36,063,100 bytes,
and has SHA-256
`eabac798fb2da2a126203d965e5de41012af7aad3b74358ee70302b8dd35606b`.
The script verifies these facts before reading it.

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

## Reproduce Table 6 and Figure 3

The Hillstrom CSV is not redistributed because no formal upstream
redistribution licence was located. The acquisition script downloads or
accepts only the exact locked file: 3,964,977 bytes with SHA-256
`0e5893329d8b93cefecc571777672028290ab69865718020c78c7284f291aece`.

```sh
make data
make hillstrom
```

If the upstream server cannot be validated by the local TLS trust store, place
the official file at `data/hillstrom_2008_raw.csv`; `make hillstrom` still
requires the locked byte count and SHA-256 before any computation.

The application compares the `No E-Mail` and `Mens E-Mail` groups, runs the
locked 9,999-repetition simultaneous multiplier procedure, reconstructs the
continuum enclosure, and checks the generated scientific files against the
committed references in `results/hillstrom/`.

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
