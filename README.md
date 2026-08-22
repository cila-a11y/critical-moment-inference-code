# Critical moment inference: reproducibility code

This repository contains the code, locked configurations, aggregate reference
inputs, and reference outputs needed to reproduce the numerical tables and
figures in *Inference for critical moment orders under tail-crossing
constraints*.

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

## Contents

- `simulation/`: simulation engine and locked `main_final_v003` design;
- `application/`: locked Hillstrom application;
- `inputs/simulation/`: aggregate inputs for the fast rebuild;
- `results/`: reference tables and figures;
- `scripts/`: verification and reproduction entry points.

The manuscript, raw data, machine-specific job files, logs, and administrative
archives are intentionally excluded.
