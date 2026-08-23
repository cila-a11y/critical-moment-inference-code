# Frozen enclosure-depth audit

This directory records the fully paired diagnostic used to assess whether
increasing the maximum continuum-enclosure depth from six to eight changes the
reported scientific conclusions. The relative-variance guard remains fixed at
depth six in both variants; all samples, multiplier streams, grids, bootstrap
replications, node limits, and remaining inferential settings are paired and
unchanged.

The audit contains six scientific cells and 254 independent pairs per cell:
1,524 paired replications, or 3,048 variant evaluations. The depth-eight outer
set is contained in the depth-six outer set in every pair. No population root
is lost among the 1,454 pairs with continuum-band coverage, all 1,524
length-bound checks pass, and the prespecified certification, isolation,
reporting, and coverage decisions agree in every pair. There are no numerical,
node-limit, or variance-limit failures.

Across the six cells, the mean per-pair increase ranges from 8.51 to 18.93
enclosure nodes and from 0.0127 to 0.0575 seconds. The largest absolute change
in mean recorded R-heap use is 0.322 percent of the corresponding reference
mean.

The study was executed with source version `2.0.0-bridge.7` and the exact
configuration in `simulation/config/enclosure_depth_v001.R`, whose SHA-256 is
`1ee6977fd82e0852f60e25b9afdc1c094d9461c8e7fe2bed8179bbfa15bfcd7e`.
The public design, data-generating, inferential, execution, reporting, and
validation source files are byte-identical to those used for the frozen run.
The public summarizer subsequently added reporting-only fields and
deduplicated tables; reaggregation therefore compares every historical field
while allowing these additional columns.

The replication-level RDS is committed under
`inputs/simulation/enclosure_depth_v001/`. Run `make verify-depth` for the
frozen scientific audit and `make reproduce-depth` to reaggregate its CSV
outputs with R 4.4.2. The historical
`sensitivity_paired_comparisons.csv` alias is omitted because it is byte for
byte identical to `paired_comparisons.csv`.
