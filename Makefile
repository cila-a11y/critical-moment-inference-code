SHELL := bash
CPUS ?= 2
# Replication-level asset stored in GitHub release v1.0.0; see README.md.
ASSET ?= main_final_v003_replication_results.rds
DATA ?= data/hillstrom_2008_raw.csv
SOURCE ?=

.PHONY: verify verify-depth rebuild reproduce reproduce-depth simulate data hillstrom audit-hillstrom figure3-jns

verify:
	python3 scripts/verify.py
	python3 scripts/verify_enclosure_depth.py

verify-depth:
	python3 scripts/verify_enclosure_depth.py

rebuild:
	bash scripts/rebuild_simulation_outputs.sh

reproduce:
	bash scripts/reproduce_from_replications.sh "$(ASSET)"

reproduce-depth:
	bash scripts/reproduce_enclosure_depth.sh

simulate:
	CPUS="$(CPUS)" bash scripts/run_full_simulation.sh all

data:
	bash scripts/acquire_hillstrom.sh "$(DATA)" "$(SOURCE)"

hillstrom:
	bash scripts/run_hillstrom.sh "$(DATA)"

audit-hillstrom:
	bash scripts/audit_hillstrom_numerics.sh "$(DATA)"

figure3-jns:
	mkdir -p generated/jns_figure03/figures
	cp results/hillstrom/moment_band.csv generated/jns_figure03/moment_band.csv
	Rscript --vanilla scripts/build_figure03_jns.R generated/jns_figure03
	test -s generated/jns_figure03/figures/figure_03_hillstrom_moment_contrast_jns.pdf
