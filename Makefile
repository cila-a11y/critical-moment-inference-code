SHELL := bash
CPUS ?= 2
ASSET ?= main_final_v003_replication_results.rds
DATA ?= data/hillstrom_2008_raw.csv
SOURCE ?=

.PHONY: verify verify-depth rebuild reproduce reproduce-depth simulate data hillstrom

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
