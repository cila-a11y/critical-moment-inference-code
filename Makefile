SHELL := bash
CPUS ?= 2
ASSET ?= main_final_v003_replication_results.rds
DATA ?= data/hillstrom_2008_raw.csv
SOURCE ?=

.PHONY: verify rebuild reproduce simulate data hillstrom

verify:
	python3 scripts/verify.py

rebuild:
	bash scripts/rebuild_simulation_outputs.sh

reproduce:
	bash scripts/reproduce_from_replications.sh "$(ASSET)"

simulate:
	CPUS="$(CPUS)" bash scripts/run_full_simulation.sh all

data:
	bash scripts/acquire_hillstrom.sh "$(DATA)" "$(SOURCE)"

hillstrom:
	bash scripts/run_hillstrom.sh "$(DATA)"
