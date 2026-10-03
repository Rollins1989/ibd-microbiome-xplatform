#!/usr/bin/env bash
# End-to-end test of the downstream analysis on simulated data (no external databases or tools needed
# beyond Python/R with the packages in workflow/envs/py_stats.yaml and r_stats.yaml).
set -euo pipefail
cd "$(dirname "$0")/.."
python tests/simulate_data.py tests/sim
rm -rf tests/sim/results/{metadata,layers,da,diversity,concordance,longitudinal,networks,ml,strain,design,figures,report,config,logs} tests/sim/results/mags/analysis tests/sim/results/16s/asv_table.decontam.tsv tests/sim/results/16s/decontam_removed_asvs.tsv
python -m pytest -q tests
snakemake -s workflow/Snakefile --configfile config/test.config.yaml --cores "${CORES:-4}" --use-conda "$@"
python tests/check_results.py tests/sim/results
