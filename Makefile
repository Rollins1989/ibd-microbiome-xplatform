.PHONY: test unit demo dryrun lint run clean
CORES ?= 8

unit:        ## unit + parser tests (seconds)
	python -m pytest -q tests

test:        ## simulated end-to-end test of the downstream analysis (minutes)
	CORES=$(CORES) bash tests/run_test.sh

demo:        ## lightweight recruiter demo: synthetic cohort, no external databases (minutes)
	python tests/simulate_data.py demo/sim --seed 7
	snakemake -s workflow/Snakefile --configfile config/demo.config.yaml --use-conda --cores 2


dryrun:      ## validate the complete DAG (all upstream + downstream rules) with the example sample sheet
	snakemake -s workflow/Snakefile -n --cores $(CORES) --config samples=config/samples.example.tsv outdir=results_dryrun

lint:
	snakemake -s workflow/Snakefile --lint --config samples=config/samples.example.tsv outdir=results_dryrun

run:         ## full run (needs databases, see docs/databases.md)
	snakemake -s workflow/Snakefile --use-conda --conda-frontend mamba --cores $(CORES) --rerun-incomplete --keep-going

clean:
	rm -rf results_dryrun tests/sim demo/sim demo/results .snakemake
