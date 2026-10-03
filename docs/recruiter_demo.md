# Recruiter demo

This is the fastest way to see the workflow run without downloading large reference databases.

## Requirements

- Python 3.11+
- Snakemake 8.x
- The Python/R packages required by the lightweight analysis environment
- `make` (or run the commands below directly)

No GTDB-Tk, CheckM2, MetaPhlAn, HUMAnN, or other large external databases are required for this demo.

## Run it

From the repository root:

```bash
make demo
```

The demo will:

1. Generate a deterministic synthetic longitudinal IBD cohort with planted signals.
2. Run the analysis-only Snakemake DAG against that cohort.
3. Produce representative diversity, concordance, differential-abundance, ML, longitudinal, network, and report outputs.

The demo uses `config/demo.config.yaml`, which deliberately uses fewer permutations/bootstrap replicates than the full validation configuration so that it is suitable for a quick repository walkthrough.

## What the demo proves

The demo is a **software/reproducibility demonstration**, not a biological result. The synthetic cohort contains planted effects so the analysis code has a known signal to recover. For the full truth-recovery suite, run:

```bash
make test
```

For the distinction between simulated validation and third-party tool validation, see [`validation.md`](validation.md).

## Full pipeline

The production workflow can additionally execute raw-read QC, 16S processing, shotgun profiling, MAG recovery, genome annotation, AMR screening, and database-backed taxonomic analysis. Those stages require the configured reference databases and substantially more compute/storage; see [`databases.md`](databases.md).
