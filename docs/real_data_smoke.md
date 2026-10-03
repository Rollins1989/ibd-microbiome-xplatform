# Real-data smoke test

The repository now includes a lightweight GitHub Actions smoke test that executes a real workflow tool on public sequencing data.

## Dataset

The smoke test downloads one paired-end sample from the public nf-core `ampliseq` test dataset:

- `1_S103_L001_R1_001.fastq.gz`
- `1_S103_L001_R2_001.fastq.gz`

These are real 16S rRNA V4 reads from the ZymoBIOMICS Microbial Community Standard. The dataset is a small CI-scale mock community with known composition and is described as real sequencing data in the upstream test-dataset workflow.

## Executed tool

GitHub Actions installs **Cutadapt** and runs the same type of paired-end primer-removal operation used by the workflow's `cutadapt_16s` rule, using the V4 primers:

- Forward: `GTGYCAGCMGCCGCGGTAA`
- Reverse: `GGACTACNVGGGTWTCTAAT`

The smoke test verifies that Cutadapt produces non-empty paired FASTQ outputs and that both compressed outputs pass `gzip -t`.

## Scope

This is intentionally a small real-data/tool smoke test rather than a full biological analysis. It demonstrates that a real public sequencing input can reach a real preprocessing tool successfully without requiring the large reference databases needed by the complete workflow.

The main end-to-end validation remains the synthetic cohort because it provides deterministic truth-recovery checks. The real-data smoke test complements it by exercising an actual sequencing-tool path on public reads.
