# What has been validated, and how

**Honest summary:** the statistical/analysis layer and all glue code were executed and checked on simulated data;
third-party bioinformatics tools were *not* executed (no databases or tool installs in the build environment). The
rules for those tools are syntactically validated (Snakemake DAG, rendered shell commands reviewed) but their first
run on real data may need parameter or version adjustments.

| Component | Status | Evidence |
|---|---|---|
| Sample-sheet validation, flare/activity derivation, CLR, taxon-name harmonisation, SparCC | **Executed and tested** | `tests/test_units.py` (hand-checked flare labels, CLR invariants, planted SparCC correlations recovered) |
| MetaPhlAn merge, CheckM2/GTDB-Tk/RGI/ABRicate/CoverM/Bracken parsers, bin collection, HMP2 sheet builder | **Executed on mimicked tool output** | `tests/test_parsers.py`; parsers raise explicit errors on format drift. Real-output formats written from tool documentation – verify on first run |
| DADA2 scripts (per-batch denoising, merge, chimeras, taxonomy) | **Executed on simulated amplicon reads** | recovered all 6 true sequences with correct taxonomy; two batches, pseudo-pooling |
| decontam, diversity (alpha/beta/UniFrac/PERMANOVA/betadisper), CLR-LMM DA, consensus, concordance, longitudinal, networks, nested-CV ML + SHAP, MAG analysis, strain analysis, figures, report | **Executed end-to-end via Snakemake** on a simulated cohort (64 subjects, 383 specimens, 281 paired) | `make test`: 25 checks of recovered planted signal (all pass); ML label-permutation control (CD vs non-IBD, 16S, logistic regression): mean AUROC 0.52 over 12 seeds |
| MaAsLin2, ANCOM-BC2, ALDEx2 wrappers | **Written to the packages' documented APIs; not executed** (packages unavailable in the build environment) | They share the tested input preparation and output schema with CLR-LMM; the first real run should be checked against `da/*.clr_lmm.tsv` |
| Rules for fastp, bowtie2, MetaPhlAn, HUMAnN, StrainPhlAn, Kraken2/Bracken, MEGAHIT, binners, DAS Tool, CheckM2, dRep, GTDB-Tk, Bakta, RGI, ABRicate, CoverM, MAFFT, FastTree, PICRUSt2 | **DAG-validated only** (`make dryrun`: 96 jobs; 113 with all optional modules) | commands rendered and inspected; conda pins are unverified – adjust if resolution fails |

Simulation caveat: simulated cohorts verify that the code is *correct and sensitive*; they say nothing about the
biology of a real cohort, where effect sizes are smaller and confounding is real.
