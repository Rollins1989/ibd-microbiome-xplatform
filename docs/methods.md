# Methods

## 1. 16S branch
Primers (optional) are removed with cutadapt (`--discard-untrimmed`). Reads are denoised **per sequencing batch**
with DADA2 (`filterAndTrim` with configured truncation and expected-error limits, run-specific error models via
`learnErrors`, `dada`, `mergePairs`); batch sequence tables are merged, short sequences and bimeras
(`removeBimeraDenovo`, consensus) are removed, and taxonomy is assigned with the naive Bayesian classifier
(`assignTaxonomy`, `minBoot` configurable; optional `addSpecies`). Sequence variants are named by abundance rank.
Reagent contaminants are flagged with `decontam` (prevalence method) using negative-control samples. The ASV
sequences are aligned (MAFFT `--auto`), a GTR+Γ tree is built (FastTree) and midpoint-rooted (phangorn), and
MetaCyc pathway abundances are predicted with PICRUSt2. Samples below `min_reads_per_sample` are excluded.
ASVs are collapsed to genus for differential abundance and platform comparison (unresolved genera become
`Unclassified_<lowest rank>`).

## 2. Shotgun branch
fastp (adapter detection, tail trimming, length/low-complexity filters) → bowtie2 `--very-sensitive-local` against
GRCh38, keeping pairs with both mates unmapped. MetaPhlAn 4 is run with `-t rel_ab_w_read_stats` (relative
abundance **and** estimated clade read counts; counts feed count-based methods) and re-run from the saved bowtie2
output with `-t rel_ab` for HUMAnN. HUMAnN uses the MetaPhlAn profile to build sample-specific pangenomes; pathway
abundance is renormalised to CPM and split into stratified/unstratified tables (UNMAPPED/UNINTEGRATED rows are
removed before analysis). StrainPhlAn markers are extracted from the MetaPhlAn SAM output and trees are built for
the configured species-level genome bins. Optional Kraken2/Bracken profiles allow a profiler-agreement check.
Samples below `min_microbial_reads` (estimated MetaPhlAn reads) are excluded.

## 3. MAGs
Per selected sample: MEGAHIT (contigs ≥ `min_contig_len`), read mapping (bowtie2) and depth, three binners
(MetaBAT2, MaxBin2, CONCOCT with contig cutting), consolidation by DAS Tool, quality by CheckM2. Bins with
completeness ≥ 50 % and contamination ≤ 10 % (MIMAG medium; ≥ 90/< 5 flagged high quality) are dereplicated with
dRep (primary 90 %, secondary 95 % ANI, quality from CheckM2), classified with GTDB-Tk, annotated with Bakta and
screened with RGI (CARD; Perfect/Strict hits) and ABRicate (NCBI, VFDB). A de novo bacterial tree comes from the
GTDB-Tk user alignment (FastTree LG+Γ, midpoint rooted). MAG abundance in **every** shotgun sample is measured
with CoverM (≥ 95 % identity, ≥ 75 % aligned, ≥ 10 % covered). A MAG is "novel" when GTDB-Tk assigns no species;
"not in the profiler" when its GTDB species name is absent from all MetaPhlAn species tables (string match after
normalisation: a conservative, approximate flag).

## 4. Statistics
**Layers.** Each analysis consumes standardised feature × sample tables: 16S genus (counts), shotgun genus and
species (estimated read counts), 16S and shotgun pathways (abundance), MAG abundance. Features are filtered to
prevalence ≥ 10 % (and mean proportion ≥ 0.01 % for abundance layers) inside each analysis.

**Alpha/beta diversity.** Observed richness, Shannon and (16S) Faith's PD. Alpha diversity ~ diagnosis + covariates
+ (1 | subject) (Satterthwaite tests). Beta diversity: Bray–Curtis, Aitchison (Euclidean on CLR), UniFrac (16S).
PERMANOVA (`adonis2`, marginal terms) on one baseline sample per subject, adjusted for covariates and batch,
overall and pairwise (Benjamini–Hochberg); `betadisper` + `permutest` for dispersion; sensitivity over repeated random
one-sample-per-subject draws. (Using `strata = subject` for a between-subject factor is invalid, and rows with
missing covariates are removed explicitly.) Rarefaction is not used; a rarefied Shannon comparison is a sensitivity
check only.

**Differential abundance.** Four methods when available, always with the same covariates:
CLR + linear mixed model (lmerTest), MaAsLin2 (CLR, LM with random subject effect), ANCOM-BC2 (random intercept),
ALDEx2 (baseline samples, Welch). Benjamini–Hochberg within method × contrast. Consensus: q < `q_threshold` in at
least `min_methods` methods (capped by the number of methods run) with identical direction. Effect sizes shown are
those of the primary method (CLR-LMM).

**Platform concordance** (specimens with both platforms). Genus names are normalised (strip rank prefixes, brackets,
GTDB suffixes, group labels; unresolvable names removed) and the covered abundance is reported. Procrustes
(`protest`) and Mantel (Spearman) on Aitchison distances, per-genus Spearman correlation, Bland–Altman on Shannon
diversity of shared genera, and per-pathway Spearman plus Mantel between PICRUSt2 and HUMAnN pathway profiles.

**Longitudinal.** Within-subject association of disease activity with CLR abundance (mixed model with diagnosis and
covariates) and trajectories aligned to flare onset (visits −3…+2 by default), pooled across events.

**Prediction.** Binary tasks (IBD vs non-IBD, CD vs non-IBD, UC vs non-IBD, CD vs UC, pre-flare). Pre-flare:
specimens that are inactive and whose next visit (≤ `max_gap_weeks`) is active are positive; those whose next visit
is inactive are negative. Random forest, XGBoost and elastic-net logistic regression; outer StratifiedGroupKFold
(subjects never split) repeated with different seeds, inner grouped CV for tuning; prevalence filter + CLR fitted per
fold. Out-of-fold predictions averaged over repeats give pooled ROC curves; 95 % intervals by subject-level bootstrap.
Feature sets: 16S genus, shotgun species, both (on the same specimens). SHAP (TreeExplainer) on the best-tuned model
refit on all data.

**Networks.** SparCC (median of Dirichlet resamples, 10 exclusion iterations); empirical FDR from a pooled
permutation null; edges need |r| ≥ `min_abs_r` and FDR ≤ `fdr`; Louvain modules; same features across groups;
comparison by edge overlap (Jaccard), sign flips, density and modularity.

## 5. Deviations from the original project outline
| Outline | Implemented | Why |
|---|---|---|
| `adonis2(..., strata = subject)` | baseline-per-subject PERMANOVA + resampling sensitivity | within-subject permutation cannot test a between-subject factor |
| ConQuR batch correction | batch as covariate in all models; run-specific DADA2 error models | ConQuR is not reliably installable via conda; modelling avoids altering data |
| SpiecEasi or SparCC | in-house SparCC with permutation-FDR, unit-tested | dependency risk; algorithm is short and verifiable |
| Kraken2/Bracken as main profiler | MetaPhlAn 4 primary; Kraken2/Bracken optional | comparable cost/benefit; optional agreement analysis |
| Prokka/Bakta | Bakta only | Prokka is no longer maintained |
| MOFA2/DIABLO, metatranscriptomics, pangenomes, LSTM/Cox, trajectory clustering | not implemented | listed as optional extensions |
