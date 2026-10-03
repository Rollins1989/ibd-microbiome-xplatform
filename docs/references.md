# Methods and software references

This workflow combines established methods for amplicon processing, shotgun profiling, contamination control, genome recovery, taxonomic classification, and microbiome association analysis. Cite the primary methods below when publishing results generated with the corresponding part of the workflow.

## Workflow and reproducibility

1. **Snakemake** — Mölder F, Jablonski KP, Letcher B, et al. (2021). *Sustainable data analysis with Snakemake*. F1000Research 10:33. https://doi.org/10.12688/f1000research.29032.3

## 16S / amplicon analysis

2. **DADA2** — Callahan BJ, McMurdie PJ, Rosen MJ, et al. (2016). *DADA2: High-resolution sample inference from Illumina amplicon data*. Nature Methods 13:581–583. https://doi.org/10.1038/nmeth.3869

3. **decontam** — Davis NM, Proctor DM, Holmes SP, Relman DA, Callahan BJ. (2018). *Simple statistical identification and removal of contaminant sequences in marker-gene and metagenomics data*. Microbiome 6:226. https://doi.org/10.1186/s40168-018-0605-2

4. **PICRUSt2** — Douglas GM, Maffei VJ, Zaneveld JR, et al. (2020). *PICRUSt2 for prediction of metagenome functions*. Nature Biotechnology 38:685–688. https://doi.org/10.1038/s41587-020-0548-6

## Shotgun profiling and association analysis

5. **MetaPhlAn 4** — Blanco-Míguez A, Beghini F, Cumbo F, et al. (2023). *Extending and improving metagenomic taxonomic profiling with uncharacterized species using MetaPhlAn 4*. Nature Biotechnology 41:1633–1644. https://doi.org/10.1038/s41587-023-01688-w

6. **MaAsLin 2** — Mallick H, Rahnavard A, McIver LJ, et al. (2021). *Multivariable association discovery in population-scale meta-omics studies*. PLOS Computational Biology 17:e1009442. https://doi.org/10.1371/journal.pcbi.1009442

7. **ANCOM-BC2** — Lin H, Peddada SD. (2023). *Multigroup analysis of compositions of microbiomes with covariate adjustments and repeated measures*. Nature Methods. https://doi.org/10.1038/s41592-023-02092-7

## MAG recovery and genome analysis

8. **MEGAHIT** — Li D, Liu C-M, Luo R, Sadakane K, Lam T-W. (2015). *MEGAHIT: an ultra-fast single-node solution for large and complex metagenomics assembly via succinct de Bruijn graph*. Bioinformatics 31:1674–1676. https://doi.org/10.1093/bioinformatics/btv033

9. **MetaBAT 2** — Kang DD, Li F, Kirton E, et al. (2019). *MetaBAT 2: an adaptive binning algorithm for robust and efficient genome reconstruction from metagenome assemblies*. PeerJ 7:e7359. https://doi.org/10.7717/peerj.7359

10. **CONCOCT** — Alneberg J, Bjarnason BS, de Bruijn I, et al. (2014). *Binning metagenomic contigs by coverage and composition*. Nature Methods 11:1144–1146. https://doi.org/10.1038/nmeth.3103

11. **DAS Tool** — Sieber CMK, Probst AJ, Sharrar A, et al. (2018). *Recovery of genomes from metagenomes via a dereplication, aggregation and scoring strategy*. Nature Microbiology 3:836–843. https://doi.org/10.1038/s41564-018-0171-1

12. **CheckM2** — Chklovski A, Parks DH, Woodcroft BJ, Tyson GW. (2023). *CheckM2: a rapid, scalable and accurate tool for assessing microbial genome quality using machine learning*. Nature Methods 20:1203–1212. https://doi.org/10.1038/s41592-023-01940-w

13. **dRep** — Olm MR, Brown CT, Brooks B, Banfield JF. (2017). *dRep: a tool for fast and accurate genomic comparisons that enables improved genome recovery from metagenomes through de-replication*. ISME Journal 11:2864–2868. https://doi.org/10.1038/ismej.2017.126

14. **GTDB-Tk** — Chaumeil P-A, Mussig AJ, Hugenholtz P, Parks DH. (2020). *GTDB-Tk: a toolkit to classify genomes with the Genome Taxonomy Database*. Bioinformatics 36:1925–1927. https://doi.org/10.1093/bioinformatics/btz848

15. **Bakta** — Schwengers O, Jelonek L, Dieckmann MA, et al. (2021). *Bakta: rapid and standardized annotation of bacterial genomes via alignment-free sequence identification*. Microbial Genomics 7:000685. https://doi.org/10.1099/mgen.0.000685

## Notes

- These references document the methods/tools used by the workflow; they do **not** imply that every third-party tool has been executed end-to-end in the repository's lightweight CI environment.
- See [`validation.md`](validation.md) for the exact distinction between simulated end-to-end validation, parser validation, and third-party tool/DAG validation.
- For software/database versioning, use the pinned Conda environment files and the workflow configuration in this repository.
