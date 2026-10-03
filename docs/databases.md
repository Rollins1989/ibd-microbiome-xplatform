# Databases and how to obtain them

All paths go into the `databases:` block of `config/config.yaml`. Commands are the tools' standard ones; check each
tool's current documentation, as download locations and versions change.

| Resource | How |
|---|---|
| HMP2 metadata + data | https://ibdmdb.org (downloads) and ENA/SRA BioProject **PRJNA398089**. Build the sheet with `workflow/scripts/build_sample_sheet_hmp2.py` (curl command in its docstring) |
| Human reference (host removal) | GRCh38 FASTA from NCBI/Ensembl → `bowtie2-build GRCh38.fa GRCh38` → `host_bowtie2_index` = prefix |
| 16S training set | DADA2-formatted SILVA/GTDB references: see https://benjjneb.github.io/dada2/training.html (SILVA) or a GTDB DADA2-formatted release. Using a **GTDB** classifier makes 16S genus names compatible with MetaPhlAn/GTDB-Tk (the concordance step harmonises residual differences and reports the coverage) |
| MetaPhlAn 4 | `metaphlan --install --index <index> --bowtie2db <dir>`; `metaphlan_index` must match the one HUMAnN/StrainPhlAn expect |
| HUMAnN | `humann_databases --download chocophlan full <dir>`; `humann_databases --download uniref uniref90_diamond <dir>`. **Check the HUMAnN release notes for the MetaPhlAn index it supports** and set `metaphlan_index` accordingly |
| Kraken2 / Bracken (optional) | `kraken2-build --standard --db <dir>` then `bracken-build -d <dir> -t <threads> -k 35 -l <read length>` |
| CheckM2 | `checkm2 database --download --path <dir>` → `checkm2` = path to the `.dmnd` file |
| GTDB-Tk | `download-db.sh` shipped with GTDB-Tk (release matching your GTDB-Tk version) → `gtdbtk` = data directory |
| Bakta | `bakta_db download --output <dir> --type full` |
| RGI / CARD | install the CARD database with `rgi load` (see RGI docs) before the first run |
| ABRicate | `abricate --setupdb` (uses the `ncbi` and `vfdb` databases) |
| StrainPhlAn clades | species-level genome bins (SGBs) are listed in the database folder of your MetaPhlAn index; put ids such as `t__SGB<number>` in `shotgun.strainphlan.clades` |

Resource planning (rough): MetaPhlAn ≈ 2 CPU-h and HUMAnN ≈ 10–20 CPU-h per sample; MEGAHIT + binning ≈ 50–150
CPU-h per sample; GTDB-Tk needs ≈ 100 GB RAM. Use `assemble` to limit assembly to a subset.
