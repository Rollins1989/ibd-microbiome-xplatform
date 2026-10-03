#!/usr/bin/env Rscript
# Merge batch sequence tables, remove chimeras, assign taxonomy, write ASV table / FASTA / taxonomy / tracking.
suppressPackageStartupMessages({library(optparse); library(dada2)})
opt <- parse_args(OptionParser(option_list = list(
  make_option("--seqtabs"), make_option("--tracks"), make_option("--train-set"), make_option("--species-set", default = ""),
  make_option("--min-boot", type = "integer", default = 50L), make_option("--min-length", type = "integer", default = 100L),
  make_option("--threads", type = "integer", default = 1L), make_option("--seed", type = "integer", default = 1L),
  make_option("--outdir"))), convert_hyphens_to_underscores = FALSE)
set.seed(opt$seed)
mt <- if (opt$threads > 1) opt$threads else FALSE
tabs <- lapply(strsplit(opt$seqtabs, ",")[[1]], readRDS)
st <- if (length(tabs) == 1) tabs[[1]] else do.call(mergeSequenceTables, c(tabs, list(repeats = "error")))
st <- st[, nchar(colnames(st)) >= opt$`min-length`, drop = FALSE]
nochim <- removeBimeraDenovo(st, method = "consensus", multithread = mt, verbose = FALSE)
message(sprintf("chimeras removed: %.1f%% of reads", 100 * (1 - sum(nochim) / sum(st))))
nochim <- nochim[, order(colSums(nochim), decreasing = TRUE), drop = FALSE]
asv_ids <- sprintf("ASV%05d", seq_len(ncol(nochim)))
tax <- assignTaxonomy(colnames(nochim), opt$`train-set`, multithread = mt, minBoot = opt$`min-boot`, tryRC = TRUE,
                      taxLevels = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"))
tax <- as.data.frame(tax, stringsAsFactors = FALSE)
tax$Species <- NA_character_
if (nzchar(opt$`species-set`)) {
  sp <- addSpecies(as.matrix(tax[, 1:6]), opt$`species-set`, tryRC = TRUE)
  tax$Species <- ifelse(is.na(sp[, "Species"]), NA, paste(sp[, "Genus"], sp[, "Species"]))
}
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)
tab <- t(nochim); rownames(tab) <- asv_ids
write.table(cbind(ASV = asv_ids, as.data.frame(tab, check.names = FALSE)), file.path(opt$outdir, "asv_table.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(c(rbind(paste0(">", asv_ids), colnames(nochim))), file.path(opt$outdir, "asv_seqs.fasta"))
write.table(cbind(ASV = asv_ids, tax), file.path(opt$outdir, "taxonomy.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
trk <- do.call(rbind, lapply(strsplit(opt$tracks, ",")[[1]], read.delim))
trk$nonchim <- rowSums(nochim)[trk$sample_id]
write.table(trk, file.path(opt$outdir, "read_tracking.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
message(sprintf("ASVs: %d; classified to genus: %.1f%% of ASVs", ncol(nochim), 100 * mean(!is.na(tax$Genus))))
