#!/usr/bin/env Rscript
# Remove likely reagent/kit contaminant ASVs with decontam (prevalence method) using negative controls
# (sample_type == "negative_control" in the sample sheet). Without controls the table passes through.
# Control columns are always removed from the output.
suppressPackageStartupMessages(library(optparse))
.sd <- dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])); source(file.path(.sd, "utils.R"))
opt <- cli_args(list(make_option("--samples"), make_option("--asv"), make_option("--config"),
                     make_option("--out"), make_option("--removed")))
cfg <- load_cfg(opt$config)
sheet <- read_meta(opt$samples)
asv <- read_layer(opt$asv)
ctrl <- intersect(sheet$sample_id[sheet$sample_type == "negative_control" & sheet$data_type == "16S"], colnames(asv))
real <- intersect(sheet$sample_id[sheet$sample_type == "sample" & sheet$data_type == "16S"], colnames(asv))
removed <- data.frame(ASV = character(), p = numeric())
if (length(ctrl) >= 2) {
  suppressPackageStartupMessages(library(decontam))
  mat <- t(asv[, c(real, ctrl)])
  neg <- c(rep(FALSE, length(real)), rep(TRUE, length(ctrl)))
  res <- isContaminant(mat, neg = neg, method = "prevalence", threshold = cfg$amplicon$decontam$threshold)
  bad <- rownames(asv)[which(res$contaminant)]
  removed <- data.frame(ASV = bad, p = res$p[match(bad, rownames(asv))])
  log_msg(sprintf("decontam: %d negative controls, %d ASVs flagged", length(ctrl), length(bad)))
} else {
  log_msg("decontam: <2 negative controls available; no ASVs removed")
}
keep <- setdiff(rownames(asv), removed$ASV)
out <- asv[keep, real, drop = FALSE]
out <- out[rowSums(out) > 0, , drop = FALSE]
write_tsv(cbind(ASV = rownames(out), as.data.frame(out, check.names = FALSE)), opt$out)
write_tsv(removed, opt$removed)
