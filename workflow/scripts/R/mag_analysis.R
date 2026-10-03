#!/usr/bin/env Rscript
# MAG-level analysis: disease association of MAG abundance (CLR mixed model), novelty relative to GTDB,
# whether the reference-based profiler (MetaPhlAn) reports the MAG's species, AMR gene burden; Figure 6.
suppressPackageStartupMessages(library(optparse))
.sd <- dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])); source(file.path(.sd, "utils.R"))
suppressPackageStartupMessages({library(ape); library(parallel)})
opt <- cli_args(list(make_option("--config"), make_option("--meta"), make_option("--layer-dir"),
                     make_option("--mag-summary"), make_option("--amr"), make_option("--tree"),
                     make_option("--outdir"), make_option("--figdir"), make_option("--threads", type = "integer", default = 1)))
cfg <- load_cfg(opt$config); design <- cfg$design; st <- cfg$stats
meta <- read_meta(opt$meta)
mags <- read_meta(opt$`mag-summary`)
amr <- read_meta(opt$amr)
A <- read_layer(file.path(opt$`layer-dir`, "mag.tsv"))          # MAG x sample, % of reads
mpa <- read_layer(file.path(opt$`layer-dir`, "mgx_species.tsv"))
tree <- read.tree(opt$tree)

# ---- 1. association with diagnosis (abundance-type CLR, subject random effect)
md <- model_data(meta[meta$data_type == "MGX", ], design, colnames(A)); d <- md$data
M <- filter_layer(A[, rownames(d), drop = FALSE], st$min_prevalence, 0)
X <- clr_mat(M, "abundance")
one <- function(i) {
  ct <- tryCatch(coef_table(fit_lmm(X[i, rownames(d)], d, c("diagnosis", md$covars))), error = function(e) NULL)
  if (is.null(ct)) return(NULL)
  ct <- ct[grepl("^diagnosis", ct$term), ]
  data.frame(mag = rownames(X)[i], contrast = paste0(sub("^diagnosis", "", ct$term), "_vs_", levels(d$diagnosis)[1]),
             estimate = ct$estimate, se = ct$se, pvalue = ct$pvalue)
}
assoc <- do.call(rbind, mclapply(seq_len(nrow(X)), one, mc.cores = opt$threads))
assoc$qvalue <- ave(assoc$pvalue, assoc$contrast, FUN = function(p) p.adjust(p, "BH"))

# ---- 2. novelty, profiler coverage, AMR burden
sp_level <- sub(".*;s__", "", mags$gtdb_classification)
mags$novel_species <- !nzchar(sp_level)          # GTDB-Tk could not assign a named species (ANI < 95% to all references)
mags$species_label <- ifelse(mags$novel_species, "novel (no GTDB species match)", sp_level)
key_mags <- norm_species <- function(x) tolower(gsub("[^A-Za-z]", "", x))
detected <- norm_species(rownames(mpa)[rowSums(mpa > 0) >= 1])
mags$in_metaphlan_profile <- !mags$novel_species & norm_species(sp_level) %in% detected
amr_n <- table(factor(amr$mag, levels = mags$mag))
mags$n_amr_genes <- as.integer(amr_n[mags$mag])
mags$genus <- sub(";s__.*", "", sub(".*;g__", "", mags$gtdb_classification))
prev_by_dx <- sapply(levels(d$diagnosis), function(l) rowMeans(M[, d$diagnosis == l, drop = FALSE] > 0))
colnames(prev_by_dx) <- paste0("prevalence_", colnames(prev_by_dx))
out <- merge(mags, cbind(mag = rownames(prev_by_dx), as.data.frame(prev_by_dx)), by = "mag", all.x = TRUE)
write_tsv(out, file.path(opt$outdir, "mag_catalogue_annotated.tsv"))
write_tsv(assoc, file.path(opt$outdir, "mag_diagnosis_association.tsv"))
summ <- data.frame(n_mags = nrow(mags), n_high_quality = sum(mags$quality == "HQ"),
                   n_novel_species = sum(mags$novel_species), n_not_in_metaphlan = sum(!mags$in_metaphlan_profile),
                   n_with_amr = sum(mags$n_amr_genes > 0),
                   n_disease_enriched = length(unique(assoc$mag[assoc$qvalue < st$q_threshold & assoc$estimate > 0])))
write_tsv(summ, file.path(opt$outdir, "mag_summary_stats.tsv"))

# ---- 3. Figure 6: tree + annotation rings
tree <- keep.tip(tree, intersect(tree$tip.label, mags$mag))
tree <- ladderize(tree)
tips <- tree$tip.label
ann <- out[match(tips, out$mag), ]
best <- do.call(rbind, lapply(split(assoc, assoc$mag), function(x) x[which.min(x$pvalue), ]))
ann$enrich <- "ns"
b <- best[match(tips, best$mag), ]
ann$enrich[!is.na(b$qvalue) & b$qvalue < st$q_threshold] <- paste0(ifelse(b$estimate[!is.na(b$qvalue) & b$qvalue < st$q_threshold] > 0, "up in ", "down in "),
                                         sub("_vs_.*", "", b$contrast[!is.na(b$qvalue) & b$qvalue < st$q_threshold]))
phy_col <- c(Firmicutes = "#88CCEE", Bacteroidota = "#CC6677", Actinobacteria = "#DDCC77", Proteobacteria = "#117733",
             Verrucomicrobiota = "#AA4499", Other = "grey60")
phy <- sub(";c__.*", "", sub(".*;p__", "", ann$gtdb_classification)); phy <- ifelse(phy %in% names(phy_col), phy, "Other")
en_col <- c(ns = "grey88", "up in CD" = "#D55E00", "up in UC" = "#009E73", "down in CD" = "#E69F00", "down in UC" = "#56B4E9")
fn <- function() {
  op <- par(mar = c(1, 1, 2, 1), xpd = NA); on.exit(par(op))
  n <- length(tips)
  layout(matrix(1:2, 1), widths = c(3, 4))
  plot(tree, type = "phylogram", cex = 0.6, label.offset = 0.01, y.lim = c(1, n), main = "MAG phylogeny (labels coloured by phylum)",
       tip.color = phy_col[phy])
  legend("topleft", legend = names(phy_col), text.col = phy_col, cex = 0.55, bty = "n")
  plot.new(); plot.window(xlim = c(0, 7), ylim = c(0.5, n + 0.5))
  rect(0.0, seq_len(n) - 0.4, 0.8, seq_len(n) + 0.4, col = en_col[ann$enrich], border = NA)
  rect(1.0, seq_len(n) - 0.4, 1.8, seq_len(n) + 0.4, col = ifelse(ann$novel_species, "#332288", "grey92"), border = NA)
  rect(2.0, seq_len(n) - 0.4, 2.8, seq_len(n) + 0.4, col = ifelse(ann$in_metaphlan_profile, "grey92", "#CC3311"), border = NA)
  rect(3.0, seq_len(n) - 0.4, 3.0 + 3.5 * ann$n_amr_genes / max(1, max(ann$n_amr_genes)), seq_len(n) + 0.4, col = "#117733", border = NA)
  text(c(0.4, 1.4, 2.4), n + 1.2, c("Disease\nassoc.", "Novel\nspecies", "Not in\nMetaPhlAn"), cex = 0.6, adj = c(0.5, 0))
  text(3.0, n + 1.2, "AMR genes\n(count)", cex = 0.6, adj = c(0, 0))
  text(3.0 + 3.5 * ann$n_amr_genes / max(1, max(ann$n_amr_genes)) + 0.1, seq_len(n), ifelse(ann$n_amr_genes > 0, ann$n_amr_genes, ""), cex = 0.5, adj = 0)
  legend("bottomright", legend = names(en_col), fill = en_col, cex = 0.5, bty = "n", title = "Association", title.adj = 0)
}
for (ext in c("pdf", "png")) {
  f <- file.path(opt$figdir, paste0("Fig6_MAG_tree.", ext)); dir.create(opt$figdir, showWarnings = FALSE, recursive = TRUE)
  if (ext == "pdf") pdf(f, width = 10, height = max(6, 0.22 * length(tips) + 1.5)) else png(f, width = 10, height = max(6, 0.22 * length(tips) + 1.5), units = "in", res = 200)
  fn(); dev.off()
}
log_msg("mag_analysis done")
