#!/usr/bin/env Rscript
# Cross-platform concordance on PAIRED specimens (same stool sample sequenced by 16S and shotgun):
#   * genus-level: Procrustes/PROTEST + Mantel on Aitchison distances, per-genus Spearman,
#     Bland-Altman for Shannon diversity
#   * function: PICRUSt2-predicted vs HUMAnN-measured MetaCyc pathways (per-pathway Spearman, Mantel)
# Genus labels are harmonised (SILVA/GTDB/MetaPhlAn naming) with norm_taxon(); the share of each
# platform's total abundance that maps to the shared genus set is reported so name mismatches are visible.
suppressPackageStartupMessages({library(optparse); library(vegan); library(patchwork)})
.sd <- dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])); source(file.path(.sd, "utils.R"))
opt <- cli_args(list(make_option("--config"), make_option("--specimens"), make_option("--meta"),
                     make_option("--layer-dir"), make_option("--outdir"), make_option("--figdir")))
cfg <- load_cfg(opt$config); set.seed(cfg$stats$seed)
sp <- read_meta(opt$specimens); dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)
L <- function(n) read_layer(file.path(opt$`layer-dir`, paste0(n, ".tsv")))
g16 <- L("16s_genus"); gmg <- L("mgx_genus"); p16 <- L("16s_pathway"); pmg <- if (file.exists(file.path(opt$`layer-dir`, "mgx_pathway.tsv"))) L("mgx_pathway") else matrix(0, 0, 0)

pairs <- sp[!is.na(sp$sample_16s) & !is.na(sp$sample_mgx) & sp$sample_16s %in% colnames(g16) & sp$sample_mgx %in% colnames(gmg), ]
stopifnot(nrow(pairs) >= 10)
message("[concordance] paired specimens: ", nrow(pairs))

# ---- genus harmonisation
k16 <- norm_taxon(rownames(g16)); kmg <- norm_taxon(rownames(gmg))
H16 <- collapse_by_key(g16[, pairs$sample_16s], k16); Hmg <- collapse_by_key(gmg[, pairs$sample_mgx], kmg)
colnames(H16) <- colnames(Hmg) <- pairs$specimen_id
shared <- intersect(rownames(H16), rownames(Hmg))
cov <- data.frame(platform = c("16S", "MGX"),
  n_genus_total = c(sum(!is.na(k16)), sum(!is.na(kmg))), n_genus_shared = length(shared),
  pct_abundance_in_shared = 100 * c(sum(H16[shared, ]) / sum(g16[, pairs$sample_16s]), sum(Hmg[shared, ]) / sum(gmg[, pairs$sample_mgx])))
write_tsv(cov, file.path(opt$outdir, "genus_harmonisation_coverage.tsv"))
H16 <- H16[shared, ]; Hmg <- Hmg[shared, ]
keep <- rowMeans(H16 > 0) >= 0.1 | rowMeans(Hmg > 0) >= 0.1
H16 <- H16[keep, ]; Hmg <- Hmg[keep, ]

# ---- Procrustes / Mantel on Aitchison distances (shared genera)
a16 <- dist(t(clr_mat(H16, "counts"))); amg <- dist(t(clr_mat(Hmg, "counts")))
k <- min(5, nrow(H16) - 1)
X <- cmdscale(a16, k = k); Y <- cmdscale(amg, k = k)
pt <- protest(X, Y, permutations = 999, symmetric = TRUE)
mt <- mantel(a16, amg, method = "spearman", permutations = 999)
write_tsv(data.frame(test = c("PROTEST", "Mantel (Spearman)"), statistic_name = c("symmetric Procrustes correlation", "Mantel r"),
                     statistic = c(sqrt(1 - pt$ss), mt$statistic), p = c(pt$signif, mt$signif), n_pairs = nrow(pairs), n_genera = nrow(H16)),
          file.path(opt$outdir, "genus_beta_concordance.tsv"))

# ---- per-genus Spearman on relative abundance
P16 <- sweep(H16, 2, colSums(H16), "/"); Pmg <- sweep(Hmg, 2, colSums(Hmg), "/")
gs <- do.call(rbind, lapply(rownames(P16), function(g) {
  ct <- suppressWarnings(cor.test(P16[g, ], Pmg[g, ], method = "spearman", exact = FALSE))
  data.frame(genus = g, rho = unname(ct$estimate), p = ct$p.value, mean_prop_16s = mean(P16[g, ]), mean_prop_mgx = mean(Pmg[g, ]),
             prev_16s = mean(P16[g, ] > 0), prev_mgx = mean(Pmg[g, ] > 0))
}))
gs$q <- p.adjust(gs$p, "BH"); gs <- gs[order(-gs$rho), ]
write_tsv(gs, file.path(opt$outdir, "genus_spearman.tsv"))

# ---- Bland-Altman on Shannon (shared genera) and on log10 relative abundance
sh16 <- vegan::diversity(t(P16)); shmg <- vegan::diversity(t(Pmg))
ba <- data.frame(mean = (sh16 + shmg) / 2, diff = sh16 - shmg)
bias <- mean(ba$diff); loa <- bias + c(-1.96, 1.96) * sd(ba$diff)
ct_sh <- cor.test(sh16, shmg, method = "spearman", exact = FALSE)
write_tsv(data.frame(metric = "Shannon (shared genera)", bias_16S_minus_MGX = bias, loa_lower = loa[1], loa_upper = loa[2],
                     spearman = unname(ct_sh$estimate), p = ct_sh$p.value, n = nrow(ba)), file.path(opt$outdir, "bland_altman_shannon.tsv"))
eps <- 1e-5
lab <- data.frame(l16 = log10(as.vector(P16) + eps), lmg = log10(as.vector(Pmg) + eps))
lab$m <- (lab$l16 + lab$lmg) / 2; lab$d <- lab$l16 - lab$lmg

# ---- pathways
ps <- intersect(rownames(p16), rownames(pmg))
if (length(ps) >= 5) {
  pr <- pairs[pairs$sample_16s %in% colnames(p16) & pairs$sample_mgx %in% colnames(pmg), ]
  Q16 <- p16[ps, pr$sample_16s]; Qmg <- pmg[ps, pr$sample_mgx]
  Q16 <- sweep(Q16, 2, colSums(Q16), "/"); Qmg <- sweep(Qmg, 2, colSums(Qmg), "/")
  ps_df <- do.call(rbind, lapply(ps, function(p) {
    ct <- suppressWarnings(cor.test(Q16[p, ], Qmg[p, ], method = "spearman", exact = FALSE))
    data.frame(pathway = p, rho = unname(ct$estimate), p = ct$p.value)
  }))
  ps_df$q <- p.adjust(ps_df$p, "BH"); ps_df <- ps_df[order(-ps_df$rho), ]
  write_tsv(ps_df, file.path(opt$outdir, "pathway_spearman.tsv"))
  mp <- mantel(dist(t(clr_mat(Q16, "abundance"))), dist(t(clr_mat(Qmg, "abundance"))), method = "spearman", permutations = 999)
  write_tsv(data.frame(test = "Mantel (Spearman), pathway Aitchison", r = mp$statistic, p = mp$signif, n_pairs = nrow(pr),
                       n_pathways = length(ps), median_pathway_rho = median(ps_df$rho)),
            file.path(opt$outdir, "pathway_beta_concordance.tsv"))
}

# ---- Figure 3
pg <- merge(data.frame(sp[, c("specimen_id", "diagnosis")]), data.frame(specimen_id = pairs$specimen_id), by = "specimen_id")
rot <- pt$Yrot; cx <- pt$X
pdat <- data.frame(x1 = cx[, 1], y1 = cx[, 2], x2 = rot[, 1], y2 = rot[, 2], specimen_id = rownames(cx))
pdat$diagnosis <- factor(sp$diagnosis[match(pdat$specimen_id, sp$specimen_id)], levels = cfg$design$diagnosis_levels)
pA <- ggplot(pdat) + geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2, colour = diagnosis), linewidth = 0.25, alpha = 0.6) +
  geom_point(aes(x1, y1, colour = diagnosis), shape = 16, size = 1) + geom_point(aes(x2, y2, colour = diagnosis), shape = 17, size = 1) +
  scale_colour_manual(values = DX_COLORS) +
  labs(title = "Procrustes: 16S (circle) vs shotgun (triangle)", x = "Dim 1", y = "Dim 2",
       subtitle = sprintf("PROTEST r = %.2f, p = %s; n = %d specimens", sqrt(1 - pt$ss), format.pval(pt$signif, digits = 2), nrow(pairs))) + theme_ibd()
top <- head(gs[order(-(gs$mean_prop_16s + gs$mean_prop_mgx)), ], 20)
top$genus <- factor(top$genus, levels = rev(top$genus[order(top$rho)]))
pB <- ggplot(top, aes(rho, genus, fill = rho)) + geom_col(width = 0.7) + scale_fill_gradient(low = "#9ecae1", high = "#08519c", guide = "none") +
  labs(title = "Per-genus agreement (top 20 by abundance)", x = "Spearman rho (16S vs shotgun)", y = NULL) + theme_ibd()
pC <- ggplot(ba, aes(mean, diff)) + geom_point(size = 0.8, alpha = 0.5) + geom_hline(yintercept = bias, colour = "#D55E00") +
  geom_hline(yintercept = loa, linetype = 2, colour = "grey40") +
  labs(title = "Bland-Altman: Shannon (shared genera)", x = "Mean of platforms", y = "16S - shotgun",
       subtitle = sprintf("bias %.2f, LoA [%.2f, %.2f]", bias, loa[1], loa[2])) + theme_ibd()
pD <- ggplot(lab, aes(m, d)) + geom_bin2d(bins = 60) + scale_fill_viridis_c(trans = "log10", name = "n") +
  geom_hline(yintercept = 0, colour = "white", linewidth = 0.3) +
  labs(title = "Bland-Altman: log10 genus abundance", x = "Mean log10 proportion", y = "16S - shotgun") + theme_ibd()
plots <- list(pA, pB, pC, pD)
if (exists("ps_df")) {
  pE <- ggplot(ps_df, aes(rho)) + geom_histogram(bins = 25, fill = "grey40", colour = "white") +
    geom_vline(xintercept = median(ps_df$rho), colour = "#D55E00") +
    labs(title = "PICRUSt2 vs HUMAnN pathways", x = "Spearman rho per pathway", y = "Pathways",
         subtitle = sprintf("median rho = %.2f; Mantel r = %.2f", median(ps_df$rho), mp$statistic)) + theme_ibd()
  plots[[5]] <- pE
}
fig <- wrap_plots(plots, ncol = 2) + plot_annotation(tag_levels = "A")
save_fig(fig, file.path(opt$figdir, "Fig3_platform_concordance"), w = 10, h = 12)
message("[concordance] done")
