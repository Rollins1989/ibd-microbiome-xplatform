#!/usr/bin/env Rscript
# Alpha / beta diversity, PERMANOVA (independent samples), dispersion, rarefaction sensitivity, Fig 2.
#
# Statistical notes
#  * Alpha diversity: linear mixed model with random subject intercept (all samples).
#  * PERMANOVA assumes exchangeable observations. Repeated visits violate that, and restricting
#    permutations *within* subject cannot test a between-subject factor. We therefore (i) test
#    one sample per subject (baseline) as the primary analysis and (ii) repeat the test on many
#    random one-sample-per-subject draws as a sensitivity analysis.
suppressPackageStartupMessages({library(vegan); library(optparse); library(patchwork)})
.sd <- dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])); source(file.path(.sd, "utils.R"))
opt <- cli_args(list(
  make_option("--config"), make_option("--meta"), make_option("--layer-dir"),
  make_option("--tree"), make_option("--outdir"), make_option("--figdir")))
cfg <- load_cfg(opt$config); design <- cfg$design; st <- cfg$stats
set.seed(st$seed)
meta <- read_meta(opt$meta)
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

faith_pd <- function(counts, tree) {
  tree <- ape::keep.tip(tree, intersect(tree$tip.label, rownames(counts)))
  counts <- counts[tree$tip.label, , drop = FALSE]
  n <- length(tree$tip.label)
  pp <- ape::prop.part(tree)           # tip sets for internal nodes (ids n+2 ...)
  E <- Matrix::Matrix(0, nrow(tree$edge), n, sparse = TRUE)
  for (i in seq_len(nrow(tree$edge))) {
    ch <- tree$edge[i, 2]
    tips <- if (ch <= n) ch else pp[[ch - n - 1]]
    E[i, tips] <- 1
  }
  P <- Matrix::Matrix((counts > 0) * 1, sparse = TRUE)
  M <- E %*% P
  as.numeric(Matrix::colSums((M > 0) * tree$edge.length))
}

platforms <- list(`16S` = list(layer = "asv_table.filtered.tsv", label = "16S rRNA"),
                  MGX   = list(layer = "mgx_species.tsv", label = "Shotgun"))
alpha_all <- list(); pcoa_plots <- list(); r2_rows <- list()
perm_all <- list(); pair_all <- list(); disp_all <- list(); sub_all <- list(); rar_all <- list(); amod_all <- list()

for (P in names(platforms)) {
  C <- read_layer(file.path(opt$`layer-dir`, platforms[[P]]$layer))
  m <- meta[meta$data_type == P & meta$sample_id %in% colnames(C), ]
  C <- C[, m$sample_id, drop = FALSE]; C <- C[rowSums(C) > 0, , drop = FALSE]
  message("[", P, "] ", ncol(C), " samples, ", nrow(C), " features")
  prop <- sweep(C, 2, colSums(C), "/")

  # ---- alpha
  al <- data.frame(sample_id = colnames(C), platform = P,
                   observed = colSums(C > 0),
                   shannon = vegan::diversity(t(C), "shannon"))
  if (P == "16S") {
    tree <- ape::multi2di(ape::read.tree(opt$tree), random = FALSE)  # UniFrac needs a strictly binary tree
    al$faith_pd <- faith_pd(C, tree)
    if (isTRUE(st$rarefaction_sensitivity)) {
      d <- min(colSums(C)); rr <- t(vegan::rrarefy(t(C), d))
      rar_all[[P]] <- data.frame(platform = P, depth = d,
        metric = c("shannon", "observed"),
        spearman_vs_unrarefied = c(cor(vegan::diversity(t(rr), "shannon"), al$shannon, method = "spearman"),
                                   cor(colSums(rr > 0), al$observed, method = "spearman")))
    }
  }
  al <- merge(al, m[, c("sample_id", "subject_id", "diagnosis", "week", "is_baseline")], by = "sample_id")
  alpha_all[[P]] <- al
  md <- model_data(m, design, al$sample_id)
  rownames(al) <- al$sample_id
  for (met in intersect(c("observed", "shannon", "faith_pd"), colnames(al))) {
    fit <- fit_lmm(al[rownames(md$data), met], md$data, c("diagnosis", md$covars))
    ct <- coef_table(fit); ct <- ct[grepl("^diagnosis", ct$term), ]
    ct$metric <- met; ct$platform <- P; ct$n_samples <- nrow(md$data)
    ct$n_subjects <- nlevels(droplevels(md$data$subject)); amod_all[[paste(P, met)]] <- ct
  }

  # ---- distances (all samples, for ordination)
  clr <- clr_mat(C, "counts")
  dists <- list(aitchison = dist(t(clr)), bray = vegdist(t(prop), "bray"))
  if (P == "16S") {
    ps <- phyloseq::phyloseq(phyloseq::otu_table(C, taxa_are_rows = TRUE),
                             phyloseq::phy_tree(ape::keep.tip(tree, rownames(C))))
    dists$unifrac_unweighted <- phyloseq::UniFrac(ps, weighted = FALSE)
    dists$unifrac_weighted <- phyloseq::UniFrac(ps, weighted = TRUE)
  }
  base_ids <- m$sample_id[m$is_baseline %in% 1]
  mdb <- model_data(m, design, base_ids)
  bd <- mdb$data
  for (dn in names(dists)) {
    D <- as.matrix(dists[[dn]])
    # PCoA plot
    pc <- cmdscale(as.dist(D), k = 2, eig = TRUE)
    ve <- 100 * pc$eig[1:2] / sum(pc$eig[pc$eig > 0])
    pdat <- data.frame(PC1 = pc$points[, 1], PC2 = pc$points[, 2], sample_id = rownames(pc$points))
    pdat$diagnosis <- m$diagnosis[match(pdat$sample_id, m$sample_id)]
    pdat$diagnosis <- factor(pdat$diagnosis, levels = design$diagnosis_levels)
    pcoa_plots[[paste(P, dn)]] <- ggplot(pdat, aes(PC1, PC2, colour = diagnosis, fill = diagnosis)) +
      geom_point(size = 1, alpha = 0.55) + stat_ellipse(geom = "polygon", alpha = 0.08, colour = NA, level = 0.8) +
      scale_colour_manual(values = DX_COLORS) + scale_fill_manual(values = DX_COLORS) +
      labs(title = paste(platforms[[P]]$label, "-", dn), x = sprintf("PCo1 (%.1f%%)", ve[1]), y = sprintf("PCo2 (%.1f%%)", ve[2])) +
      theme_ibd()
    # PERMANOVA on baseline samples
    Db <- as.dist(D[rownames(bd), rownames(bd)])
    f <- as.formula(paste("Db ~", rhs(c(mdb$covars, "diagnosis"))))
    ad <- adonis2(f, data = bd, by = "margin", permutations = st$permanova$permutations)
    adt <- data.frame(term = rownames(ad), ad, row.names = NULL, check.names = FALSE)
    adt$platform <- P; adt$distance <- dn; adt$n <- nrow(bd)
    perm_all[[paste(P, dn)]] <- adt
    # pairwise
    lv <- levels(bd$diagnosis); pw <- list()
    for (i in seq_along(lv)) for (j in seq_along(lv)) if (i < j) {
      sub <- bd[bd$diagnosis %in% c(lv[i], lv[j]), ]; sub$diagnosis <- droplevels(sub$diagnosis)
      Ds <- as.dist(D[rownames(sub), rownames(sub)])
      a2 <- adonis2(as.formula(paste("Ds ~", rhs(c(mdb$covars, "diagnosis")))), data = sub, by = "margin",
                    permutations = st$permanova$permutations)
      pw[[paste(i, j)]] <- data.frame(platform = P, distance = dn, contrast = paste(lv[j], "vs", lv[i]),
                                      R2 = a2["diagnosis", "R2"], p = a2["diagnosis", "Pr(>F)"], n = nrow(sub))
    }
    pw <- do.call(rbind, pw); pw$p_adj <- p.adjust(pw$p, "BH"); pair_all[[paste(P, dn)]] <- pw
    # dispersion
    bdisp <- betadisper(Db, bd$diagnosis)
    pt <- permutest(bdisp, permutations = st$permanova$permutations)
    disp_all[[paste(P, dn)]] <- data.frame(platform = P, distance = dn, F = pt$tab[1, "F"], p = pt$tab[1, "Pr(>F)"])
    # subsampling sensitivity: random 1 sample / subject
    sm <- split(m$sample_id, m$subject_id); res <- NULL
    for (k in seq_len(st$permanova$subsample_draws)) {
      ids <- vapply(sm, function(x) x[sample.int(length(x), 1)], character(1))
      dd <- model_data(m, design, ids)$data
      Dk <- as.dist(D[rownames(dd), rownames(dd)])
      cv <- model_data(m, design, ids)$covars
      a3 <- adonis2(as.formula(paste("Dk ~", rhs(c(cv, "diagnosis")))), data = dd, by = "margin", permutations = 199)
      res <- rbind(res, data.frame(R2 = a3["diagnosis", "R2"], p = a3["diagnosis", "Pr(>F)"]))
    }
    sub_all[[paste(P, dn)]] <- data.frame(platform = P, distance = dn, draws = nrow(res), median_R2 = median(res$R2),
                                          median_p = median(res$p), frac_p_lt_0.05 = mean(res$p < 0.05))
  }
}

alpha <- rbindlist(alpha_all, fill = TRUE)
write_tsv(alpha, file.path(opt$outdir, "alpha_diversity.tsv"))
write_tsv(rbindlist(amod_all), file.path(opt$outdir, "alpha_models.tsv"))
write_tsv(rbindlist(perm_all, fill = TRUE), file.path(opt$outdir, "permanova_baseline.tsv"))
write_tsv(rbindlist(pair_all), file.path(opt$outdir, "permanova_pairwise.tsv"))
write_tsv(rbindlist(disp_all), file.path(opt$outdir, "betadisper.tsv"))
write_tsv(rbindlist(sub_all), file.path(opt$outdir, "permanova_subsample_sensitivity.tsv"))
if (length(rar_all)) write_tsv(rbindlist(rar_all), file.path(opt$outdir, "rarefaction_sensitivity.tsv"))

# ---------------- Figure 2
alpha$diagnosis <- factor(alpha$diagnosis, levels = design$diagnosis_levels)
alpha$platform <- factor(alpha$platform, levels = c("16S", "MGX"), labels = c("16S rRNA", "Shotgun"))
pA <- ggplot(alpha, aes(diagnosis, shannon, fill = diagnosis)) +
  geom_boxplot(outlier.size = 0.4, linewidth = 0.3, alpha = 0.85) + facet_wrap(~platform, scales = "free_y") +
  scale_fill_manual(values = DX_COLORS) + labs(title = "Shannon diversity", x = NULL, y = "Shannon index") +
  theme_ibd() + theme(legend.position = "none")
pr <- rbindlist(perm_all, fill = TRUE)[term == "diagnosis"]
pr$label <- paste(pr$platform, pr$distance)
pE <- ggplot(pr, aes(reorder(label, R2), R2)) + geom_col(fill = "grey45", width = 0.7) +
  geom_text(aes(label = ifelse(`Pr(>F)` < 0.001, "p<0.001", sprintf("p=%.3f", `Pr(>F)`))), hjust = -0.1, size = 2.6) +
  coord_flip(clip = "off") + scale_y_continuous(expand = expansion(mult = c(0, 0.35))) +
  labs(title = "PERMANOVA (baseline, 1 sample/subject)", x = NULL, y = expression(R^2)) + theme_ibd()
pp <- pcoa_plots
fig <- (pA | pE) / (pp[["16S aitchison"]] | pp[["MGX aitchison"]]) / (pp[["16S unifrac_weighted"]] | pp[["MGX bray"]]) +
  plot_layout(guides = "collect") + plot_annotation(tag_levels = "A")
save_fig(fig, file.path(opt$figdir, "Fig2_diversity"), w = 9, h = 11)
message("[diversity] done")
