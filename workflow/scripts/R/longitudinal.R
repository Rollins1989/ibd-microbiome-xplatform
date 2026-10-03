#!/usr/bin/env Rscript
# Longitudinal analyses (repeated measures -> mixed models, random subject intercept):
#   1. alpha diversity ~ diagnosis * time + covariates + (1|subject)
#   2. within-IBD association of each feature with disease activity (active vs inactive)
#   3. flare-aligned trajectories: visits relative to the first active visit after an inactive one,
#      with a paired within-event test of the pre-flare visit (-1) against earlier quiescent visits
suppressPackageStartupMessages({library(optparse); library(patchwork)})
.sd <- dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])); source(file.path(.sd, "utils.R"))
opt <- cli_args(list(make_option("--config"), make_option("--meta"), make_option("--specimens"), make_option("--events"),
                     make_option("--alpha"), make_option("--layer-dir"), make_option("--layers", default = "16s_genus,mgx_species"),
                     make_option("--outdir"), make_option("--figdir"), make_option("--top", type = "integer", default = 6)))
cfg <- load_cfg(opt$config); design <- cfg$design; set.seed(cfg$stats$seed)
meta <- read_meta(opt$meta); sp <- read_meta(opt$specimens); ev <- read_meta(opt$events)
alpha <- read_meta(opt$alpha); dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)
layers <- strsplit(opt$layers, ",")[[1]]
plat_of <- function(l) if (grepl("^16s", l)) "16S" else "MGX"

# ---- 1. alpha ~ diagnosis x time
rows <- list()
for (P in unique(alpha$platform)) {
  a <- alpha[alpha$platform == P, ]
  mm <- merge(a, meta[, setdiff(colnames(meta), c("subject_id", "diagnosis", "week", "is_baseline"))], by = "sample_id")
  mm$subject_id <- a$subject_id[match(mm$sample_id, a$sample_id)]; mm$diagnosis <- a$diagnosis[match(mm$sample_id, a$sample_id)]
  mm$week <- a$week[match(mm$sample_id, a$sample_id)]
  md <- model_data(mm, design); d <- md$data
  d$week_c <- as.numeric(scale(d$week, scale = FALSE))
  for (met in intersect(c("shannon", "observed", "faith_pd"), colnames(d))) {
    if (all(is.na(d[[met]])) || isTRUE(sd(d[[met]], na.rm = TRUE) == 0)) next  # undefined (Faith PD on shotgun) or constant
    d$.y <- as.numeric(scale(d[[met]]))       # standardised response -> comparable, better-conditioned fits
    f <- as.formula(paste(".y ~ diagnosis * week_c +", rhs(md$covars), "+ (1|subject)"))
    fit <- suppressMessages(lmerTest::lmer(f, data = d))
    an <- anova(fit)
    rows[[paste(P, met)]] <- data.frame(platform = P, metric = met, term = rownames(an), F = an[["F value"]], p = an[["Pr(>F)"]], n = nrow(d))
  }
}
write_tsv(do.call(rbind, rows), file.path(opt$outdir, "alpha_time_models.tsv"))

# ---- 2. activity association within IBD
act_rows <- list(); clr_store <- list()
for (l in layers) {
  M <- read_layer(file.path(opt$`layer-dir`, paste0(l, ".tsv")))
  m <- meta[meta$data_type == plat_of(l) & meta$sample_id %in% colnames(M) & meta$diagnosis != design$reference_level & !is.na(meta$active), ]
  md <- model_data(m, design); d <- md$data
  d$active <- factor(ifelse(d$active == 1, "active", "inactive"), levels = c("inactive", "active"))
  M <- filter_layer(M[, rownames(d), drop = FALSE], cfg$stats$min_prevalence, cfg$stats$min_abundance_prop)
  X <- clr_mat(M, "counts"); clr_store[[l]] <- clr_mat(read_layer(file.path(opt$`layer-dir`, paste0(l, ".tsv")))[rownames(M), , drop = FALSE], "counts")
  for (i in seq_len(nrow(X))) {
    fit <- try(fit_lmm(X[i, rownames(d)], d, c("active", "diagnosis", md$covars)), silent = TRUE)
    if (inherits(fit, "try-error")) next
    ct <- coef_table(fit); ct <- ct[ct$term == "activeactive", ]
    act_rows[[paste(l, i)]] <- data.frame(layer = l, feature = rownames(X)[i], estimate = ct$estimate, se = ct$se, p = ct$pvalue,
                                          n_samples = nrow(d), n_subjects = nlevels(droplevels(d$subject)))
  }
}
act <- do.call(rbind, act_rows); act$q <- ave(act$p, act$layer, FUN = function(p) p.adjust(p, "BH"))
act <- act[order(act$q), ]
write_tsv(act, file.path(opt$outdir, "activity_associations.tsv"))
message("[longitudinal] activity hits q<", cfg$stats$q_threshold, ": ", sum(act$q < cfg$stats$q_threshold))

# ---- 3. flare-aligned trajectories
sel <- head(act[order(act$q, -abs(act$estimate)), ], opt$top)
traj <- list(); shift <- list()
for (i in seq_len(nrow(sel))) {
  l <- sel$layer[i]; f <- sel$feature[i]; col <- paste0("sample_", tolower(plat_of(l)))
  e <- merge(ev, sp[, c("specimen_id", col)], by = "specimen_id"); e$sample_id <- e[[col]]
  e <- e[!is.na(e$sample_id) & e$sample_id %in% colnames(clr_store[[l]]), ]
  e$value <- clr_store[[l]][f, e$sample_id]; e$layer <- l; e$feature <- f
  traj[[i]] <- e[, c("layer", "feature", "event_id", "rel_visit", "value")]
  w <- tapply(e$value, list(e$event_id, ifelse(e$rel_visit == -1, "pre", ifelse(e$rel_visit < -1, "early", "other"))), mean)
  if (all(c("pre", "early") %in% colnames(w))) {
    ok <- complete.cases(w[, c("pre", "early")]); dd <- w[ok, "pre"] - w[ok, "early"]
    if (length(dd) >= 5) {
      wt <- suppressWarnings(wilcox.test(dd, exact = FALSE))
      shift[[i]] <- data.frame(layer = l, feature = f, n_events = length(dd), median_shift = median(dd), p = wt$p.value)
    }
  }
}
traj <- do.call(rbind, traj)
if (length(shift)) { shift <- do.call(rbind, shift); shift$q <- p.adjust(shift$p, "BH"); write_tsv(shift, file.path(opt$outdir, "preflare_shift.tsv")) }
write_tsv(traj, file.path(opt$outdir, "flare_trajectories_long.tsv"))

# ---- Figure 7
alpha$diagnosis <- factor(alpha$diagnosis, levels = design$diagnosis_levels)
alpha$platform <- factor(alpha$platform, levels = c("16S", "MGX"), labels = c("16S rRNA", "Shotgun"))
pA <- ggplot(alpha, aes(week, shannon, colour = diagnosis, fill = diagnosis)) +
  geom_point(size = 0.5, alpha = 0.25) + geom_smooth(method = "loess", formula = y ~ x, se = TRUE, alpha = 0.15, linewidth = 0.7) +
  facet_wrap(~platform, scales = "free_y") + scale_colour_manual(values = DX_COLORS) + scale_fill_manual(values = DX_COLORS) +
  labs(title = "Alpha diversity over time", x = "Week", y = "Shannon") + theme_ibd()
traj$label <- paste0(traj$feature, "\n(", gsub("_", " ", traj$layer), ")")
sm <- aggregate(value ~ label + rel_visit, traj, function(v) c(m = mean(v), se = sd(v) / sqrt(length(v))))
sm <- do.call(data.frame, sm); names(sm)[3:4] <- c("m", "se")
pB <- ggplot() + geom_line(data = traj, aes(rel_visit, value, group = event_id), colour = "grey70", linewidth = 0.2, alpha = 0.6) +
  geom_ribbon(data = sm, aes(rel_visit, ymin = m - se, ymax = m + se), fill = "#D55E00", alpha = 0.25) +
  geom_line(data = sm, aes(rel_visit, m), colour = "#D55E00", linewidth = 0.9) + geom_point(data = sm, aes(rel_visit, m), colour = "#D55E00", size = 1.2) +
  geom_vline(xintercept = -0.5, linetype = 2, colour = "grey30") + facet_wrap(~label, scales = "free_y", ncol = 3) +
  scale_x_continuous(breaks = seq(design$flare$window[1], design$flare$window[2])) +
  labs(title = "Flare-aligned trajectories (0 = first active visit)", x = "Visit relative to flare onset", y = "CLR abundance",
       caption = "Dashed line: transition from quiescence to flare. Grey = individual flare events; red = mean +/- SE.") + theme_ibd()
save_fig(pA / pB + plot_layout(heights = c(1, 1.8)) + plot_annotation(tag_levels = "A"), file.path(opt$figdir, "Fig7_longitudinal"), w = 9, h = 9)
message("[longitudinal] done")
