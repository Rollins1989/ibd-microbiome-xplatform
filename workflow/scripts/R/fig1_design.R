#!/usr/bin/env Rscript
# Figure 1: study design and sample overview.
suppressPackageStartupMessages(library(optparse))
.sd <- dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])); source(file.path(.sd, "utils.R"))
suppressPackageStartupMessages(library(patchwork))
opt <- cli_args(list(make_option("--config"), make_option("--specimens"), make_option("--figdir"), make_option("--outdir")))
cfg <- load_cfg(opt$config)
sp <- read_meta(opt$specimens)
sp$diagnosis <- factor(sp$diagnosis, levels = cfg$design$diagnosis_levels)
sp$platforms <- ifelse(!is.na(sp$sample_16s) & !is.na(sp$sample_mgx), "16S + shotgun",
                ifelse(!is.na(sp$sample_16s), "16S only", "Shotgun only"))
subj <- unique(sp[, c("subject_id", "diagnosis")])
tab <- as.data.frame(table(diagnosis = sp$diagnosis, platforms = sp$platforms))
write_tsv(data.frame(n_subjects = nrow(subj), n_specimens = nrow(sp),
                     n_paired = sum(sp$platforms == "16S + shotgun")), file.path(opt$outdir, "design_summary.tsv"))
write_tsv(tab, file.path(opt$outdir, "specimens_by_diagnosis_platform.tsv"))
pA <- ggplot(tab, aes(diagnosis, Freq, fill = platforms)) + geom_col(width = 0.7) +
  scale_fill_manual(values = c("16S + shotgun" = "#332288", "16S only" = "#88CCEE", "Shotgun only" = "#DDCC77")) +
  labs(title = "Specimens by platform", x = NULL, y = "Specimens", fill = NULL) + theme_ibd() + theme(legend.position = "bottom")
nvis <- as.data.frame(table(subject_id = sp$subject_id)); nvis$diagnosis <- subj$diagnosis[match(nvis$subject_id, subj$subject_id)]
pB <- ggplot(nvis, aes(Freq, fill = diagnosis)) + geom_histogram(binwidth = 1, colour = "white", linewidth = 0.2) +
  facet_wrap(~diagnosis, ncol = 1) + scale_fill_manual(values = DX_COLORS) +
  labs(title = "Visits per subject", x = "Visits", y = "Subjects") + theme_ibd() + theme(legend.position = "none")
ord <- sp[order(sp$diagnosis, -sp$week), ]
sp$subject_id <- factor(sp$subject_id, levels = unique(sp$subject_id[order(sp$diagnosis, sp$subject_id)]))
sp$state <- ifelse(is.na(sp$active), "no activity score", ifelse(sp$active == 1, "active", "inactive"))
pC <- ggplot(sp, aes(week, subject_id, colour = state, shape = platforms)) + geom_point(size = 1.2) +
  scale_colour_manual(values = c(active = "#D55E00", inactive = "#009E73", "no activity score" = "grey60")) +
  facet_grid(diagnosis ~ ., scales = "free_y", space = "free_y") +
  labs(title = "Sampling timeline", x = "Study week", y = NULL, colour = NULL, shape = NULL) +
  theme_ibd(8) + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.position = "bottom")
save_fig((pA / pB | pC) + plot_layout(widths = c(1, 1.6)) + plot_annotation(tag_levels = "A"),
         file.path(opt$figdir, "Fig1_study_design"), w = 10, h = 7)
