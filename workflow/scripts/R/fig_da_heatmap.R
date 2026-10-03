#!/usr/bin/env Rscript
# Figure 4: consensus differential-abundance heatmap across layers and contrasts.
suppressPackageStartupMessages(library(optparse))
.sd <- dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1]))
source(file.path(.sd, "utils.R"))
opt <- cli_args(list(make_option("--config"), make_option("--inputs"), make_option("--figdir"),
                     make_option("--top", type = "integer", default = 20)))
cfg <- load_cfg(opt$config)
files <- strsplit(opt$inputs, ",")[[1]]
d <- rbindlist(lapply(files, fread), fill = TRUE)
d[, consensus := as.logical(consensus)]
hits <- d[consensus == TRUE]
if (nrow(hits) == 0) {
  p <- ggplot() + annotate("text", x = 0, y = 0, label = "No consensus features") + theme_void()
  save_fig(p, file.path(opt$figdir, "Fig4_differential_abundance"), 5, 3); quit(save = "no")
}
keep <- hits[, .(m = max(abs(effect))), by = .(layer, feature)][order(-m)][, head(.SD, opt$top), by = layer]
plt <- merge(d, keep[, .(layer, feature)], by = c("layer", "feature"))
plt[, label := ifelse(consensus, "*", "")]
ord <- plt[contrast == plt$contrast[1], .(layer, feature, e = effect)]
ord <- ord[order(layer, -e)]
plt[, feature := factor(feature, levels = rev(unique(ord$feature)))]
lim <- max(abs(plt$effect), na.rm = TRUE)
plt[, layer_label := gsub("_", " ", layer)]
p <- ggplot(plt, aes(contrast, feature, fill = effect)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = label), size = 4, vjust = 0.75) +
  facet_grid(layer_label ~ ., scales = "free_y", space = "free_y") +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", limits = c(-lim, lim),
                       name = "CLR\neffect", na.value = "grey85") +
  labs(title = "Differentially abundant features (consensus)", x = NULL, y = NULL,
       caption = sprintf("* significant (q < %.2f) in >= %d methods with concordant direction", cfg$stats$q_threshold, cfg$stats$min_methods)) +
  theme_ibd() + theme(strip.text.y = element_text(angle = 0), panel.grid = element_blank())
save_fig(p, file.path(opt$figdir, "Fig4_differential_abundance"), w = 6.5, h = max(4, 0.17 * nrow(keep) + 2))
