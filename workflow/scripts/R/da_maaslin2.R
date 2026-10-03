#!/usr/bin/env Rscript
# MaAsLin2: CLR-normalised linear mixed model, subject random effect.
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])), "utils.R"))
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])), "da_common.R"))
opt <- cli_args(da_options())
suppressPackageStartupMessages(library(Maaslin2))
S <- da_setup(opt)
d <- S$md$data
feat <- as.data.frame(t(S$M))                      # samples x features (safe ids)
md <- d[, c("diagnosis", S$md$covars, "subject"), drop = FALSE]
use_re <- anyDuplicated(d$subject) > 0
ref <- paste0("diagnosis,", S$ref)
for (v in S$md$covars) if (is.factor(d[[v]]) && nlevels(d[[v]]) > 2) ref <- c(ref, paste0(v, ",", levels(d[[v]])[1]))
wd <- file.path(tempdir(), "maaslin2")
fit <- Maaslin2(input_data = feat, input_metadata = md, output = wd, min_abundance = 0, min_prevalence = 0,
                normalization = "CLR", transform = "NONE", analysis_method = "LM", max_significance = 0.1,
                fixed_effects = c("diagnosis", S$md$covars), random_effects = if (use_re) "subject" else NULL,
                reference = ref, standardize = FALSE, cores = opt$threads, plot_heatmap = FALSE, plot_scatter = FALSE)
r <- fit$results[fit$results$metadata == "diagnosis", ]
res <- data.frame(id = r$feature, contrast = paste0(r$value, "_vs_", S$ref), estimate = r$coef, se = r$stderr, pvalue = r$pval)
da_finish(res, S, "maaslin2", opt)
