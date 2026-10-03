#!/usr/bin/env Rscript
# ANCOM-BC2 (counts only) with a subject random intercept.
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])), "utils.R"))
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])), "da_common.R"))
opt <- cli_args(da_options())
suppressPackageStartupMessages({library(ANCOMBC); library(phyloseq)})
S <- da_setup(opt)
stopifnot(S$type == "counts")
d <- S$md$data
md <- d[, c("diagnosis", S$md$covars, "subject"), drop = FALSE]
ps <- phyloseq(otu_table(round(S$M), taxa_are_rows = TRUE), sample_data(md))
use_re <- anyDuplicated(d$subject) > 0
out <- ancombc2(data = ps, fix_formula = rhs(c("diagnosis", S$md$covars)),
                rand_formula = if (use_re) "(1 | subject)" else NULL, p_adj_method = "BH",
                prv_cut = 0, lib_cut = 0, group = "diagnosis", struc_zero = FALSE, neg_lb = FALSE,
                alpha = 0.1, n_cl = opt$threads, verbose = FALSE)
r <- out$res
res <- do.call(rbind, lapply(S$contrasts, function(l)
  data.frame(id = r$taxon, contrast = paste0(l, "_vs_", S$ref), estimate = r[[paste0("lfc_diagnosis", l)]],
             se = r[[paste0("se_diagnosis", l)]], pvalue = r[[paste0("p_diagnosis", l)]])))
da_finish(res, S, "ancombc2", opt)
