#!/usr/bin/env Rscript
# CLR + linear mixed model with a subject random intercept (Satterthwaite df; lmerTest).
# Valid for repeated measures: the between-subject diagnosis effect is tested with the
# appropriate degrees of freedom instead of treating visits as independent.
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])), "utils.R"))
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])), "da_common.R"))
opt <- cli_args(da_options())
suppressPackageStartupMessages(library(parallel))
S <- da_setup(opt)
clr <- clr_mat(S$M, S$type)
d <- S$md$data
terms <- c("diagnosis", S$md$covars)
one <- function(i) {
  ct <- tryCatch(coef_table(fit_lmm(clr[i, rownames(d)], d, terms)), error = function(e) NULL)
  if (is.null(ct)) return(NULL)
  ct <- ct[grepl("^diagnosis", ct$term), ]
  data.frame(id = rownames(clr)[i], contrast = paste0(sub("^diagnosis", "", ct$term), "_vs_", S$ref),
             estimate = ct$estimate, se = ct$se, pvalue = ct$pvalue)
}
res <- do.call(rbind, mclapply(seq_len(nrow(clr)), one, mc.cores = opt$threads))
da_finish(res, S, "clr_lmm", opt)
