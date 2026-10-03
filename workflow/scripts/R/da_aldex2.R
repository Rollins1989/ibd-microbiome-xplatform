#!/usr/bin/env Rscript
# ALDEx2 (Welch t on CLR Monte-Carlo instances; counts only). ALDEx2 supports neither random effects
# nor covariates, so it runs on BASELINE samples (1 per subject => independent), one contrast at a time.
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])), "utils.R"))
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])), "da_common.R"))
opt <- cli_args(da_options())
suppressPackageStartupMessages(library(ALDEx2))
S <- da_setup(opt, baseline_only = TRUE)
stopifnot(S$type == "counts")
d <- S$md$data
clr <- clr_mat(S$M, "counts")
res <- list()
for (l in S$contrasts) {
  keep <- d$diagnosis %in% c(S$ref, l)
  cond <- as.character(d$diagnosis[keep])
  ax <- aldex(round(S$M[, keep, drop = FALSE]), cond, mc.samples = 128, test = "t", effect = FALSE,
              denom = "all", verbose = FALSE)
  # sign from our own CLR medians (test minus reference): unambiguous, comparable across methods
  est <- apply(clr[, keep, drop = FALSE], 1, function(x) median(x[cond == l]) - median(x[cond == S$ref]))
  res[[l]] <- data.frame(id = rownames(ax), contrast = paste0(l, "_vs_", S$ref), estimate = est[rownames(ax)],
                         se = NA_real_, pvalue = ax$we.ep)
}
da_finish(do.call(rbind, res), S, "aldex2", opt)
