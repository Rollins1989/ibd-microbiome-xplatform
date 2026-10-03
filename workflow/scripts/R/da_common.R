# Shared setup for the differential-abundance scripts (sourced, not executed).
da_options <- function() list(
  optparse::make_option("--config"), optparse::make_option("--meta"), optparse::make_option("--layer"),
  optparse::make_option("--layer-name"), optparse::make_option("--type"), optparse::make_option("--out"),
  optparse::make_option("--threads", type = "integer", default = 1))

da_setup <- function(opt, baseline_only = FALSE) {
  cfg <- load_cfg(opt$config)
  meta <- read_meta(opt$meta)
  M <- read_layer(opt$layer)
  ids <- intersect(colnames(M), meta$sample_id)
  if (baseline_only) ids <- intersect(ids, meta$sample_id[meta$is_baseline %in% 1])
  md <- model_data(meta, cfg$design, ids)
  M <- M[, rownames(md$data), drop = FALSE]
  M <- filter_layer(M, cfg$stats$min_prevalence, cfg$stats$min_abundance_prop)
  M <- M[rowSums(M) > 0, , drop = FALSE]
  # Safe identifiers: tools mangle names such as "[Eubacterium] rectale group" or "PWY-1: foo".
  key <- data.frame(id = sprintf("F%05d", seq_len(nrow(M))), feature = rownames(M), stringsAsFactors = FALSE)
  rownames(M) <- key$id
  message(sprintf("[DA] layer=%s samples=%d subjects=%d features=%d covariates=[%s]", opt$`layer-name`,
                  ncol(M), nlevels(droplevels(md$data$subject)), nrow(M), paste(md$covars, collapse = ",")))
  list(cfg = cfg, md = md, M = M, key = key, type = opt$type,
       contrasts = levels(md$data$diagnosis)[-1], ref = levels(md$data$diagnosis)[1])
}

#' Standardised result table; BH q-values within contrast.
da_finish <- function(res, S, method, opt) {
  res <- as.data.frame(res)
  res$feature <- S$key$feature[match(res$id, S$key$id)]
  res$id <- NULL
  res$qvalue <- ave(res$pvalue, res$contrast, FUN = function(p) p.adjust(p, "BH"))
  res$method <- method; res$layer <- opt$`layer-name`
  res$prevalence <- rowMeans(S$M[S$key$id[match(res$feature, S$key$feature)], , drop = FALSE] > 0)
  res$n_samples <- ncol(S$M); res$n_subjects <- nlevels(droplevels(S$md$data$subject))
  cols <- c("layer", "method", "feature", "contrast", "estimate", "se", "pvalue", "qvalue", "prevalence", "n_samples", "n_subjects")
  write_tsv(res[, cols], opt$out)
  message("[DA] ", method, ": ", sum(res$qvalue < S$cfg$stats$q_threshold, na.rm = TRUE), " hits at q<", S$cfg$stats$q_threshold)
}
