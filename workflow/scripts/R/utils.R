# Shared helpers for all R analysis scripts (sourced, never executed directly).
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(yaml)
})

`%||%` <- function(a, b) if (is.null(a)) b else a
log_msg <- function(...) message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...)))

# ---- CLI / config -----------------------------------------------------------
cli_args <- function(option_list, description = "") {
  suppressPackageStartupMessages(library(optparse))
  parse_args(OptionParser(option_list = option_list, description = description),
             convert_hyphens_to_underscores = FALSE)
}
load_cfg <- function(path) yaml::read_yaml(path)

# ---- I/O --------------------------------------------------------------------
read_layer <- function(path) {
  dt <- fread(path, sep = "\t", data.table = FALSE, check.names = FALSE, na.strings = c("NA", ""))
  m <- as.matrix(dt[, -1, drop = FALSE])
  rownames(m) <- as.character(dt[[1]])
  storage.mode(m) <- "double"
  m
}
read_meta <- function(path) {
  d <- as.data.frame(fread(path, sep = "\t", data.table = FALSE, check.names = FALSE, na.strings = c("NA", "")))
  if ("is_baseline" %in% names(d)) d$is_baseline <- d$is_baseline %in% c(TRUE, "True", "TRUE", 1)
  d
}
write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  fwrite(as.data.frame(x), path, sep = "\t", na = "NA", quote = FALSE)
}

# ---- Compositional helpers ----------------------------------------------------
to_prop <- function(m) sweep(m, 2, pmax(colSums(m), 1e-300), "/")
clr_mat <- function(m, type = c("counts", "abundance")) {
  type <- match.arg(type)
  if (type == "counts") {
    x <- m + 0.5
  } else {
    p <- to_prop(m)
    x <- p + min(p[p > 0]) / 2
  }
  x <- sweep(x, 2, colSums(x), "/")
  lx <- log(x)
  sweep(lx, 2, colMeans(lx), "-")
}
filter_layer <- function(m, min_prev, min_prop = 0) {
  keep <- rowMeans(m > 0) >= min_prev
  if (!is.null(min_prop) && min_prop > 0) keep <- keep & rowMeans(to_prop(m)) >= min_prop
  m[keep, , drop = FALSE]
}

# ---- Modelling ------------------------------------------------------------------
rhs <- function(terms) paste(terms, collapse = " + ")

#' Analysis frame for the samples `ids`. `diagnosis` becomes a factor with the reference level first,
#' numeric covariates are z-scored, categorical ones become factors. Covariates that are too sparse,
#' constant or single-level are dropped (and reported), then incomplete rows are removed.
model_data <- function(meta, design, ids = NULL, require = character()) {
  d <- if (is.null(ids)) meta else meta[meta$sample_id %in% ids, , drop = FALSE]
  d$diagnosis <- factor(d$diagnosis, levels = design$diagnosis_levels)
  d <- d[!is.na(d$diagnosis), , drop = FALSE]
  covn <- intersect(unlist(design$covariates_numeric), names(d))
  covc <- intersect(unlist(design$covariates_categorical), names(d))
  if (isTRUE(design$include_batch) && "batch" %in% names(d)) covc <- c(covc, "batch")
  for (v in c(covn, covc)) {
    if (mean(is.na(d[[v]])) > design$max_covariate_missing) {
      log_msg("covariate '", v, "' dropped (>", 100 * design$max_covariate_missing, "% missing)")
      covn <- setdiff(covn, v); covc <- setdiff(covc, v)
    }
  }
  d <- d[complete.cases(d[, c("diagnosis", "subject_id", covn, covc, require), drop = FALSE]), , drop = FALSE]
  d$diagnosis <- droplevels(d$diagnosis)
  for (v in covn) {
    d[[v]] <- as.numeric(d[[v]])
    if (sd(d[[v]]) == 0) { log_msg("covariate '", v, "' dropped (constant)"); covn <- setdiff(covn, v) }
    else d[[v]] <- as.numeric(scale(d[[v]]))
  }
  for (v in covc) {
    d[[v]] <- droplevels(factor(d[[v]]))
    if (nlevels(d[[v]]) < 2) { log_msg("covariate '", v, "' dropped (single level)"); covc <- setdiff(covc, v) }
  }
  d$subject <- factor(d$subject_id)
  rownames(d) <- d$sample_id
  list(data = d, covars = c(covn, covc))
}

#' y ~ terms + (1 | subject) via lmerTest (Satterthwaite df); plain lm when no subject repeats.
fit_lmm <- function(y, data, terms) {
  suppressPackageStartupMessages(library(lmerTest))
  data$.y <- y
  if (anyDuplicated(data$subject) > 0) {
    f <- as.formula(paste(".y ~", rhs(terms), "+ (1 | subject)"))
    suppressWarnings(suppressMessages(
      lmerTest::lmer(f, data, REML = TRUE, control = lme4::lmerControl(check.conv.singular = "ignore"))))
  } else {
    lm(as.formula(paste(".y ~", rhs(terms))), data)
  }
}
coef_table <- function(fit) {
  co <- as.data.frame(summary(fit)$coefficients)
  data.frame(term = rownames(co), estimate = co[["Estimate"]], se = co[["Std. Error"]],
             statistic = co[["t value"]], pvalue = co[[grep("^Pr\\(", colnames(co))]], row.names = NULL)
}

# ---- Plot style -------------------------------------------------------------------
DX_COLORS <- c(nonIBD = "#0072B2", CD = "#D55E00", UC = "#009E73")   # Okabe-Ito, colour-blind safe
theme_ibd <- function(base = 10) {
  theme_bw(base_size = base) +
    theme(panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = "grey95", colour = NA),
          legend.key.size = unit(0.4, "cm"), plot.title = element_text(face = "bold", size = base + 1))
}
save_fig <- function(p, stem, w = 7, h = 5) {
  dir.create(dirname(stem), recursive = TRUE, showWarnings = FALSE)
  ggsave(paste0(stem, ".pdf"), p, width = w, height = h)
  ggsave(paste0(stem, ".png"), p, width = w, height = h, dpi = 300)
  invisible(stem)
}

# ---- Taxon-name harmonisation (mirror of mbstats.normalize_taxon_name) ----------------
norm_taxon <- function(x) {
  vapply(as.character(x), function(s) {
    if (is.na(s)) return(NA_character_)
    s <- trimws(s)
    s <- sub("^[a-z]__", "", s)
    s <- gsub("[][]", "", s)
    if (grepl("^(unclassified|uncultured|unknown|incertae|ggb[0-9]+|sgb[0-9]+)", s, ignore.case = TRUE) || !nzchar(s))
      return(NA_character_)
    if (grepl("^[A-Za-z]+[ _/-]", s)) s <- strsplit(s, "[ _/-]")[[1]][1]
    s <- sub("_[A-Z]{1,2}$", "", s)
    s <- tolower(gsub("[^A-Za-z]", "", s))
    if (nzchar(s)) s else NA_character_
  }, character(1), USE.NAMES = FALSE)
}
#' Sum rows that share the same (non-NA) key.
collapse_by_key <- function(m, keys) {
  ok <- !is.na(keys)
  rowsum(m[ok, , drop = FALSE], keys[ok])
}
