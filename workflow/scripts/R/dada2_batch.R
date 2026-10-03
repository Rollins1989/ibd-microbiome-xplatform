#!/usr/bin/env Rscript
# DADA2 denoising of ONE sequencing batch/run (error models are run-specific).
# Output: sequence table (RDS) + per-sample read tracking + error-model plots.
suppressPackageStartupMessages({library(optparse); library(dada2)})
opt <- parse_args(OptionParser(option_list = list(
  make_option("--samples"), make_option("--r1"), make_option("--r2"),        # comma-separated, same order
  make_option("--filtdir"), make_option("--trunc-len", default = "240,200"), make_option("--max-ee", default = "2,2"),
  make_option("--trunc-q", type = "integer", default = 2L), make_option("--min-overlap", type = "integer", default = 12L),
  make_option("--max-mismatch", type = "integer", default = 0L), make_option("--pool", default = "false"),
  make_option("--threads", type = "integer", default = 1L), make_option("--seed", type = "integer", default = 1L),
  make_option("--out-seqtab"), make_option("--out-track"), make_option("--out-errplot"))),
  convert_hyphens_to_underscores = FALSE)
set.seed(opt$seed)
sid <- strsplit(opt$samples, ",")[[1]]; r1 <- strsplit(opt$r1, ",")[[1]]; r2 <- strsplit(opt$r2, ",")[[1]]
stopifnot(length(sid) == length(r1), length(r1) == length(r2))
num2 <- function(x) as.numeric(strsplit(x, ",")[[1]])
pool <- switch(opt$pool, "false" = FALSE, "true" = TRUE, "pseudo" = "pseudo", stop("--pool must be false|true|pseudo"))
mt <- if (opt$threads > 1) opt$threads else FALSE
dir.create(opt$filtdir, recursive = TRUE, showWarnings = FALSE)
fF <- file.path(opt$filtdir, paste0(sid, "_F.fastq.gz")); fR <- file.path(opt$filtdir, paste0(sid, "_R.fastq.gz"))
names(fF) <- names(fR) <- sid
flt <- filterAndTrim(r1, fF, r2, fR, truncLen = num2(opt$`trunc-len`), maxN = 0, maxEE = num2(opt$`max-ee`),
                     truncQ = opt$`trunc-q`, rm.phix = TRUE, compress = TRUE, multithread = mt)
rownames(flt) <- sid
ok <- file.exists(fF) & file.exists(fR) & flt[, "reads.out"] > 0
if (sum(ok) < 1) stop("no reads survived filtering; relax --trunc-len/--max-ee")
if (any(!ok)) message("samples with no reads after filtering (dropped): ", paste(sid[!ok], collapse = ", "))
fF <- fF[ok]; fR <- fR[ok]
errF <- learnErrors(fF, multithread = mt, verbose = FALSE)
errR <- learnErrors(fR, multithread = mt, verbose = FALSE)
pdf(opt$`out-errplot`, width = 8, height = 8); print(plotErrors(errF, nominalQ = TRUE)); print(plotErrors(errR, nominalQ = TRUE)); dev.off()
dF <- dada(fF, err = errF, pool = pool, multithread = mt, verbose = FALSE)
dR <- dada(fR, err = errR, pool = pool, multithread = mt, verbose = FALSE)
if (inherits(dF, "dada")) { dF <- list(dF); dR <- list(dR); names(dF) <- names(dR) <- names(fF) }
mg <- mergePairs(dF, fF, dR, fR, minOverlap = opt$`min-overlap`, maxMismatch = opt$`max-mismatch`, verbose = FALSE)
if (is.data.frame(mg)) { mg <- list(mg); names(mg) <- names(fF) }
st <- makeSequenceTable(mg)
getN <- function(x) sum(getUniques(x))
track <- data.frame(sample_id = sid, input = flt[, "reads.in"], filtered = flt[, "reads.out"],
                    denoised_fwd = NA_real_, denoised_rev = NA_real_, merged = NA_real_, batch_seqtab_reads = NA_real_)
m <- match(names(fF), sid)
track$denoised_fwd[m] <- vapply(dF, getN, numeric(1)); track$denoised_rev[m] <- vapply(dR, getN, numeric(1))
track$merged[m] <- vapply(mg, function(x) sum(x$abundance), numeric(1)); track$batch_seqtab_reads[m] <- rowSums(st)[names(fF)]
saveRDS(st, opt$`out-seqtab`)
write.table(track, opt$`out-track`, sep = "\t", quote = FALSE, row.names = FALSE)
message("batch done: ", nrow(st), " samples, ", ncol(st), " sequence variants")
