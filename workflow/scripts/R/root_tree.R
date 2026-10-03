#!/usr/bin/env Rscript
# Midpoint-root an unrooted FastTree tree (UniFrac / Faith's PD need a rooted tree).
suppressPackageStartupMessages({library(optparse); library(ape); library(phangorn)})
opt <- parse_args(OptionParser(option_list = list(make_option("--tree"), make_option("--out"))))
tr <- read.tree(opt$tree)
tr <- midpoint(multi2di(tr, random = FALSE))
tr$edge.length[tr$edge.length < 0] <- 0
write.tree(tr, opt$out)
message("rooted tree with ", length(tr$tip.label), " tips")
