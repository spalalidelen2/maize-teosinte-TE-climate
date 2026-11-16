#!/usr/bin/env Rscript
# Build the maize TE-insertion frequency matrix (14 populations) from the
# PoPoolationTE2 calls, with the same procedure as the teosinte matrix
# (scripts/teosinte/10_build_te_matrix.R): satellite classes removed, calls of
# the same family within 100 bp clustered into one locus, undetected
# populations scored 0, polymorphic loci = variance > 0 and detected in at
# least 2 populations.
#   Rscript scripts/maize/06_build_te_matrix.R [project_dir]
# Output: maize/te_out/te_freq_matrix_all.tsv and maize/te_out/te_freq_matrix.tsv
suppressMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
proj <- if (length(args) >= 1) args[1] else "."
ti <- fread(file.path(proj, "maize/te_out/maize.teinsertions"), header = FALSE)
setnames(ti, c("sample", "chr", "pos", "dot", "family", "order", "support", "strand", "freq"))
po <- fread(file.path(proj, "maize/ppileup/pop_order.txt"), header = FALSE)
po[, pop := sub("\\.bam$", "", basename(V2))]
popmap <- setNames(po$pop, as.character(po$V1))
ti[, pop := popmap[as.character(sample)]]
ti <- ti[!order %in% c("knob", "Cent", "rDNA", "subtelomere")]
setorder(ti, chr, family, pos)
ti[, gap := pos - shift(pos)]
ti[, newcl := is.na(gap) | chr != shift(chr) | family != shift(family) | gap > 100]
ti[, locus := cumsum(newcl)]
loci <- ti[, .(chr = chr[1], family = family[1], order = order[1], pos = min(pos), nsig = .N), by = locus]
lp <- ti[, .(freq = mean(freq)), by = .(locus, pop)]
popcols <- po$pop
mat <- dcast(lp, locus ~ pop, value.var = "freq", fill = 0)
mat <- merge(loci, mat, by = "locus")
m <- as.matrix(mat[, ..popcols])
mat[, nDetected := rowSums(m > 0)]
mat[, variance := apply(m, 1, var)]
matf <- mat[variance > 0 & nDetected >= 2]
fwrite(mat,  file.path(proj, "maize/te_out/te_freq_matrix_all.tsv"), sep = "\t")
fwrite(matf, file.path(proj, "maize/te_out/te_freq_matrix.tsv"),     sep = "\t")
cat("populations:", length(popcols), "\n")
cat("clustered loci:", nrow(mat), "| polymorphic:", nrow(matf),
    "| fully observed:", sum(matf$nDetected == length(popcols)), "\n")
print(sort(table(matf$order), decreasing = TRUE))
