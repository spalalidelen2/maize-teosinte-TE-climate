#!/usr/bin/env Rscript
# Build the teosinte TE-insertion frequency matrix from the PoPoolationTE2
# insertion calls (hilo.teinsertions, one row per sample x insertion).
#   - sample index -> population from ppileup/pop_order.txt
#   - satellite and tandem repeat classes (knob, Cent, rDNA, subtelomere) removed
#   - calls of the same TE family within 100 bp clustered into one locus
#   - frequency per locus and population (mean over clustered calls);
#     populations without a call are assigned frequency 0
#   - polymorphic loci retained: variance > 0 and detected in at least 2 populations
#   Rscript scripts/teosinte/10_build_te_matrix.R [project_dir]
# Output: te_out/te_freq_matrix_all.tsv (all loci) and te_out/te_freq_matrix.tsv (polymorphic)
suppressMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
proj <- if (length(args) >= 1) args[1] else "."
ti <- fread(file.path(proj, "te_out/hilo.teinsertions"), header = FALSE)
setnames(ti, c("sample", "chr", "pos", "dot", "family", "order", "support", "strand", "freq"))
po <- fread(file.path(proj, "ppileup/pop_order.txt"), header = FALSE)     # V1 = index, V2 = BAM path
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
fwrite(mat,  file.path(proj, "te_out/te_freq_matrix_all.tsv"), sep = "\t")
fwrite(matf, file.path(proj, "te_out/te_freq_matrix.tsv"),     sep = "\t")
cat("populations:", length(popcols), "\n")
cat("clustered loci:", nrow(mat), "| polymorphic:", nrow(matf),
    "| fully observed:", sum(matf$nDetected == length(popcols)), "\n")
print(sort(table(matf$order), decreasing = TRUE))
