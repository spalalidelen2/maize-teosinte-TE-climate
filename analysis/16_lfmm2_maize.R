#!/usr/bin/env Rscript
# Latent-factor mixed models (LFMM2, R package LEA) of maize SNP and TE
# frequencies on mean annual temperature, K = 1..5, genomic-control
# calibration (Table S5).  Same procedure as analysis/12_lfmm2_teosinte.R.
#
# Run from the repository root:  Rscript analysis/16_lfmm2_maize.R
# Inputs: data/snp_matrix_maize.tsv, data/te_maize/te_freq_matrix.tsv,
#         data/metadata/mexicana_ancestry_maize_perpop.tsv (locality, popN, bio1)
suppressMessages(library(LEA))
set.seed(1)
args <- commandArgs(trailingOnly = FALSE)
root <- normalizePath(file.path(dirname(sub("--file=", "", args[grep("--file=", args)])), ".."))
rd <- function(f) read.delim(file.path(root, f), check.names = FALSE, stringsAsFactors = FALSE)

cov <- rd("data/metadata/mexicana_ancestry_maize_perpop.tsv")
bio1v <- setNames(as.numeric(cov$bio1), cov$population)
popN2loc <- setNames(cov$population, as.character(cov$popN))

## TE matrix, fully observed loci (nDetected == 14)
te <- rd("data/te_maize/te_freq_matrix.tsv")
metacols <- c("locus","chr","family","order","pos","nsig","nDetected","variance")
tepops <- setdiff(colnames(te), metacols)
te <- te[te$nDetected == length(tepops), ]
pops <- tepops
Yte <- t(as.matrix(te[, pops]))
X <- matrix(bio1v[pops], ncol = 1)

## SNP matrix (columns: chr, pos, pop360..pop374)
snp <- rd("data/snp_matrix_maize.tsv")
snp_loc <- unname(popN2loc[sub("pop", "", colnames(snp)[3:ncol(snp)])])
chrom <- as.character(snp[[1]]); pos <- as.numeric(snp[[2]])
Ysnp <- t(as.matrix(snp[, 3:ncol(snp)])); Ysnp[Ysnp < 0] <- 0; Ysnp[Ysnp > 1] <- 1
rownames(Ysnp) <- snp_loc; Ysnp <- Ysnp[pops, ]
keepS <- apply(Ysnp, 2, sd) > 1e-6; Ysnp <- Ysnp[, keepS]; chrom <- chrom[keepS]; pos <- pos[keepS]
keepT <- apply(Yte, 2, sd) > 1e-6;  Yte <- Yte[, keepT]
chr9 <- (chrom == "9" & pos >= 120e6 & pos < 150e6)
cat(sprintf("maize: SNP %d loci, TE %d loci, %d pops (aligned = %s)\n\n",
            ncol(Ysnp), ncol(Yte), nrow(Ysnp), all(rownames(Ysnp) == pops)))

run <- function(Y, label, snp = FALSE) {
  cat(sprintf("=== %s ===\n", label))
  for (K in 1:5) {
    m <- lfmm2(input = Y, env = X, K = K)
    tt <- lfmm2.test(object = m, input = Y, env = X, linear = TRUE, genomic.control = TRUE)
    p <- as.numeric(tt$pvalues); z <- as.numeric(tt$zscores); g <- as.numeric(tt$gif)
    if (is.na(g)) g <- median(z^2, na.rm = TRUE) / qchisq(0.5, 1)
    q <- p.adjust(p, "BH"); ns <- sum(q < 0.05, na.rm = TRUE)
    if (snp) {
      cat(sprintf("  K=%d: GIF=%.2f FDR<0.05=%d (chr9=%d) min_p=%.1e\n", K, g, ns, sum(q < 0.05 & chr9, na.rm = TRUE), min(p, na.rm = TRUE)))
    } else {
      cat(sprintf("  K=%d: GIF=%.2f FDR<0.05=%d min_p=%.1e\n", K, g, ns, min(p, na.rm = TRUE)))
    }
  }
  cat("\n")
}
run(Ysnp, "maize SNP ~ temperature (LFMM2)", snp = TRUE)
run(Yte,  "maize TE (fully observed) ~ temperature (LFMM2)")
