#!/usr/bin/env Rscript
# Latent-factor mixed models (LFMM2, R package LEA) of teosinte SNP and TE
# frequencies on mean annual temperature, K = 1..5 latent factors, with
# genomic-control calibration.  Reports the genomic inflation factor (GIF),
# the number of loci at FDR < 0.05 and, for SNPs, how many of those fall in
# the chr9:120-150 Mb candidate region (Table S3).
#
# Run from the repository root:  Rscript analysis/12_lfmm2_teosinte.R
# Inputs: data/snp_matrix_teosinte.tsv, data/te_teosinte/te_freq_matrix.tsv,
#         data/metadata/climate_table_hilo.tsv
suppressMessages(library(LEA))
set.seed(1)
args <- commandArgs(trailingOnly = FALSE)
root <- normalizePath(file.path(dirname(sub("--file=", "", args[grep("--file=", args)])), ".."))
rd <- function(f) read.delim(file.path(root, f), check.names = FALSE, stringsAsFactors = FALSE)
nrm <- function(x) gsub(" ", "_", x)

## climate (the K latent factors absorb ancestry and structure, so no explicit covariate)
cl <- rd("data/metadata/climate_table_hilo.tsv")
bio1v <- setNames(as.numeric(cl$bio1), nrm(cl$population))

## SNP matrix (columns: chr, pos, pop18..pop35; Figshare population numbers)
pn2loc <- c('18'='Nabogame','19'='Penjamillo','20'='Malinalco','21'='El_Porvenir',
            '22'='Amecameca','23'='Tlapala','24'='Jicaltepec','25'='Cocotilan',
            '26'='Opopeo','27'='Amatlan','28'='Puruandiro','29'='Ixtlan','30'='San_Pedro',
            '31'='Tenango_del_Aire','33'='Puerta_Encantada','34'='Santa_Clara','35'='Xochimilco')
snp <- rd("data/snp_matrix_teosinte.tsv")
pops <- unname(pn2loc[sub("pop", "", colnames(snp)[3:ncol(snp)])])
chrom <- as.character(snp[[1]]); pos <- as.numeric(snp[[2]])
Ysnp <- t(as.matrix(snp[, 3:ncol(snp)])); Ysnp[Ysnp < 0] <- 0; Ysnp[Ysnp > 1] <- 1
rownames(Ysnp) <- pops
X <- matrix(bio1v[pops], ncol = 1)
keepS <- apply(Ysnp, 2, sd) > 1e-6          # drop loci without among-population variance
Ysnp <- Ysnp[, keepS]; chrom <- chrom[keepS]; pos <- pos[keepS]
cat(sprintf("SNP: %d of %d loci retained after the variance filter\n", sum(keepS), length(keepS)))

## TE matrix, fully observed loci (nDetected == 17)
te <- rd("data/te_teosinte/te_freq_matrix.tsv")
metacols <- c("locus","chr","family","order","pos","nsig","nDetected","variance")
tepops <- setdiff(colnames(te), metacols)
te <- te[te$nDetected == length(tepops), ]
Yte <- t(as.matrix(te[, tepops])); rownames(Yte) <- tepops
Yte <- Yte[pops, ]
keepT <- apply(Yte, 2, sd) > 1e-6; Yte <- Yte[, keepT]
cat(sprintf("SNP %d loci x %d pops | TE %d loci | pops aligned = %s\n\n",
            ncol(Ysnp), nrow(Ysnp), ncol(Yte), all(rownames(Ysnp) == rownames(Yte))))
chr9 <- (chrom == "9" & pos >= 120e6 & pos < 150e6)

run_lfmm <- function(Y, label, is_snp = FALSE) {
  cat(sprintf("=== %s ===\n", label))
  for (K in 1:5) {
    mod <- lfmm2(input = Y, env = X, K = K)
    tst <- lfmm2.test(object = mod, input = Y, env = X, linear = TRUE, genomic.control = TRUE)
    p <- as.numeric(tst$pvalues); z <- as.numeric(tst$zscores)
    gif <- as.numeric(tst$gif)
    if (is.na(gif)) gif <- median(z^2, na.rm = TRUE) / qchisq(0.5, 1)
    q <- p.adjust(p, method = "BH")
    nsig <- sum(q < 0.05, na.rm = TRUE)
    if (is_snp) {
      cat(sprintf("  K=%d: GIF=%.2f  FDR<0.05=%d (chr9:120-150Mb=%d)  min p=%.1e\n",
                  K, gif, nsig, sum(q < 0.05 & chr9, na.rm = TRUE), min(p, na.rm = TRUE)))
    } else {
      cat(sprintf("  K=%d: GIF=%.2f  FDR<0.05=%d  min p=%.1e\n", K, gif, nsig, min(p, na.rm = TRUE)))
    }
  }
  cat("\n")
}
run_lfmm(Ysnp, "teosinte SNP ~ temperature (LFMM2)", is_snp = TRUE)
run_lfmm(Yte,  "teosinte TE (fully observed) ~ temperature (LFMM2)")
