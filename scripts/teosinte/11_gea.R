#!/usr/bin/env Rscript
# Initial genotype-environment association of the teosinte TE-insertion
# frequency matrix on the cluster (R packages vegan and LEA); produces
# te_out/gea_results.tsv (data/te_teosinte/gea_results.tsv).
#   - redundancy analysis of the 2,871 polymorphic insertions on the first two
#     climate PCs, conditioned on maize ancestry (vegan::rda, anova.cca with
#     999 permutations); candidate loci = |RDA1 loading| > 3 SD
#   - LFMM2 (LEA) of the same matrix on climate PC1 with K = 2 latent factors;
#     candidate loci = Benjamini-Hochberg adjusted p < 0.1
#   Rscript scripts/teosinte/11_gea.R [project_dir]
# Inputs: te_out/te_freq_matrix.tsv, ppileup/pop_order.txt,
#         metadata/climate_table_hilo.tsv, metadata/maize_ancestry_perpop.tsv
suppressMessages({library(data.table); library(vegan); library(LEA)})
args <- commandArgs(trailingOnly = TRUE)
proj <- if (length(args) >= 1) args[1] else "."
mat <- fread(file.path(proj, "te_out/te_freq_matrix.tsv"))
po  <- fread(file.path(proj, "ppileup/pop_order.txt"), header = FALSE)
po[, pop := sub("\\.bam$", "", basename(V2))]
popcols <- po$pop
clim <- fread(file.path(proj, "metadata/climate_table_hilo.tsv"))
anc  <- fread(file.path(proj, "metadata/maize_ancestry_perpop.tsv"))
clim[, population := gsub(" ", "_", population)]
anc[,  population := gsub(" ", "_", population)]
clim <- clim[match(popcols, population)]; anc <- anc[match(popcols, population)]
stopifnot(all(clim$population == popcols), all(anc$population == popcols))

Y <- t(as.matrix(mat[, ..popcols])); rownames(Y) <- popcols
cat("Y (populations x loci):", paste(dim(Y), collapse = " x "), "\n")

bio <- scale(as.matrix(clim[, paste0("bio", 1:19), with = FALSE]))
pca <- prcomp(bio); climPC <- pca$x[, 1:2]
cat("climate PCA, variance explained (PC1-4):", paste(round(summary(pca)$importance[2, 1:4], 3), collapse = " "), "\n")
ancv <- anc$mean_maize_ancestry

df <- data.frame(PC1 = climPC[, 1], PC2 = climPC[, 2], anc = ancv)
rda1 <- rda(Y ~ PC1 + PC2 + Condition(anc), data = df)
print(rda1)
set.seed(1); print(anova.cca(rda1, permutations = 999))
load1 <- scores(rda1, choices = 1, display = "species", scaling = 3)[, 1]
cand_rda <- abs(load1 - mean(load1)) > 3 * sd(load1)
cat("RDA candidate loci (|RDA1 loading| > 3 SD):", sum(cand_rda), "\n")

K <- 2
mod <- lfmm2(input = Y, env = as.matrix(climPC[, 1]), K = K)
tt  <- lfmm2.test(mod, input = Y, env = as.matrix(climPC[, 1]), full = TRUE)
p   <- as.numeric(tt$pvalues); padj <- p.adjust(p, "BH")
cand_lfmm <- padj < 0.1
cat("LFMM2 candidate loci (BH < 0.1, climate PC1, K = 2):", sum(cand_lfmm, na.rm = TRUE), "\n")
cand_both <- cand_rda & cand_lfmm
cat("candidates in both methods:", sum(cand_both, na.rm = TRUE), "\n")

res <- mat[, .(chr, pos, family, order, nDetected)]
res[, rda_load := load1][, lfmm_p := p][, lfmm_padj := padj]
res[, cand_rda := cand_rda][, cand_lfmm := cand_lfmm][, cand_both := cand_both]
fwrite(res, file.path(proj, "te_out/gea_results.tsv"), sep = "\t")
cat("written te_out/gea_results.tsv\n")
