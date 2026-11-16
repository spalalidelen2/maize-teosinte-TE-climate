# Supplementary tables

Values are produced by the scripts named in each caption; the tables are provided as
`TableS1.tsv` to `TableS6.tsv` in this directory.

**Table S1** (`analysis/06_snp_noisematch.py`). Noise-matched SNP control, teosinte. SNP
frequencies were degraded to TE-like measurement (10-fold binomial resampling and a 0.20
detection floor) and the partial RDA (bio1, conditioned on ancestry) and the per-locus FDR test
re-run over 20 replicates. Columns: SNP data; RDA F; RDA p (median over replicates); replicates
with p < 0.05; per-locus FDR < 0.05 SNPs (median). The reference rows are the undegraded SNPs and
the fully observed TE insertions.

**Table S2** (`analysis/08_contrast_under_structure.py`, `analysis/07_structure_lambda.py`).
Partial RDA (pseudo-F, 1,999 permutations) of temperature (bio1) on each teosinte frequency
matrix under three conditioning schemes. Background structure axes (bgPC1, bgPC2) are the leading
principal components of the SNP matrix with the chr9:120 to 150 Mb region excluded. The SNP
count-matched column is the mean p of five random 2,871-SNP subsets. lambda_GC is the genomic
inflation factor of the per-locus SNP test (permuted-phenotype calibration 0.92).

**Table S3** (`analysis/12_lfmm2_teosinte.R`). LFMM2 (LEA) of the teosinte SNP and TE frequency
matrices on temperature, K = 1 to 5 latent factors, genomic-control calibration; SNPs filtered to
the 78,258 loci with among-population variance. GIF, genomic inflation factor.

**Table S4** (`analysis/14_maize_noisematch.py`). Noise-matched SNP control, maize. SNP
frequencies were degraded to maize-TE measurement quality (8-fold binomial resampling and a 0.20
detection floor) and the unconditioned RDA on temperature re-run over 20 replicates.

**Table S5** (`analysis/16_lfmm2_maize.R`). LFMM2 (LEA) of the maize SNP and TE frequency
matrices on temperature, K = 1 to 5 latent factors, genomic-control calibration.

**Table S6** (`analysis/13_maize_gea.py`). Partial RDA (pseudo-F, 1,999 permutations) of
temperature on the maize TE and SNP frequency matrices under three conditioning schemes;
covariate is mexicana introgression. The SNP count-matched column is the mean p of five random
subsets matched to the maize TE locus count.
