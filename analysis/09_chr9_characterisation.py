"""Characterisation of the chr9:120-150 Mb candidate region from the SNP data.

(1) internal linkage of the block (mean absolute pairwise correlation of SNP
frequencies among populations) compared with control windows and the genome;
(2) the leading principal component of the block and its correlation with
temperature; (3) population ordering along that axis (gradient versus a
two-haplotype split); (4) overlap of the block axis with the genome-wide
background structure axes.
"""
import numpy as np
from common import SNPMatrix, climate, env_vectors

snp = SNPMatrix("teosinte")
pops = snp.pops
bio1, _, _ = env_vectors(pops, climate())
Y, chrom, pos = snp.Y, snp.chrom, snp.pos
print("chr9 SNPs span %d-%d Mb (B73 RefGen_V4)" % (pos[chrom == "9"].min() / 1e6, pos[chrom == "9"].max() / 1e6))


def mean_abs_corr(M, cap=400, seed=0):
    rng = np.random.default_rng(seed)
    idx = rng.choice(M.shape[1], min(cap, M.shape[1]), replace=False)
    X = M[:, idx]
    Xs = X - X.mean(0)
    s = Xs.std(0)
    keep = s > 1e-9
    Xs = Xs[:, keep] / s[keep]
    C = np.corrcoef(Xs.T)
    iu = np.triu_indices(C.shape[0], 1)
    return np.nanmean(np.abs(C[iu]))


blk = snp.chr9_block
ctrl9 = (chrom == "9") & (pos >= 80e6) & (pos < 110e6)
ctrl1 = (chrom == "1") & (pos >= 120e6) & (pos < 150e6)
print("\n(1) internal linkage, mean |r| among SNP pairs (400 SNPs subsampled)")
print("    chr9:120-150 Mb (candidate): |r| = %.3f (n = %d)" % (mean_abs_corr(Y[:, blk]), blk.sum()))
print("    chr9:80-110 Mb (control):    |r| = %.3f (n = %d)" % (mean_abs_corr(Y[:, ctrl9]), ctrl9.sum()))
print("    chr1:120-150 Mb (control):   |r| = %.3f (n = %d)" % (mean_abs_corr(Y[:, ctrl1]), ctrl1.sum()))
print("    genome-wide random:          |r| = %.3f" % mean_abs_corr(Y))

Xb = Y[:, blk]
Xs = (Xb - Xb.mean(0)) / (Xb.std(0) + 1e-9)
U, S, Vt = np.linalg.svd(Xs, full_matrices=False)
ev = S ** 2 / np.sum(S ** 2)
b1 = U[:, 0]
print("\n(2) block PC1 explains %.1f%% of variance (PC2 %.1f%%); r(block PC1, bio1) = %+.3f"
      % (ev[0] * 100, ev[1] * 100, np.corrcoef(b1, bio1)[0, 1]))
print("\n(3) populations ordered along block PC1")
for i in np.argsort(b1):
    print("    %-18s PC1 = %+.3f  bio1 = %.1f" % (pops[i], b1[i], bio1[i]))
bg = snp.background_pcs(k=3)
print("\n(4) overlap of the block axis with background structure axes (chr9 block excluded)")
for k in range(3):
    print("    |r(block PC1, bgPC%d)| = %.3f" % (k + 1, abs(np.corrcoef(b1, bg[:, k])[0, 1])))
