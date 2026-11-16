"""Robustness of the teosinte SNP-climate association.

(A) all SNPs outside the chr9:120-150 Mb candidate region, (B) the candidate
region alone, and (C) leave-one-population-out re-tests of the RDA
(bio1 | ancestry) for the genome-wide set and for the chr9 region.
"""
import numpy as np
from common import SNPMatrix, climate, teosinte_ancestry, rda

snp = SNPMatrix("teosinte")
pops = snp.pops
clim = climate()
bio1 = np.array([clim[p]["bio1"] for p in pops])
anc = teosinte_ancestry()
a = np.array([anc[p] for p in pops])
blk = snp.chr9_block


def test(Y, idx=None, seed=1):
    if idx is not None:
        Y, x, z = Y[idx], bio1[idx], a[idx]
    else:
        x, z = bio1, a
    _, F, p = rda(Y, x.reshape(-1, 1), z.reshape(-1, 1), nperm=2999, seed=seed)
    return F, p


print("(A) all SNPs excluding chr9:120-150 Mb (n = %d): F = %.2f  p = %.4f" % ((~blk).sum(), *test(snp.Y[:, ~blk])))
print("(B) chr9:120-150 Mb only (n = %d):               F = %.2f  p = %.4f" % (blk.sum(), *test(snp.Y[:, blk])))
F, p = test(snp.Y)
print("(C) leave-one-population-out, all SNPs; all %d populations: F = %.2f  p = %.4f" % (len(pops), F, p))
res = []
for i, pop in enumerate(pops):
    idx = np.array([j for j in range(len(pops)) if j != i])
    res.append((pop,) + test(snp.Y, idx, seed=10 + i))
res.sort(key=lambda z: -z[2])
print("    five populations whose removal weakens the signal most:")
for pop, F, p in res[:5]:
    print("      without %-18s F = %.2f  p = %.4f" % (pop, F, p))
print("    chr9 region, populations whose removal gives p > 0.01:")
weak = False
for i, pop in enumerate(pops):
    idx = np.array([j for j in range(len(pops)) if j != i])
    F, p = test(snp.Y[:, blk], idx, seed=50 + i)
    if p > 0.01:
        weak = True
        print("      without %-18s F = %.2f  p = %.4f" % (pop, F, p))
if not weak:
    print("      none (p <= 0.01 in every leave-one-out subset)")
