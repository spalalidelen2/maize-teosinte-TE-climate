"""Count-matched SNP control.

Repeats the teosinte SNP RDA (bio1 | ancestry) on five random subsets of 2,871
SNPs, matching the number of polymorphic TE insertions, to show that the
SNP-climate signal is not a feature-count effect.
"""
import numpy as np
from common import SNPMatrix, TEMatrix, climate, teosinte_ancestry, rda

snp = SNPMatrix("teosinte")
n_te = TEMatrix("teosinte").Y.shape[1]
pops = snp.pops
clim = climate()
bio1 = np.array([clim[p]["bio1"] for p in pops]).reshape(-1, 1)
anc = teosinte_ancestry()
A = np.array([anc[p] for p in pops]).reshape(-1, 1)

print("SNP RDA (bio1 | ancestry) on random subsets of %d SNPs, five replicates:" % n_te)
ps = []
for s in range(5):
    idx = np.random.default_rng(s).choice(snp.Y.shape[1], n_te, replace=False)
    _, F, p = rda(snp.Y[:, idx], bio1, A, nperm=2999, seed=100 + s)
    ps.append(p)
    print("  replicate %d: F = %.2f  p = %.4f" % (s + 1, F, p))
print("  mean p = %.4f" % np.mean(ps))
