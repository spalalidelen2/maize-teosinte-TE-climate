"""Teosinte TE-insertion frequencies versus the maize-ancestry gradient.

Tests whether TE frequencies are structured by maize introgression (RDA on
ancestry, all and fully observed loci), whether per-population TE load tracks
ancestry, and the per-locus correlation with ancestry (FDR).
"""
import numpy as np
from scipy import stats
from common import TEMatrix, climate, teosinte_ancestry, rda, bh_fdr

te = TEMatrix("teosinte")
pops = te.pops
anc = np.array([teosinte_ancestry()[p] for p in pops])
bio1 = np.array([climate()[p]["bio1"] for p in pops])
print("Maize ancestry range %.3f-%.3f; r(ancestry, bio1) = %.3f" % (anc.min(), anc.max(), np.corrcoef(anc, bio1)[0, 1]))

for label, m in [("all polymorphic (%d)" % te.Y.shape[1], np.ones(te.Y.shape[1], bool)),
                 ("fully observed (%d)" % te.clean.sum(), te.clean)]:
    r2, F, p = rda(te.Y[:, m], anc.reshape(-1, 1), None, seed=3)
    print("  RDA TE ~ ancestry, %s: R2 = %.3f  F = %.2f  p = %.4f" % (label, r2, F, p))

Yc = te.Y[:, te.clean]
print("  r(population mean TE frequency, ancestry) = %.2f (p = %.3f)" % stats.pearsonr(Yc.mean(1), anc))
print("  r(number of detected insertions, ancestry) = %.2f (p = %.3f)" % stats.pearsonr((te.Y > 0).sum(1), anc))
ps = []
for j in range(Yc.shape[1]):
    y = Yc[:, j]
    if y.std() < 1e-9:
        ps.append(1.0)
        continue
    r = np.corrcoef(y, anc)[0, 1]
    t = r * np.sqrt((len(y) - 2) / max(1e-9, 1 - r * r))
    ps.append(2 * stats.t.sf(abs(t), len(y) - 2))
q = bh_fdr(np.array(ps))
print("  per-locus correlation with ancestry (fully observed): min p = %.2g, FDR < 0.10: %d, FDR < 0.05: %d"
      % (min(ps), (q < 0.10).sum(), (q < 0.05).sum()))
