"""Teosinte TE insertions versus climate: partial redundancy analysis.

Tests the 2,871 polymorphic TE insertions (undetected populations scored 0) and
the 929 insertions measured in all 17 populations against the first two climate
PCs and against mean annual temperature (bio1), conditioning on per-population
maize ancestry.  Also reports the per-locus partial correlation with
temperature (Benjamini-Hochberg FDR) and the unconditioned tests.
"""
import numpy as np
from common import TEMatrix, climate, climate_pcs, teosinte_ancestry, rda, perlocus_p, bh_fdr

te = TEMatrix("teosinte")
pops = te.pops
clim = climate()
pcs, ve = climate_pcs(pops, clim)
bio1 = np.array([clim[p]["bio1"] for p in pops])
anc = teosinte_ancestry()
A = np.array([anc[p] for p in pops]).reshape(-1, 1)
print("Climate PCA, variance explained by PC1-3: %.3f %.3f %.3f" % tuple(ve[:3]))

for label, mask in [("All polymorphic insertions (undetected = 0)", np.ones(te.Y.shape[1], bool)),
                    ("Fully observed insertions (nDetected = 17)", te.clean)]:
    Y = te.Y[:, mask]
    print("\n== %s: %d loci ==" % (label, mask.sum()))
    r2, F, p = rda(Y, pcs[:, :2], A)
    print("  RDA ~ climate PC1+PC2 | ancestry:  R2c = %.3f  F = %.2f  p = %.4f" % (r2, F, p))
    r2, F, p = rda(Y, bio1.reshape(-1, 1), A)
    print("  RDA ~ bio1 | ancestry:             R2c = %.3f  F = %.2f  p = %.4f" % (r2, F, p))

print("\n== Per-locus tests and unconditioned RDA, fully observed insertions ==")
Y = te.Y[:, te.clean]
p = perlocus_p(Y, bio1, A)
q = bh_fdr(p)
print("  per-locus partial correlation with bio1 | ancestry: min p = %.4g; FDR < 0.05: %d; FDR < 0.10: %d"
      % (p.min(), (q < 0.05).sum(), (q < 0.10).sum()))
r2, F, pv = rda(Y, bio1.reshape(-1, 1), None, seed=2)
print("  unconditioned RDA ~ bio1:     R2 = %.3f  F = %.2f  p = %.4f" % (r2, F, pv))
r2, F, pv = rda(Y, pcs[:, :2], None, seed=2)
print("  unconditioned RDA ~ PC1+PC2:  R2 = %.3f  F = %.2f  p = %.4f" % (r2, F, pv))

print("\n== Detection status of the top-ranked loci in the full set ==")
from collections import Counter
from common import read_tsv, DATA
pall = perlocus_p(te.Y, bio1, A)
top = np.argsort(pall)[:20]
print("  per-locus partial correlation (bio1 | ancestry), top 20: nDetected =",
      dict(sorted(Counter(te.nDetected[top]).items())))
gea = DATA / "te_teosinte" / "gea_results.tsv"
if gea.exists():
    h, rows = read_tsv(gea)
    ci = {c: i for i, c in enumerate(h)}
    rows.sort(key=lambda r: float(r[ci["lfmm_p"]]))
    print("  LFMM (scripts/teosinte/11_gea.R, full set), top 20 by p: nDetected =",
          dict(sorted(Counter(int(r[ci["nDetected"]]) for r in rows[:20]).items())))
    print("  LFMM candidates (cand_lfmm): nDetected =",
          dict(sorted(Counter(int(r[ci["nDetected"]]) for r in rows if r[ci["cand_lfmm"]] == "TRUE").items())))
