"""Teosinte SNP benchmark: the same RDA framework applied to genome-wide SNPs.

Positive control (SNP frequencies structured by maize ancestry), SNP ~ climate
PC1+PC2 and SNP ~ bio1 conditioned on ancestry, the unconditioned test, and the
per-SNP partial correlation with temperature (FDR).
"""
import numpy as np
from common import SNPMatrix, climate, climate_pcs, teosinte_ancestry, rda, perlocus_p, bh_fdr

snp = SNPMatrix("teosinte")
pops = snp.pops
clim = climate()
pcs, _ = climate_pcs(pops, clim)
bio1 = np.array([clim[p]["bio1"] for p in pops])
anc = teosinte_ancestry()
A = np.array([anc[p] for p in pops]).reshape(-1, 1)
print("SNP matrix: %d SNPs x %d populations\n" % (snp.Y.shape[1], snp.Y.shape[0]))

r2, F, p = rda(snp.Y, A, None)
print("Positive control, SNP ~ maize ancestry (unconditioned): R2 = %.3f  F = %.2f  p = %.4f" % (r2, F, p))
r2, F, p = rda(snp.Y, pcs[:, :2], A)
print("SNP ~ climate PC1+PC2 | ancestry: R2c = %.3f  F = %.2f  p = %.4f" % (r2, F, p))
r2, F, p = rda(snp.Y, bio1.reshape(-1, 1), A)
print("SNP ~ bio1 | ancestry:            R2c = %.3f  F = %.2f  p = %.4f" % (r2, F, p))
r2, F, p = rda(snp.Y, bio1.reshape(-1, 1), None)
print("SNP ~ bio1 (unconditioned):       R2 = %.3f  F = %.2f  p = %.4f" % (r2, F, p))

pv = perlocus_p(snp.Y, bio1, A)
q = bh_fdr(pv)
print("\nPer-SNP partial correlation with bio1 | ancestry: min p = %.2g; FDR < 0.10: %d; FDR < 0.05: %d"
      % (pv.min(), (q < 0.10).sum(), (q < 0.05).sum()))
