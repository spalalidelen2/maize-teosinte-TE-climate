"""Population structure and genomic inflation in the teosinte SNP-climate test.

(1) correlation of temperature with ancestry, background structure axes and
geography; (2) genomic inflation factor lambda_GC of the per-SNP test under
different conditioning schemes, with a permuted-phenotype calibration;
(3) structure-conditioned RDA of all SNPs and (4) of the chr9:120-150 Mb
candidate region.  Background structure axes are the leading population PCs
of the SNP matrix with the chr9 candidate region excluded.
"""
import numpy as np
from common import SNPMatrix, climate, teosinte_ancestry, env_vectors, rda1, perlocus_p, lambda_gc

snp = SNPMatrix("teosinte")
pops = snp.pops
n = len(pops)
clim = climate()
bio1, lat, lon = env_vectors(pops, clim)
anc = np.array([teosinte_ancestry()[p] for p in pops])
A = anc.reshape(-1, 1)
bg = snp.background_pcs(k=4)
Y = snp.Y
blk = snp.chr9_block

print("(1) correlation of temperature (bio1) with structure and geography")
for name, v in [("ancestry", anc), ("bgPC1", bg[:, 0]), ("bgPC2", bg[:, 1]), ("bgPC3", bg[:, 2]),
                ("latitude", lat), ("longitude", lon)]:
    print("    r(bio1, %-9s) = %+.3f" % (name, np.corrcoef(bio1, v)[0, 1]))

print("\n(2) genomic inflation lambda_GC of the per-SNP temperature test")
print("    bio1 | ancestry:                 %.3f" % lambda_gc(perlocus_p(Y, bio1, A)))
print("    bio1 unconditioned:              %.3f" % lambda_gc(perlocus_p(Y, bio1, None)))
print("    bio1 | ancestry + bgPC1 + bgPC2: %.3f" % lambda_gc(perlocus_p(Y, bio1, np.column_stack([anc, bg[:, :2]]))))
rng = np.random.default_rng(7)
lp = [lambda_gc(perlocus_p(Y, bio1[rng.permutation(n)], A)) for _ in range(10)]
print("    permuted bio1 (calibration):     %.3f (mean of 10 permutations)" % np.mean(lp))

schemes = [("ancestry", A),
           ("ancestry + bgPC1", np.column_stack([anc, bg[:, 0]])),
           ("ancestry + bgPC1 + bgPC2", np.column_stack([anc, bg[:, :2]])),
           ("ancestry + latitude + longitude", np.column_stack([anc, lat, lon]))]
print("\n(3) structure-conditioned RDA, all SNPs ~ bio1")
for name, Z in schemes:
    F, p = rda1(Y, bio1, Z)
    print("    bio1 | %-32s F = %.2f  p = %.4f" % (name + ":", F, p))
print("\n(4) structure-conditioned RDA, chr9:120-150 Mb ~ bio1 (%d SNPs)" % blk.sum())
for name, Z in [schemes[0], schemes[2], schemes[3]]:
    F, p = rda1(Y[:, blk], bio1, Z)
    print("    bio1 | %-32s F = %.2f  p = %.4f" % (name + ":", F, p))
