"""SNP versus TE contrast under matched structure control (teosinte).

Applies the same conditioning schemes (ancestry; ancestry + background
structure PCs; ancestry + geography) to the genome-wide SNP matrix, the
chr9:120-150 Mb SNP block, the fully observed TE insertions and count-matched
random SNP subsets.  Background PCs are computed from SNPs outside the chr9
block.  This is the source of Table S2.
"""
import numpy as np
from common import SNPMatrix, TEMatrix, climate, teosinte_ancestry, env_vectors, rda1

snp = SNPMatrix("teosinte")
pops = snp.pops
te = TEMatrix("teosinte")
Yte = te.reorder(pops)[:, te.clean]
n_te_all = te.Y.shape[1]
clim = climate()
bio1, lat, lon = env_vectors(pops, clim)
anc = np.array([teosinte_ancestry()[p] for p in pops])
bg = snp.background_pcs(k=2)
blk = snp.chr9_block
print("SNP %dx%d | TE fully observed %dx%d" % (snp.Y.shape[0], snp.Y.shape[1], Yte.shape[0], Yte.shape[1]))

schemes = {"ancestry": anc.reshape(-1, 1),
           "ancestry + bgPC1 + bgPC2": np.column_stack([anc, bg]),
           "ancestry + lat + lon": np.column_stack([anc, lat, lon])}


def count_matched(Z, reps=5):
    rng = np.random.default_rng(3)
    ps = []
    for _ in range(reps):
        idx = rng.choice(snp.Y.shape[1], n_te_all, replace=False)
        ps.append(rda1(snp.Y[:, idx], bio1, Z, nperm=999, seed=rng.integers(1e6))[1])
    return np.mean(ps)


print("\n%-26s | %-16s | %-16s | %-16s | %-14s" % ("conditioning", "SNP genome-wide", "SNP chr9 block", "TE (%d)" % Yte.shape[1], "SNP count-matched"))
print("-" * 100)
for name, Z in schemes.items():
    Fa, pa = rda1(snp.Y, bio1, Z)
    Fc, pc = rda1(snp.Y[:, blk], bio1, Z)
    Ft, pt = rda1(Yte, bio1, Z)
    pcm = count_matched(Z)
    print("%-26s | F=%.2f p=%.4f | F=%.2f p=%.4f | F=%.2f p=%.4f | p=%.4f" % (name, Fa, pa, Fc, pc, Ft, pt, pcm))
