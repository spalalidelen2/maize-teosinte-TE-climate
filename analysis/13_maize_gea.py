"""Replication in domesticated maize: the same GEA framework on 14 sympatric
maize populations, with mexicana introgression as the covariate.

(1) positive control and the SNP versus TE contrast under three conditioning
schemes (unconditioned; introgression; introgression + background structure),
including count-matched SNP subsets (Table S6); (2) per-locus TE tests (FDR);
(3) the full versus fully observed TE sets; (4) leave-one-population-out
robustness of the one nominally significant TE test; (5) genomic inflation of
the per-SNP temperature test.
"""
import numpy as np
from common import SNPMatrix, TEMatrix, climate, maize_introgression, rda, rda1, perlocus_p, bh_fdr, lambda_gc

te = TEMatrix("maize")
pops = te.pops                      # population order of the TE matrix (alphabetical)
n = len(pops)
Yall = te.Y
Yte = Yall[:, te.clean]
snp = SNPMatrix("maize")
Ysnp = snp.reorder(pops)
clim = climate()
bio1 = np.array([clim[p]["bio1"] for p in pops])
intro = maize_introgression()
anc = np.array([intro[p][1] for p in pops])
A = anc.reshape(-1, 1)
none = np.zeros((n, 0))
blk = snp.chr9_block
bg = snp.background_pcs(k=2, pops=pops)
print("maize: %d populations | TE fully observed %d of %d | SNP %d (chr9 block %d)"
      % (n, Yte.shape[1], Yall.shape[1], Ysnp.shape[1], blk.sum()))
print("r(bio1, mexicana introgression) = %.3f" % np.corrcoef(bio1, anc)[0, 1])

_, F, p = rda(Ysnp, A, None, nperm=1999)
print("\npositive control, SNP ~ mexicana introgression (unconditioned): F = %.2f  p = %.4f" % (F, p))


def count_matched(Z, reps=5):
    rng = np.random.default_rng(3)
    ps = []
    for _ in range(reps):
        idx = rng.choice(Ysnp.shape[1], min(Yte.shape[1], Ysnp.shape[1]), replace=False)
        ps.append(rda1(Ysnp[:, idx], bio1, Z, nperm=999, seed=rng.integers(1e6))[1])
    return np.mean(ps)


schemes = [("unconditioned", none), ("introgression", A), ("introgression + bgPC1 + bgPC2", np.column_stack([anc, bg]))]
print("\n(1) SNP versus TE contrast, RDA on bio1")
print("%-30s | %-16s | %-16s | %-16s | %-14s" % ("conditioning", "SNP genome-wide", "SNP chr9 block", "TE (%d)" % Yte.shape[1], "SNP count-matched"))
print("-" * 104)
for name, Z in schemes:
    Fa, pa = rda1(Ysnp, bio1, Z)
    Fc, pc = rda1(Ysnp[:, blk], bio1, Z)
    Ft, pt = rda1(Yte, bio1, Z)
    print("%-30s | F=%.2f p=%.4f | F=%.2f p=%.4f | F=%.2f p=%.4f | p=%.4f" % (name, Fa, pa, Fc, pc, Ft, pt, count_matched(Z)))

print("\n(2) per-locus TE tests, fully observed insertions")
for name, Z in [("bio1 unconditioned", none), ("bio1 | introgression", A)]:
    p = perlocus_p(Yte, bio1, Z)
    q = bh_fdr(p)
    print("    %-22s FDR < 0.05: %d  FDR < 0.10: %d  min p = %.2g" % (name, (q < 0.05).sum(), (q < 0.10).sum(), p.min()))
p = perlocus_p(Ysnp, bio1, none)
print("    SNPs, bio1 unconditioned: FDR < 0.05: %d" % (bh_fdr(p) < 0.05).sum())

print("\n(3) all polymorphic versus fully observed TE sets, bio1 | introgression")
for name, Y in [("all %d" % Yall.shape[1], Yall), ("fully observed %d" % Yte.shape[1], Yte)]:
    F, p = rda1(Y, bio1, A, nperm=2999)
    print("    %-22s F = %.2f  p = %.4f" % (name, F, p))

print("\n(4) leave-one-population-out, TE ~ bio1 | introgression + bgPC1 + bgPC2")
Z = np.column_stack([anc, bg])
print("    all %d populations: p = %.4f" % (n, rda1(Yte, bio1, Z, nperm=2999)[1]))
res = []
for i in range(n):
    keep = np.array([j for j in range(n) if j != i])
    res.append((pops[i], rda1(Yte[keep], bio1[keep], Z[keep], nperm=1999, seed=5 + i)[1]))
res.sort(key=lambda z: -z[1])
print("    populations whose removal weakens the test most:")
for pop, pv in res[:4]:
    print("      without %-18s p = %.4f" % (pop, pv))

print("\n(5) genomic inflation lambda_GC of the per-SNP temperature test")
print("    bio1 unconditioned:                 %.3f" % lambda_gc(perlocus_p(Ysnp, bio1, none)))
print("    bio1 | introgression:               %.3f" % lambda_gc(perlocus_p(Ysnp, bio1, A)))
print("    bio1 | introgression + bgPC1 + 2:   %.3f" % lambda_gc(perlocus_p(Ysnp, bio1, Z)))
rng = np.random.default_rng(7)
lp = [lambda_gc(perlocus_p(Ysnp, bio1[rng.permutation(n)], A)) for _ in range(10)]
print("    permuted bio1 (calibration):        %.3f" % np.mean(lp))
