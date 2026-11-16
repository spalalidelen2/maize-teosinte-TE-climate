"""Noise-matched SNP control (maize).

Degrades the maize SNP frequencies to the measurement quality of the maize TE
calls (binomial resampling at 8-fold depth, the maize PoPoolationTE2
subsampling target, and a 0.20 detection floor) and re-runs the unconditioned
RDA on temperature over 20 replicates (Table S4).
"""
import numpy as np
from common import SNPMatrix, climate, rda1

snp = SNPMatrix("maize")
pops = snp.pops
bio1 = np.array([climate()[p]["bio1"] for p in pops])
Y0 = snp.Y
F, p = rda1(Y0, bio1, None, nperm=2999, seed=1)
print("baseline maize SNP ~ bio1 (unconditioned): F = %.2f  p = %.4f\n" % (F, p))


def run(label, depth, floor, reps=20):
    Fs, ps = [], []
    for rep in range(reps):
        rng = np.random.default_rng(1000 + rep)
        Fdeg = rng.binomial(depth, Y0) / depth
        if floor:
            Fdeg = np.where(Fdeg < 0.20, 0.0, Fdeg)
        keep = Fdeg.std(0) > 1e-9
        f, pv = rda1(Fdeg[:, keep], bio1, None, nperm=999, seed=2000 + rep)
        Fs.append(f)
        ps.append(pv)
    ps = np.array(ps)
    print("%-44s F = %.2f  p median = %.4f  (p < 0.05 in %d/%d replicates)" % (label, np.mean(Fs), np.median(ps), (ps < 0.05).sum(), reps))


print("degraded SNP frequencies, 20 replicates each:")
run("depth 8, no floor (sampling noise only):", 8, False)
run("depth 8 + 0.20 floor (full TE-like measurement):", 8, True)
