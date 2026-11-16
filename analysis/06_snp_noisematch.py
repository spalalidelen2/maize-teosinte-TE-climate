"""Noise-matched SNP control (teosinte).

Degrades the SNP frequencies to the measurement quality of the pooled TE
calls (binomial resampling at 10-fold depth, the PoPoolationTE2 subsampling
target, and a 0.20 detection floor) and re-runs the partial RDA
(bio1 | ancestry) and the per-locus FDR test over 20 replicates.
"""
import numpy as np
from common import SNPMatrix, climate, teosinte_ancestry, rda1, perlocus_p, bh_fdr, read_tsv, TABLES

snp = SNPMatrix("teosinte")
pops = snp.pops
n = len(pops)
clim = climate()
bio1 = np.array([clim[p]["bio1"] for p in pops])
anc = teosinte_ancestry()
A = np.array([anc[p] for p in pops]).reshape(-1, 1)
h, rows = read_tsv(TABLES / "Table1_population_metadata.tsv")
ti = {c: i for i, c in enumerate(h)}
cov = np.array([{r[ti["Population"]].replace(" ", "_"): float(r[ti["Pooled_cov_x"]]) for r in rows}[p] for p in pops])
Y0 = snp.Y
print("populations = %d, SNPs = %d, pooled coverage %.1f-%.1f" % (n, Y0.shape[1], cov.min(), cov.max()))


def fdr_count(Y):
    return (bh_fdr(perlocus_p(Y, bio1, A)) < 0.05).sum()


F, p = rda1(Y0, bio1, A, nperm=2999, seed=1)
print("\nBaseline (undegraded SNPs): F = %.2f  p = %.4f  FDR < 0.05: %d\n" % (F, p, fdr_count(Y0)))


def run(label, depth_mode, floor, reps=20):
    Fs, ps, fdrs, nkeep = [], [], [], []
    for rep in range(reps):
        rng = np.random.default_rng(1000 + rep)
        d = np.full(n, 10) if depth_mode == "fixed10" else np.maximum(np.round(cov).astype(int), 1)
        Fdeg = rng.binomial(d[:, None], Y0) / d[:, None]
        if floor:
            Fdeg = np.where(Fdeg < 0.20, 0.0, Fdeg)
        keep = Fdeg.std(0) > 1e-9
        Yk = Fdeg[:, keep]
        nkeep.append(keep.sum())
        F, p = rda1(Yk, bio1, A, nperm=999, seed=2000 + rep)
        Fs.append(F)
        ps.append(p)
        fdrs.append(fdr_count(Yk))
    ps = np.array(ps)
    print("%-48s F = %.2f  p median = %.4f  mean = %.4f  (p < 0.05 in %d/%d replicates)"
          % (label, np.mean(Fs), np.median(ps), ps.mean(), (ps < 0.05).sum(), reps))
    print("%-48s FDR < 0.05 SNPs: median %d (range %d-%d); polymorphic SNPs retained about %d"
          % ("", int(np.median(fdrs)), min(fdrs), max(fdrs), int(np.median(nkeep))))


print("Degraded SNP frequencies, 20 replicates each:")
run("depth 10, no floor (sampling noise only):", "fixed10", False)
run("actual pooled coverage + 0.20 floor (floor only):", "coverage", True)
run("depth 10 + 0.20 floor (full TE-like measurement):", "fixed10", True)
