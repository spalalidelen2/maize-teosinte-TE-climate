"""Fig. S2: Manhattan plot of the per-SNP association with temperature in
teosinte (partial correlation controlling for maize ancestry)."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from common import SNPMatrix, climate, teosinte_ancestry, perlocus_p, bh_fdr, FIGURES

plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False, "figure.dpi": 150})
snp = SNPMatrix("teosinte")
bio1 = np.array([climate()[p]["bio1"] for p in snp.pops])
A = np.array([teosinte_ancestry()[p] for p in snp.pops]).reshape(-1, 1)
p = perlocus_p(snp.Y, bio1, A)
q = bh_fdr(p)
sig = q < 0.05
thr = p[sig].max() if sig.any() else 0
chrs = [str(c) for c in range(1, 11)]
off = 0
X = np.zeros(len(p))
ticks = []
for c in chrs:
    mk = snp.chrom == c
    if mk.sum() == 0:
        continue
    X[mk] = snp.pos[mk] + off
    ticks.append(off + snp.pos[mk].mean())
    off += snp.pos[mk].max() + 2e7
fig, ax = plt.subplots(figsize=(11, 3.6))
for i, c in enumerate(chrs):
    mk = snp.chrom == c
    ax.scatter(X[mk], -np.log10(p[mk]), s=4, c=("#34495e" if i % 2 == 0 else "#95a5a6"), rasterized=True)
ax.scatter(X[sig], -np.log10(p[sig]), s=12, c="#e74c3c", label="FDR<0.05 (n=%d)" % sig.sum(), zorder=3)
if thr > 0:
    ax.axhline(-np.log10(thr), ls="--", c="#e74c3c", lw=1)
ax.set_xticks(ticks)
ax.set_xticklabels(["chr%d" % i for i in range(1, 11)], fontsize=8)
ax.set_ylabel("$-\\log_{10}(p)$  SNP~temperature")
ax.set_xlabel("Genomic position (B73 RefGen_V4)")
ax.legend(frameon=False, loc="upper right")
plt.tight_layout()
plt.savefig(FIGURES / "FigS2_SNP_manhattan.png", bbox_inches="tight")
print("FDR < 0.05 SNPs: %d (p threshold %.2g); written FigS2" % (sig.sum(), thr))
