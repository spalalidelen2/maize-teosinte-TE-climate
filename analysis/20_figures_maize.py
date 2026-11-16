"""Maize supplementary figures: Fig. S3 (TE resource), Fig. S4 (RDA framework),
Fig. S5 (QQ plot) and Fig. S6 (Manhattan plot)."""
from collections import Counter
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from common import TEMatrix, SNPMatrix, climate, maize_introgression, rda1, perlocus_p, bh_fdr, lambda_gc, FIGURES

plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False, "figure.dpi": 150})
clim = climate()
intro = maize_introgression()
te = TEMatrix("maize")
snp = SNPMatrix("maize")
pops = te.pops                      # population order of the TE matrix (alphabetical)
n = len(pops)
Yte = te.Y[:, te.clean]
Ysnp = snp.reorder(pops)
bio1 = np.array([clim[p]["bio1"] for p in pops])
anc = np.array([intro[p][1] for p in pops])

# Fig. S3: maize TE resource
oc = Counter(te.order)
ords = ["LTR", "DNA", "MITE", "LINE"]
vals = [oc.get(o, 0) for o in ords]
mf = Yte.mean(0)
fig, (a1, a2) = plt.subplots(1, 2, figsize=(8, 3.4))
a1.bar(ords, vals, color=["#2c7fb8", "#7fcdbb", "#c7e9b4", "#edf8b1"], edgecolor="k", lw=0.5)
for i, v in enumerate(vals):
    a1.text(i, v + 30, str(v), ha="center", fontsize=9)
a1.set_ylabel("Polymorphic insertions")
a1.set_title("A", loc="left", fontweight="bold")
a2.hist(mf, bins=25, color="#2c7fb8", edgecolor="k", lw=0.4)
a2.set_xlabel("Mean population frequency")
a2.set_ylabel("Insertions")
a2.set_title("B", loc="left", fontweight="bold")
plt.tight_layout()
plt.savefig(FIGURES / "FigS3_maize_TE_landscape.png", bbox_inches="tight")
plt.close()

# Fig. S4: RDA framework (unconditioned)
tests = [("SNP~introg\n(pos.control)", rda1(Ysnp, anc, None)[1]),
         ("SNP~temp", rda1(Ysnp, bio1, None)[1]),
         ("TE~temp", rda1(Yte, bio1, None)[1])]
print("RDA p-values shown in Fig. S4B:", [(t[0].split("\n")[0], round(t[1], 4)) for t in tests])
Xs = Ysnp - Ysnp.mean(0)
xr = bio1 - bio1.mean()
B = (xr @ Xs) / (xr @ xr)
rda1_score = (np.outer(xr, B) @ B) / (np.linalg.norm(B) + 1e-9)
fig, (a1, a2) = plt.subplots(1, 2, figsize=(8.4, 3.6))
a1.scatter(rda1_score, bio1, c=bio1, cmap="coolwarm", s=80, edgecolor="k", lw=0.5)
a1.set_xlabel("RDA1 (SNP constrained by temperature)")
a1.set_ylabel("Mean annual temperature (C)")
a1.set_title("A", loc="left", fontweight="bold")
a2.bar(range(3), [-np.log10(t[1]) for t in tests], color=["#7f8c8d", "#c0392b", "#2c7fb8"], edgecolor="k", lw=0.5)
a2.axhline(-np.log10(0.05), ls="--", c="k", lw=0.8)
a2.set_xticks(range(3))
a2.set_xticklabels([t[0] for t in tests], fontsize=8)
a2.set_ylabel("$-\\log_{10}(p)$")
a2.set_title("B", loc="left", fontweight="bold")
for i, t in enumerate(tests):
    y = max(-np.log10(t[1]) + 0.05, -np.log10(0.05) + 0.08)   # keep labels clear of the threshold line
    a2.text(i, y, "p=%.4g" % t[1], ha="center", fontsize=8)
a2.text(-0.42, -np.log10(0.05) + 0.04, "p=0.05", fontsize=8, ha="left")
a2.set_ylim(0, a2.get_ylim()[1] * 1.08)
plt.tight_layout()
plt.savefig(FIGURES / "FigS4_maize_RDA_framework.png", bbox_inches="tight")
plt.close()

# Fig. S5: QQ plot (unconditioned per-locus correlation with temperature)
ps = perlocus_p(Ysnp, bio1, None)
pt = perlocus_p(Yte, bio1, None)
print("lambda_GC (unconditioned): SNP %.2f, TE %.2f" % (lambda_gc(ps), lambda_gc(pt)))


def qq(p):
    p = np.sort(p[~np.isnan(p)])
    m = len(p)
    return -np.log10((np.arange(1, m + 1) - 0.5) / m), -np.log10(p)


fig, ax = plt.subplots(figsize=(4.2, 4.2))
ex, ox = qq(ps)
ax.scatter(ex, ox, s=6, c="#c0392b", label="SNPs (n = %s)" % format(len(ps), ","), rasterized=True)
et, ot = qq(pt)
ax.scatter(et, ot, s=8, c="#7f8c8d", label="TE insertions (n = %s)" % format(len(pt), ","))
m = max(ex.max(), ox.max())
ax.plot([0, m], [0, m], "--", c="k", lw=0.8)
ax.set_xlabel("Expected $-\\log_{10}(p)$")
ax.set_ylabel("Observed $-\\log_{10}(p)$")
ax.legend(frameon=False, fontsize=8, loc="upper left")
plt.tight_layout()
plt.savefig(FIGURES / "FigS5_maize_QQ.png", bbox_inches="tight")
plt.close()

# Fig. S6: Manhattan plot
q = bh_fdr(ps)
sig = q < 0.05
thr = ps[sig].max() if sig.any() else 0
chrs = [str(c) for c in range(1, 11)]
off = 0
X = np.zeros(len(ps))
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
    ax.scatter(X[mk], -np.log10(ps[mk]), s=4, c=("#34495e" if i % 2 == 0 else "#95a5a6"), rasterized=True)
ax.scatter(X[sig], -np.log10(ps[sig]), s=12, c="#e74c3c", label="FDR<0.05 (n=%d)" % sig.sum(), zorder=3)
if thr > 0:
    ax.axhline(-np.log10(thr), ls="--", c="#e74c3c", lw=1)
ax.set_xticks(ticks)
ax.set_xticklabels(["chr%d" % i for i in range(1, 11)], fontsize=8)
ax.set_ylabel("$-\\log_{10}(p)$  maize SNP~temperature")
ax.set_xlabel("Genomic position (B73 RefGen_V4)")
ax.legend(frameon=False, loc="upper right")
plt.tight_layout()
plt.savefig(FIGURES / "FigS6_maize_manhattan.png", bbox_inches="tight")
plt.close()
print("maize FDR < 0.05 SNPs: %d; written FigS3, FigS4, FigS5, FigS6" % sig.sum())
