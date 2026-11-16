"""Teosinte figures: Fig. 1 (climate gradient), Fig. 2 (TE resource),
Fig. 3 (QQ plot, SNP versus TE) and Fig. S1 (RDA framework)."""
from collections import Counter
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from common import (TEMatrix, SNPMatrix, climate, climate_pcs, teosinte_ancestry, env_vectors,
                    resid, rda1, perlocus_p, FIGURES)

plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False, "figure.dpi": 150})
clim = climate()
anc = teosinte_ancestry()

te = TEMatrix("teosinte")
bio1_t, lat_t, lon_t = env_vectors(te.pops, clim)
PC_t, ve = climate_pcs(te.pops, clim)
anc_t = np.array([anc[p] for p in te.pops])
Ytec = te.Y[:, te.clean]

snp = SNPMatrix("teosinte")
bio1_s, _, _ = env_vectors(snp.pops, clim)
PC_s, _ = climate_pcs(snp.pops, clim)
anc_s = np.array([anc[p] for p in snp.pops])

p_te = perlocus_p(Ytec, bio1_t, anc_t.reshape(-1, 1))
p_snp = perlocus_p(snp.Y, bio1_s, anc_s.reshape(-1, 1))
print("per-locus minimum p: SNP %.2g, TE %.2g" % (p_snp.min(), p_te.min()))

# Fig. 1: geography and climate PCA
fig, ax = plt.subplots(1, 2, figsize=(10, 4.2))
sc = ax[0].scatter(lon_t, lat_t, c=bio1_t, s=95, cmap="RdYlBu_r", edgecolor="k", lw=0.5)
ax[0].set_xlabel("Longitude (°)")
ax[0].set_ylabel("Latitude (°)")
ax[0].set_title("A", loc="left", fontweight="bold")
plt.colorbar(sc, ax=ax[0], label="bio1 (°C)")
sc2 = ax[1].scatter(PC_t[:, 0], PC_t[:, 1], c=bio1_t, s=95, cmap="RdYlBu_r", edgecolor="k", lw=0.5)
ax[1].set_xlabel("Climate PC1 (%d%%)" % round(ve[0] * 100))
ax[1].set_ylabel("Climate PC2 (%d%%)" % round(ve[1] * 100))
ax[1].set_title("B", loc="left", fontweight="bold")
plt.colorbar(sc2, ax=ax[1], label="bio1 (°C)")
plt.tight_layout()
plt.savefig(FIGURES / "Fig1_climate_gradient.png", bbox_inches="tight")
plt.close()

# Fig. 2: TE-insertion resource (order composition, frequency spectrum)
oco = Counter(te.order)
ks = ["LTR", "DNA", "MITE", "LINE"]
vs = [oco.get(k, 0) for k in ks]
fig, ax = plt.subplots(1, 2, figsize=(10, 4))
ax[0].bar(ks, vs, color=["#34495e", "#16a085", "#e67e22", "#8e44ad"], edgecolor="k")
for i, v in enumerate(vs):
    ax[0].text(i, v + 15, str(v), ha="center", fontsize=9)
ax[0].set_ylabel("Polymorphic insertions")
ax[0].set_title("A", loc="left", fontweight="bold")
af = Ytec.flatten()
af = af[af > 0]
ax[1].hist(af, bins=30, color="#34495e", edgecolor="w")
ax[1].set_xlabel("Per-population insertion frequency")
ax[1].set_ylabel("Count")
ax[1].set_title("B", loc="left", fontweight="bold")
plt.tight_layout()
plt.savefig(FIGURES / "Fig2_TE_landscape.png", bbox_inches="tight")
plt.close()


# Fig. 3: QQ plot of per-locus p-values
def qq(p):
    p = np.sort(p)
    m = len(p)
    return -np.log10((np.arange(1, m + 1) - 0.5) / m), -np.log10(np.clip(p, 1e-300, 1))


fig, ax = plt.subplots(figsize=(5, 5))
ex, oy = qq(p_snp)
ax.scatter(ex, oy, s=6, c="#c0392b", label="SNPs (n = %s)" % format(len(p_snp), ","), alpha=0.5)
ex, oy = qq(p_te)
ax.scatter(ex, oy, s=16, c="#7f8c8d", label="TE insertions (n = %s)" % format(len(p_te), ","), alpha=0.85)
lim = max(ax.get_xlim()[1], ax.get_ylim()[1])
ax.plot([0, lim], [0, lim], "k--", lw=1)
ax.set_xlabel("Expected $-\\log_{10}(p)$")
ax.set_ylabel("Observed $-\\log_{10}(p)$")
ax.legend(frameon=False)
plt.tight_layout()
plt.savefig(FIGURES / "Fig3_QQ_SNP_vs_TE.png", bbox_inches="tight")
plt.close()


# Fig. S1: RDA ordination of SNPs and significance of the three tests
def rda_sites(Y, X, Z):
    Yr, Xr = resid(Y, Z), resid(X, Z)
    B = np.linalg.lstsq(Xr, Yr, rcond=None)[0]
    U, S, _ = np.linalg.svd(Xr @ B, full_matrices=False)
    return U[:, :2] * S[:2]


A_s = anc_s.reshape(-1, 1)
sites = rda_sites(snp.Y, PC_s[:, :2], A_s)
# p-values as reported in the text (positive control, 02_teosinte_snp_rda.py) and
# in Table S2 (08_contrast_under_structure.py: SNP and TE rows re-ordered to the SNP matrix)
Yte_s = te.reorder(snp.pops)[:, te.clean]
tests = [("SNP~ancestry\n(pos.control)", rda1(snp.Y, anc_s, None, nperm=4999)[1]),
         ("SNP~climate", rda1(snp.Y, bio1_s, A_s, nperm=1999)[1]),
         ("TE~climate", rda1(Yte_s, bio1_s, A_s, nperm=1999)[1])]
print("RDA p-values shown in Fig. S1B:", [(t[0].split("\n")[0], round(t[1], 4)) for t in tests])
fig, ax = plt.subplots(1, 2, figsize=(10, 4.3))
sc = ax[0].scatter(sites[:, 0], sites[:, 1], c=bio1_s, s=115, cmap="RdYlBu_r", edgecolor="k", lw=0.5)
for i, p in enumerate(snp.pops):
    ax[0].annotate(p.replace("_", " "), (sites[i, 0], sites[i, 1]), textcoords="offset points", xytext=(0, 9), ha="center", fontsize=6, alpha=0.85)
ax[0].set_xlabel("RDA1")
ax[0].set_ylabel("RDA2")
ax[0].set_title("A", loc="left", fontweight="bold")
plt.colorbar(sc, ax=ax[0], label="bio1 (°C)")
labels = [t[0] for t in tests]
pv = [t[1] for t in tests]
ax[1].bar(labels, [-np.log10(x) for x in pv], color=["#2980b9", "#27ae60", "#7f8c8d"], edgecolor="k")
ax[1].axhline(-np.log10(0.05), ls="--", c="r", lw=1)
ax[1].set_ylabel("$-\\log_{10}(p)$ (RDA, bio1)")
ax[1].set_title("B", loc="left", fontweight="bold")
for i, x in enumerate(pv):
    ax[1].text(i, -np.log10(x) + 0.05, "p=%.4g" % x, ha="center", fontsize=8)
ax[1].text(2.05, -np.log10(0.05) + 0.05, "p=0.05", color="r", fontsize=8)
plt.tight_layout()
plt.savefig(FIGURES / "FigS1_RDA_contrast.png", bbox_inches="tight")
plt.close()
print("written Fig1, Fig2, Fig3, FigS1 to", FIGURES)
