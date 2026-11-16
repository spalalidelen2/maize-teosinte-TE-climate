"""Purifying selection: TE-insertion frequency by distance to the nearest gene.

Each insertion (B73 RefGen_V5 coordinates) is assigned the distance to the
nearest protein-coding gene from the Ensembl Plants release 57 annotation of
Zm-B73-REFERENCE-NAM-5.0 (downloaded on first use into data/annotation/).
Insertions are classed as genic, within 1 kb, 1-5 kb or more than 5 kb from a
gene, and the mean population frequency of each insertion is compared between
genic or near-gene insertions (within 1 kb) and insertions more than 5 kb away
with a one-sided Mann-Whitney test, overall and within LTR retrotransposons
alone.  Produces Figure 4.
"""
import gzip
import urllib.request
import numpy as np
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from common import TEMatrix, ANNOT, FIGURES

plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False, "figure.dpi": 150})
ANNOT.mkdir(parents=True, exist_ok=True)
GFF = ANNOT / "Zea_mays.Zm-B73-REFERENCE-NAM-5.0.57.gff3.gz"
GFF_URL = "https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-57/gff3/zea_mays/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.57.gff3.gz"
if not GFF.exists():
    print("downloading", GFF_URL)
    urllib.request.urlretrieve(GFF_URL, GFF)

# gene intervals per chromosome, merged where they overlap
genes = {}
with gzip.open(GFF, "rt") as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "gene" or not f[0].isdigit():
            continue
        genes.setdefault(f[0], []).append((int(f[3]), int(f[4])))
merged = {}
for c, iv in genes.items():
    iv.sort()
    m = []
    for s, e in iv:
        if m and s <= m[-1][1]:
            m[-1] = (m[-1][0], max(m[-1][1], e))
        else:
            m.append((s, e))
    merged[c] = (np.array([x[0] for x in m]), np.array([x[1] for x in m]))
print("genes on chr1-10: %d" % sum(len(v) for v in genes.values()))


def dist_to_gene(chrom, pos):
    c = chrom.replace("chr", "")
    if c not in merged:
        return 1e9
    st, en = merged[c]
    idx = np.searchsorted(st, pos, "right") - 1
    if idx >= 0 and en[idx] >= pos:
        return 0
    d = 1e18
    if idx >= 0:
        d = min(d, pos - en[idx])
    if idx + 1 < len(st):
        d = min(d, st[idx + 1] - pos)
    return d


classes = [("Genic", lambda d: d == 0), ("<1 kb", lambda d: (d > 0) & (d < 1000)),
           ("1-5 kb", lambda d: (d >= 1000) & (d < 5000)), (">5 kb", lambda d: d >= 5000)]
fig, axes = plt.subplots(1, 2, figsize=(8.4, 3.8), sharey=True)
for ax, (taxon, letter) in zip(axes, [("teosinte", "A"), ("maize", "B")]):
    te = TEMatrix(taxon)
    mf = te.Y.mean(0)                                   # mean population frequency per insertion
    d = np.array([dist_to_gene(c, p) for c, p in zip(te.chrom, te.pos)])
    print("\n== %s (%d insertions) ==" % (taxon, len(mf)))
    for name, sel in classes:
        m = sel(d)
        print("  %-7s n = %4d  mean frequency = %.3f" % (name, m.sum(), mf[m].mean()))
    near, far = mf[d < 1000], mf[d >= 5000]
    _, pu = stats.mannwhitneyu(near, far, alternative="less")
    print("  genic or within 1 kb (%.3f, n = %d) vs more than 5 kb (%.3f, n = %d): Mann-Whitney p = %.2g"
          % (near.mean(), len(near), far.mean(), len(far), pu))
    L = te.order == "LTR"
    _, pl = stats.mannwhitneyu(mf[L & (d < 1000)], mf[L & (d >= 5000)], alternative="less")
    print("  LTR retrotransposons only: %.3f (n = %d) vs %.3f (n = %d): p = %.2g"
          % (mf[L & (d < 1000)].mean(), (L & (d < 1000)).sum(), mf[L & (d >= 5000)].mean(), (L & (d >= 5000)).sum(), pl))
    for o in ["LTR", "DNA", "MITE", "LINE"]:
        om = te.order == o
        if om.sum() > 20:
            print("  %-4s within 1 kb of a gene: %.0f%% (n = %d)" % (o, 100 * (d[om] < 1000).mean(), om.sum()))
    c = te.clean
    _, pc = stats.mannwhitneyu(mf[c & (d < 1000)], mf[c & (d >= 5000)], alternative="less")
    print("  fully observed insertions only: %.3f (n = %d) vs %.3f (n = %d): p = %.2g"
          % (mf[c & (d < 1000)].mean(), (c & (d < 1000)).sum(), mf[c & (d >= 5000)].mean(), (c & (d >= 5000)).sum(), pc))
    data = [mf[sel(d)] for _, sel in classes]
    bp = ax.boxplot(data, labels=[n for n, _ in classes], showfliers=False, patch_artist=True, widths=0.6)
    for patch, col in zip(bp["boxes"], ["#e74c3c", "#e67e22", "#95a5a6", "#2c7fb8"]):
        patch.set_facecolor(col)
        patch.set_alpha(0.8)
    ax.set_title(letter, loc="left", fontweight="bold")
    ax.set_xlabel("Distance to nearest gene")
    ax.tick_params(axis="x", labelrotation=20)
axes[0].set_ylabel("Mean population insertion frequency")
plt.tight_layout()
plt.savefig(FIGURES / "Fig4_purifying_selection.png", bbox_inches="tight")
print("\nwritten", FIGURES / "Fig4_purifying_selection.png")
