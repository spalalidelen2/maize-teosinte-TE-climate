"""Shared data loaders and statistics for the genotype-environment association analyses.

All analysis scripts import this module and resolve data paths relative to the
repository root, so they can be run from any working directory:

    python3 analysis/01_teosinte_te_rda.py
"""
from pathlib import Path
import csv
import numpy as np
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
META = DATA / "metadata"
TABLES = DATA / "tables"
FIGURES = ROOT / "figures"
ANNOT = DATA / "annotation"          # downloaded gene annotations (not tracked)

TE_META_COLS = {"locus", "chr", "family", "order", "pos", "nsig", "nDetected", "variance"}

# Figshare population numbers (Calfee et al. 2021) of the 17 mexicana populations.
TEOSINTE_POPN = {18: "Nabogame", 19: "Penjamillo", 20: "Malinalco", 21: "El_Porvenir",
                 22: "Amecameca", 23: "Tlapala", 24: "Jicaltepec", 25: "Cocotilan",
                 26: "Opopeo", 27: "Amatlan", 28: "Puruandiro", 29: "Ixtlan",
                 30: "San_Pedro", 31: "Tenango_del_Aire", 33: "Puerta_Encantada",
                 34: "Santa_Clara", 35: "Xochimilco"}

CHR9_BLOCK = (120e6, 150e6)          # chr9 candidate region (B73 RefGen_V4 coordinates)


def read_tsv(path):
    with open(path) as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
    return rows[0], rows[1:]


def norm(name):
    return name.replace(" ", "_")


# ---------------------------------------------------------------- environment
def climate():
    """Locality -> dict of lat, lon and bio1..bio19 (WorldClim v2.1)."""
    h, rows = read_tsv(META / "climate_table_hilo.tsv")
    ci = {c: i for i, c in enumerate(h)}
    out = {}
    for r in rows:
        d = {"lat": float(r[ci["lat"]]), "lon": float(r[ci["lon"]])}
        for k in range(1, 20):
            d["bio%d" % k] = float(r[ci["bio%d" % k]])
        out[norm(r[ci["population"]])] = d
    return out


def climate_pcs(pops, clim=None):
    """PCA of the 19 standardised bioclimatic variables for the given populations.

    Returns (scores, variance_explained) with scores aligned to `pops`.
    """
    clim = clim or climate()
    bio = np.array([[clim[p]["bio%d" % k] for k in range(1, 20)] for p in pops])
    bz = (bio - bio.mean(0)) / bio.std(0)
    U, S, Vt = np.linalg.svd(bz - bz.mean(0), full_matrices=False)
    return bz @ Vt.T, S ** 2 / (S ** 2).sum()


def teosinte_ancestry():
    """Locality -> mean genome-wide maize ancestry of the mexicana population."""
    h, rows = read_tsv(META / "maize_ancestry_perpop.tsv")
    ai = {c: i for i, c in enumerate(h)}
    return {norm(r[ai["population"]]): float(r[ai["mean_maize_ancestry"]]) for r in rows}


def maize_introgression():
    """Locality -> (popN, mean genome-wide mexicana introgression) of the maize population."""
    h, rows = read_tsv(META / "mexicana_ancestry_maize_perpop.tsv")
    ai = {c: i for i, c in enumerate(h)}
    return {norm(r[ai["population"]]): (r[ai["popN"]], float(r[ai["mean_mexicana_ancestry"]]))
            for r in rows}


def elevation():
    """Locality -> SRTM 30 m elevation (m) at the sampling coordinates."""
    h, rows = read_tsv(META / "elevation_srtm30m.tsv")
    ei = {c: i for i, c in enumerate(h)}
    return {norm(r[ei["population"]]): float(r[ei["elevation_m"]]) for r in rows}


# ---------------------------------------------------------------- TE matrices
class TEMatrix:
    """Population-level TE-insertion frequency matrix (populations x loci)."""

    def __init__(self, taxon):
        path = DATA / ("te_teosinte" if taxon == "teosinte" else "te_maize") / "te_freq_matrix.tsv"
        h, rows = read_tsv(path)
        ci = {c: i for i, c in enumerate(h)}
        self.pops = [c for c in h if c not in TE_META_COLS]
        self.Y = np.array([[float(r[ci[p]]) for p in self.pops] for r in rows]).T
        self.nDetected = np.array([int(r[ci["nDetected"]]) for r in rows])
        self.order = np.array([r[ci["order"]] for r in rows])
        self.family = np.array([r[ci["family"]] for r in rows])
        self.chrom = np.array([r[ci["chr"]] for r in rows])
        self.pos = np.array([int(r[ci["pos"]]) for r in rows])
        self.npop = len(self.pops)

    @property
    def clean(self):
        """Loci with a measured (non-zero) frequency in every population."""
        return self.nDetected == self.npop

    def reorder(self, pops):
        """Rows re-ordered to the given population order."""
        return self.Y[[self.pops.index(p) for p in pops]]


# --------------------------------------------------------------- SNP matrices
class SNPMatrix:
    """Genome-wide SNP non-reference allele-frequency matrix (populations x SNPs)."""

    def __init__(self, taxon):
        if taxon == "teosinte":
            path = DATA / "snp_matrix_teosinte.tsv"
            popn2loc = {str(k): v for k, v in TEOSINTE_POPN.items()}
        else:
            path = DATA / "snp_matrix_maize.tsv"
            popn2loc = {v[0]: k for k, v in maize_introgression().items()}
        with open(path) as fh:
            h = fh.readline().rstrip("\n").split("\t")
            rows = [l.rstrip("\n").split("\t") for l in fh]
        self.pops = [popn2loc[c[3:]] for c in h[2:]]          # locality per column
        self.chrom = np.array([r[0] for r in rows])
        self.pos = np.array([int(r[1]) for r in rows])
        Y = np.array([[float(r[i]) for i in range(2, len(h))] for r in rows]).T
        # the published frequencies carry a small numerical pad outside [0, 1]
        self.Y = np.clip(Y, 0.0, 1.0)
        self.npop = len(self.pops)

    def reorder(self, pops):
        """Rows re-ordered to the given population order."""
        return self.Y[[self.pops.index(p) for p in pops]]

    @property
    def chr9_block(self):
        lo, hi = CHR9_BLOCK
        return (self.chrom == "9") & (self.pos >= lo) & (self.pos < hi)

    def background_pcs(self, k=2, exclude=None, pops=None):
        """Leading population PCs of the SNP matrix, excluding the chr9 candidate block.
        Rows follow `pops` when given, otherwise the matrix column order."""
        mask = ~self.chr9_block if exclude is None else ~exclude
        Y = self.Y if pops is None else self.reorder(pops)
        X = Y[:, mask]
        Xs = (X - X.mean(0)) / (X.std(0) + 1e-9)
        U, S, Vt = np.linalg.svd(Xs, full_matrices=False)
        return U[:, :k]


def env_vectors(pops, clim=None):
    """bio1, latitude and longitude aligned to a population order."""
    clim = clim or climate()
    return (np.array([clim[p]["bio1"] for p in pops]),
            np.array([clim[p]["lat"] for p in pops]),
            np.array([clim[p]["lon"] for p in pops]))


# ------------------------------------------------------------------ statistics
def resid(Y, Z):
    """Residuals of the columns of Y after regression on Z (with intercept)."""
    if Z is None or Z.shape[1] == 0:
        Z1 = np.ones((len(Y), 1))
    else:
        Z1 = np.column_stack([np.ones(len(Z)), Z])
    return Y - Z1 @ np.linalg.lstsq(Z1, Y, rcond=None)[0]


def rda(Y, X, Z=None, nperm=4999, seed=1):
    """Partial redundancy analysis with a permutation test (general form).

    Y: response matrix (populations x loci); X: predictors (populations x p);
    Z: conditioning covariates (populations x q) or None.  Significance is
    assessed by permuting the rows of the residual response matrix.
    Returns (constrained R2, pseudo-F, p).
    """
    n = Y.shape[0]
    if Z is None or Z.shape[1] == 0:
        Yr, Xr, q = Y - Y.mean(0), X - X.mean(0), 0
    else:
        Yr, Xr, q = resid(Y, Z), resid(X, Z), Z.shape[1]
    p = Xr.shape[1]
    Q = np.linalg.qr(Xr)[0]            # orthonormal basis of the predictor space

    def css(Ym):
        # constrained sum of squares = squared norm of the projection of Ym on Q
        return np.sum((Q.T @ Ym) ** 2), np.sum(Ym ** 2)

    cs, ts = css(Yr)
    F = (cs / p) / ((ts - cs) / (n - p - q - 1))
    rng = np.random.default_rng(seed)
    count = 1
    for _ in range(nperm):
        c2, t2 = css(Yr[rng.permutation(n)])
        if (c2 / p) / ((t2 - c2) / (n - p - q - 1)) >= F - 1e-12:
            count += 1
    return cs / ts, F, count / (nperm + 1)


def rda1(Y, x, Z=None, nperm=1999, seed=1):
    """Partial RDA with a single predictor, vectorised permutation test.

    Permuting the rows of the residual response matrix is equivalent to
    permuting the residual predictor, which allows a fast implementation.
    Returns (pseudo-F, p).
    """
    n = Y.shape[0]
    q = 0 if Z is None else Z.shape[1]
    Yr = resid(Y, Z)
    xr = resid(x.reshape(-1, 1), Z).ravel()
    xn = (xr * xr).sum()
    ts = (Yr * Yr).sum()
    df = n - q - 2

    def Fs(xv):
        g = Yr.T @ xv
        cs = (g * g).sum() / xn
        return cs / ((ts - cs) / df)

    F = Fs(xr)
    rng = np.random.default_rng(seed)
    count = 1
    for _ in range(nperm):
        if Fs(xr[rng.permutation(n)]) >= F - 1e-12:
            count += 1
    return F, count / (nperm + 1)


def perlocus_p(Y, x, Z=None):
    """Two-sided p-value of the per-locus partial correlation between each column
    of Y and x, controlling for the covariates Z."""
    n = Y.shape[0]
    q = 0 if Z is None else Z.shape[1]
    df = n - q - 2
    Yr = resid(Y, Z)
    xr = resid(x.reshape(-1, 1), Z).ravel()
    sd = Yr.std(0)
    r = np.where(sd > 1e-9,
                 (Yr * xr[:, None]).sum(0) / (np.sqrt((Yr ** 2).sum(0)) * np.sqrt((xr * xr).sum()) + 1e-12),
                 0)
    t = r * np.sqrt(df / np.clip(1 - r * r, 1e-9, 1))
    return 2 * stats.t.sf(np.abs(t), df)


def bh_fdr(p):
    """Benjamini-Hochberg adjusted p-values."""
    p = np.asarray(p)
    o = np.argsort(p)
    m = len(p)
    q = np.empty(m)
    prev = 1.0
    for i, idx in enumerate(o[::-1]):
        rank = m - i
        prev = min(prev, p[idx] * m / rank)
        q[idx] = prev
    return q


def lambda_gc(p):
    """Genomic inflation factor from a vector of p-values (1 df chi-square)."""
    chi = stats.chi2.isf(np.clip(p, 1e-300, 1), 1)
    return np.median(chi) / 0.4549364
