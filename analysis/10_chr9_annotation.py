"""Genomic context of the chr9:120-150 Mb candidate region.

Uses the B73 RefGen_V4 gene annotation (Ensembl Plants release 50, matching
the B73 RefGen_V4 coordinates of the SNP data) to compute gene density along
chromosome 9 in 10 Mb windows, the gene count and density of the candidate
region relative to the genome-wide mean, and compares the region with the
known inversion coordinates of the HILO project (Calfee et al. 2021).
The annotation and inversion table are downloaded on first use into
data/annotation/ (not tracked).
"""
import gzip
import urllib.request
import numpy as np
from common import ANNOT, CHR9_BLOCK

ANNOT.mkdir(parents=True, exist_ok=True)
GFF = ANNOT / "Zea_mays.B73_RefGen_v4.50.gff3.gz"
GFF_URL = "https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-50/gff3/zea_mays/Zea_mays.B73_RefGen_v4.50.gff3.gz"
INV = ANNOT / "knownInv_v4_coord.txt"
INV_URL = "https://raw.githubusercontent.com/ecalfee/hilo/master/data/refMaize/inversions/knownInv_v4_coord.txt"
for path, url in [(GFF, GFF_URL), (INV, INV_URL)]:
    if not path.exists():
        print("downloading", url)
        urllib.request.urlretrieve(url, path)

genes = []   # (chrom, start, end, has_description)
with gzip.open(GFF, "rt") as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "gene" or not f[0].isdigit():
            continue
        genes.append((f[0], int(f[3]), int(f[4]), "description=" in f[8]))
chrom = np.array([g[0] for g in genes])
start = np.array([g[1] for g in genes])
desc = np.array([g[3] for g in genes])
chr_len = {c: start[chrom == c].max() for c in sorted(set(chrom), key=int)}
genome_mb = sum(chr_len.values()) / 1e6
print("protein-coding genes on chr1-10: %d; genome-wide density %.1f genes/Mb" % (len(genes), len(genes) / genome_mb))

print("\nchr9 gene density in 10 Mb windows (genes/Mb):")
c9 = chrom == "9"
for lo in range(0, int(chr_len["9"] / 1e6) + 1, 10):
    n = ((start[c9] >= lo * 1e6) & (start[c9] < (lo + 10) * 1e6)).sum()
    print("  %3d-%3d Mb: %5.1f" % (lo, lo + 10, n / 10))

lo, hi = CHR9_BLOCK
reg = c9 & (start >= lo) & (start < hi)
print("\ncandidate region chr9:%d-%d Mb: %d genes, %.1f genes/Mb, %d with a functional description"
      % (lo / 1e6, hi / 1e6, reg.sum(), reg.sum() / ((hi - lo) / 1e6), (reg & desc).sum()))

print("\nknown inversions (HILO project, B73 RefGen_V4) overlapping the candidate region:")
with open(INV) as fh:
    for line in fh:
        name, c, s, e, ln = line.split()
        s, e = int(s), int(e)
        if c == "9" and s < hi and e > lo:
            ov = min(e, hi) - max(s, lo)
            print("  %s chr%s:%s-%s (%.1f Mb), %.0f%% of the inversion lies inside the candidate region"
                  % (name, c, format(s, ","), format(e, ","), (e - s) / 1e6, 100 * ov / (e - s)))
