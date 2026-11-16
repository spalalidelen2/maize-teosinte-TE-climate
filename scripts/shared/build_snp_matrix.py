"""Build the genome-wide SNP allele-frequency matrices from the published
per-population allele frequencies of Calfee et al. (2021)
(Figshare doi:10.6084/m9.figshare.16641799; ANGSD .mafs.gz files, B73 RefGen_V4).

For each taxon the sites present in every population are taken, expressed as
the non-reference allele frequency (see --polarise), and a random subset of
80,000 sites is written.  The maize matrix is restricted to sites polymorphic across
populations (variance > 0); the teosinte matrix was built without this filter
(--keep-monomorphic), so 1,742 of its 80,000 sites are monomorphic and are
removed by the variance filters inside the analyses.  To reproduce exactly the
subsets analysed in the study, pass the site lists distributed in data/
(chr, pos), in which case those sites are written instead of a random draw.

    python3 scripts/shared/build_snp_matrix.py teosinte --keep-monomorphic --sites data/snp_sites_teosinte.tsv
    python3 scripts/shared/build_snp_matrix.py maize --sites data/snp_sites_maize.tsv

Input : data/figshare/mafs/pop<N>.mafs.gz (downloaded on first use)
Output: data/snp_matrix_teosinte.tsv or data/snp_matrix_maize.tsv
        (columns: chr, pos, pop<N> ...)
"""
import argparse
import gzip
import json
import random
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
MAFS = ROOT / "data" / "figshare" / "mafs"
POPS = {"teosinte": [18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 33, 34, 35],
        "maize": [360, 361, 362, 363, 365, 366, 367, 368, 369, 370, 371, 372, 373, 374]}
FIGSHARE = "https://api.figshare.com/v2/articles/16641799"

ap = argparse.ArgumentParser()
ap.add_argument("taxon", choices=POPS)
ap.add_argument("--sites", help="site list (chr, pos) to write instead of a random subset")
ap.add_argument("--n", type=int, default=80000)
ap.add_argument("--keep-monomorphic", action="store_true",
                help="do not apply the among-population polymorphism filter (used for the teosinte matrix)")
ap.add_argument("--polarise", choices=["major", "minor"], default=None,
                help="rule for the non-reference frequency: 'major' = phat if the major allele is the "
                     "reference allele, else 1 - phat (teosinte); 'minor' = 1 - phat if the minor allele "
                     "is the reference allele, else phat (maize).  The two rules differ only at sites "
                     "where the reference allele is neither the major nor the minor allele.")
args = ap.parse_args()
polarise = args.polarise or ("major" if args.taxon == "teosinte" else "minor")
pops = POPS[args.taxon]
MAFS.mkdir(parents=True, exist_ok=True)

missing = [p for p in pops if not (MAFS / f"pop{p}.mafs.gz").exists()]
if missing:
    with urllib.request.urlopen(FIGSHARE, timeout=60) as fh:
        files = {f["name"]: f["download_url"] for f in json.load(fh)["files"]}
    for p in missing:
        name = f"pop{p}.mafs.gz"
        print("downloading", name)
        urllib.request.urlretrieve(files[name], MAFS / name)


def read_mafs(p):
    d = {}
    with gzip.open(MAFS / f"pop{p}.mafs.gz", "rt") as fh:
        h = fh.readline().rstrip().split("\t")
        ci = {c: i for i, c in enumerate(h)}
        for line in fh:
            c = line.rstrip().split("\t")
            phat = float(c[ci["phat"]])
            if polarise == "major":
                nonref = phat if c[ci["major"]] == c[ci["ref"]] else (1.0 - phat)
            else:
                nonref = (1.0 - phat) if c[ci["minor"]] == c[ci["ref"]] else phat
            d[(c[ci["chromo"]], c[ci["position"]])] = nonref
    return d


data = {p: read_mafs(p) for p in pops}
for p in pops:
    print("  pop%d: %d sites" % (p, len(data[p])))
common = set(data[pops[0]])
for p in pops[1:]:
    common &= set(data[p])
print("sites present in all %d populations: %d" % (len(pops), len(common)))

rows = []
for key in common:
    vals = [data[p][key] for p in pops]
    mn = sum(vals) / len(vals)
    if args.keep_monomorphic or sum((v - mn) ** 2 for v in vals) / len(vals) > 1e-9:
        rows.append((key[0], int(key[1]), vals))
print("sites retained: %d" % len(rows))
rows.sort(key=lambda r: (r[0], r[1]))          # deterministic order before sampling

if args.sites:
    want = set()
    with open(args.sites) as fh:
        next(fh)
        for line in fh:
            c, p = line.split()[:2]
            want.add((c, int(p)))
    rows = [r for r in rows if (r[0], r[1]) in want]
    print("sites from %s: %d" % (args.sites, len(rows)))
elif len(rows) > args.n:
    random.seed(42)
    rows = random.sample(rows, args.n)


def chrkey(r):
    try:
        return (int(r[0]), r[1])
    except ValueError:
        return (999, r[1])


rows.sort(key=chrkey)
out = ROOT / "data" / f"snp_matrix_{args.taxon}.tsv"
with open(out, "w") as o:
    o.write("chr\tpos\t" + "\t".join(f"pop{p}" for p in pops) + "\n")
    for chrom, pos, vals in rows:
        o.write(f"{chrom}\t{pos}\t" + "\t".join(f"{v:.4f}" for v in vals) + "\n")
print("written %s: %d SNPs x %d populations" % (out, len(rows), len(pops)))
