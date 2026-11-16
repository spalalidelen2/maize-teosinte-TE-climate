"""Per-population ancestry covariates from the published HILO ancestry
estimates of Calfee et al. (2021) (Figshare doi:10.6084/m9.figshare.16641799).

  data/metadata/maize_ancestry_perpop.tsv
      mean genome-wide maize ancestry of each mexicana population
      (ind_metadata_and_genomewide_ancestry.csv, mexicana individuals, by locality)
  data/metadata/mexicana_ancestry_maize_perpop.tsv
      mean genome-wide mexicana introgression (K3 admixture, alpha_mex) of each
      sympatric maize population (K3_alphas_by_symp_pop.txt), with the
      Figshare population number and bio1 of the locality.

    python3 scripts/shared/ancestry_covariates.py
The two Figshare files are downloaded into data/figshare/ on first use.
"""
import csv
import json
import urllib.request
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
FIG = ROOT / "data" / "figshare"
META = ROOT / "data" / "metadata"
FIG.mkdir(parents=True, exist_ok=True)
FILES = ["ind_metadata_and_genomewide_ancestry.csv", "K3_alphas_by_symp_pop.txt"]
if any(not (FIG / f).exists() for f in FILES):
    with urllib.request.urlopen("https://api.figshare.com/v2/articles/16641799", timeout=60) as fh:
        urls = {f["name"]: f["download_url"] for f in json.load(fh)["files"]}
    for f in FILES:
        if not (FIG / f).exists():
            urllib.request.urlretrieve(urls[f], FIG / f)


def norm(s):
    return s.replace(" ", "_")


with open(META / "climate_table_hilo.tsv") as fh:
    bio1 = {norm(r["population"]): float(r["bio1"]) for r in csv.DictReader(fh, delimiter="\t")}

# mexicana populations: mean maize ancestry per locality
acc = defaultdict(list)
with open(FIG / "ind_metadata_and_genomewide_ancestry.csv") as fh:
    rdr = csv.DictReader(fh)
    cols = rdr.fieldnames
    zea_col = [c for c in cols if c.lower() == "zea"][0]
    loc_col = [c for c in cols if c.upper() == "LOCALITY"][0]
    anc_col = [c for c in cols if "maize" in c.lower() and "ancestry" in c.lower()][0]
    for r in rdr:
        if r[zea_col] == "mexicana" and r[anc_col] not in ("", "NA"):
            acc[norm(r[loc_col])].append(float(r[anc_col]))
with open(META / "maize_ancestry_perpop.tsv", "w") as o:
    o.write("population\tmean_maize_ancestry\tbio1\n")
    for loc in sorted(acc):
        if loc in bio1:
            o.write("%s\t%.4f\t%s\n" % (loc, sum(acc[loc]) / len(acc[loc]), bio1[loc]))
print("mexicana populations with ancestry:", sum(1 for l in acc if l in bio1))

# maize populations: K3 alpha_mex by Figshare population number
loc_of = {}
with open(META / "sample_sheet_hilo_maize.tsv") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        loc_of[r["popN"]] = norm(r["population"])
amex = {}
with open(FIG / "K3_alphas_by_symp_pop.txt") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        amex[r["popN"]] = float(r["alpha_mex"])
with open(META / "mexicana_ancestry_maize_perpop.tsv", "w") as o:
    o.write("population\tpopN\tmean_mexicana_ancestry\tbio1\n")
    for popN, loc in sorted(loc_of.items(), key=lambda x: int(x[0])):
        if popN in amex:
            o.write("%s\t%s\t%.4f\t%s\n" % (loc, popN, amex[popN], bio1[loc]))
print("maize populations with introgression:", sum(1 for p in loc_of if p in amex))
