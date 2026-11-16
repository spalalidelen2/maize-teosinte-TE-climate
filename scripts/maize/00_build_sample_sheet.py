"""Build the sympatric-maize sample sheet (129 runs, 14 populations).

Joins the ENA run table of BioProject PRJNA657016 (run accession, sample
alias, organism, FASTQ URLs) with the HILO sample metadata of Calfee et al.
(2021) (locality, Figshare population number, coordinates, elevation and
per-individual coverage), keeping the maize (Zea mays ssp. mays) individuals.
Requires internet access.
    python3 scripts/maize/00_build_sample_sheet.py [metadata_dir]
Output: <metadata_dir>/sample_sheet_hilo_maize.tsv
"""
import sys
import urllib.request
from collections import defaultdict
from pathlib import Path

metadir = Path(sys.argv[1] if len(sys.argv) > 1 else "data/metadata")
ENA = ("https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA657016&result=read_run"
       "&fields=run_accession,sample_alias,scientific_name,fastq_ftp&format=tsv")
HILO = "https://raw.githubusercontent.com/ecalfee/hilo/master/samples/HILO_meta.txt"

meta = {}
with urllib.request.urlopen(HILO, timeout=60) as fh:
    lines = fh.read().decode().splitlines()
h = lines[0].split("\t")
idx = {c: i for i, c in enumerate(h)}
for line in lines[1:]:
    c = line.split("\t")
    meta[c[idx["ID"]]] = c

rows = []
with urllib.request.urlopen(ENA, timeout=120) as fh:
    lines = fh.read().decode().splitlines()
eh = lines[0].split("\t")
ei = {c: i for i, c in enumerate(eh)}
for line in lines[1:]:
    c = line.split("\t")
    if "mexicana" in c[ei["scientific_name"]]:
        continue                                       # maize individuals only
    m = meta.get(c[ei["sample_alias"]])
    if not m:
        continue
    rows.append((m[idx["LOCALITY"]], m[idx["popN"]], c[ei["run_accession"]], c[ei["sample_alias"]],
                 m[idx["LAT"]], m[idx["LON"]], m[idx["ELEVATION"]], m[idx["est_coverage"]],
                 m[idx["zea"]], c[ei["fastq_ftp"]]))
out = metadir / "sample_sheet_hilo_maize.tsv"
with open(out, "w") as o:
    o.write("population\tpopN\trun\thilo_id\tlat\tlon\televation\test_coverage\tzea\tfastq_ftp\n")
    for r in rows:
        o.write("\t".join(r) + "\n")
pop = defaultdict(lambda: [0, 0.0])
for r in rows:
    pop[r[0]][0] += 1
    pop[r[0]][1] += float(r[7])
print("written %s: %d runs, %d populations" % (out, len(rows), len(pop)))
for p in sorted(pop):
    print("  %-20s runs = %3d  pooled coverage = %5.1f" % (p, pop[p][0], pop[p][1]))
