"""Elevation of each sampling locality from the SRTM 30 m digital elevation
model, queried through the Open Topo Data API (dataset srtm30m).

Writes data/metadata/elevation_srtm30m.tsv (population, lat, lon, elevation_m)
for the 17 localities of the climate table.  Run from the repository root:
    python3 scripts/shared/site_elevation.py
"""
import csv
import json
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
META = ROOT / "data" / "metadata"
with open(META / "climate_table_hilo.tsv") as fh:
    rows = list(csv.DictReader(fh, delimiter="\t"))
locs = "|".join("%s,%s" % (r["lat"], r["lon"]) for r in rows)
url = "https://api.opentopodata.org/v1/srtm30m?locations=" + locs
with urllib.request.urlopen(url, timeout=60) as resp:
    res = json.load(resp)["results"]
assert len(res) == len(rows)
out = META / "elevation_srtm30m.tsv"
with open(out, "w") as o:
    o.write("population\tlat\tlon\televation_m\n")
    for r, e in zip(rows, res):
        o.write("%s\t%s\t%s\t%.0f\n" % (r["population"], r["lat"], r["lon"], e["elevation"]))
        print("%-20s %9s %10s %6.0f m" % (r["population"], r["lat"], r["lon"], e["elevation"]))
print("written", out)
