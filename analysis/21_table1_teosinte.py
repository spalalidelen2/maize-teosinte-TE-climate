"""Table 1: sampling localities and summary statistics of the 17 teosinte
(Zea mays ssp. mexicana) populations.

n and pooled coverage are taken from the sample sheet (one row per sequencing
run; pooled coverage is the sum of the per-run coverage estimates, i.e. the
depth of the merged population pool), elevation from the SRTM 30 m extraction,
climate from WorldClim v2.1 and maize ancestry from the Figshare ancestry
table.  Rows are ordered by increasing mean annual temperature.
"""
from collections import defaultdict
from common import read_tsv, climate, teosinte_ancestry, elevation, META, TABLES, norm

# individuals per population (one biosample per individual; some individuals have two runs)
h, rows = read_tsv(META / "sample_sheet_hilo_mexicana.tsv")
si = {c: i for i, c in enumerate(h)}
runs, cov, ind = defaultdict(int), defaultdict(float), defaultdict(set)
for r in rows:
    p = norm(r[si["population"]])
    runs[p] += 1
    cov[p] += float(r[si["coverageX"]])
    ind[p].add(r[si["biosample"]])
clim = climate()
anc = teosinte_ancestry()
elev = elevation()
pops = sorted(runs, key=lambda p: (clim[p]["bio1"], p))   # ties broken by name
out = TABLES / "Table1_population_metadata.tsv"
with open(out, "w") as o:
    o.write("Population\tLatitude\tLongitude\tElevation_m\tn\tPooled_cov_x\tbio1_MAT_C\tbio12_precip_mm\tMaize_ancestry\n")
    for p in pops:
        o.write("%s\t%.4f\t%.4f\t%d\t%d\t%.1f\t%.1f\t%.0f\t%.4f\n"
                % (p.replace("_", " "), clim[p]["lat"], clim[p]["lon"], round(elev[p]), len(ind[p]), cov[p],
                   clim[p]["bio1"], clim[p]["bio12"], anc[p]))
print("written", out)
print("%d individuals, %d populations, %d sequencing runs; elevation %d-%d m; pooled coverage %.1f-%.1f; bio1 %.1f-%.1f; ancestry %.3f-%.3f"
      % (sum(len(v) for v in ind.values()), len(pops), sum(runs.values()),
         min(elev[p] for p in pops), max(elev[p] for p in pops), min(cov.values()), max(cov.values()),
         min(clim[p]["bio1"] for p in pops), max(clim[p]["bio1"] for p in pops), min(anc.values()), max(anc.values())))
