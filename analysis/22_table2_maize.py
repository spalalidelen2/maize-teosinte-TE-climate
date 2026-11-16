"""Table 2: sampling localities and summary statistics of the 14 sympatric
domesticated-maize (Zea mays ssp. mays) populations.

Columns as in Table 1; the final column is the mean genome-wide mexicana
introgression proportion (K3 admixture).  Elevation is the same SRTM 30 m
extraction as Table 1 (the maize and mexicana populations share localities).
"""
from collections import defaultdict
from common import read_tsv, climate, maize_introgression, elevation, META, TABLES, norm

h, rows = read_tsv(META / "sample_sheet_hilo_maize.tsv")
si = {c: i for i, c in enumerate(h)}
n, cov = defaultdict(int), defaultdict(float)
for r in rows:
    p = norm(r[si["population"]])
    n[p] += 1
    cov[p] += float(r[si["est_coverage"]])
clim = climate()
intro = maize_introgression()
elev = elevation()
pops = sorted(n, key=lambda p: (clim[p]["bio1"], p))   # ties broken by name
out = TABLES / "Table2_maize_population_metadata.tsv"
with open(out, "w") as o:
    o.write("Population\tLatitude\tLongitude\tElevation_m\tn\tPooled_cov_x\tbio1_MAT_C\tbio12_precip_mm\tMexicana_introgression\n")
    for p in pops:
        o.write("%s\t%.4f\t%.4f\t%d\t%d\t%.1f\t%.1f\t%.0f\t%.3f\n"
                % (p.replace("_", " "), clim[p]["lat"], clim[p]["lon"], round(elev[p]), n[p], cov[p],
                   clim[p]["bio1"], clim[p]["bio12"], intro[p][1]))
print("written", out)
import numpy as np
b = np.array([clim[p]["bio1"] for p in pops]); a = np.array([intro[p][1] for p in pops])
print("%d individuals, %d populations; pooled coverage %.1f-%.1f; bio1 %.1f-%.1f; introgression %.3f-%.3f; r(bio1, introgression) = %.2f"
      % (sum(n.values()), len(pops), min(cov.values()), max(cov.values()), b.min(), b.max(), a.min(), a.max(), np.corrcoef(b, a)[0, 1]))
