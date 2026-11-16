#!/usr/bin/env bash
# Build the teosinte sample sheet from the metadata retrieved by 01_find_data.sh:
# runs of BioProject PRJNA657016 (HILO) with organism Zea mays subsp. mexicana
# and BioSample coordinates.  Per-run coverage is estimated as sequenced
# bases / 2.3 Gb (B73 genome size).
#   bash scripts/teosinte/02_build_sample_table.sh [metadata_dir]
# Output: <metadata_dir>/sample_sheet_hilo_mexicana.tsv
#   columns: population  run  biosample  lat  lon  coverageX
set -euo pipefail
MET="${1:-metadata}"
cd "$MET"
RUNINFO="teosinte_wgs_runinfo.csv"
COORDS="teosinte_biosample_WITHCOORDS.tsv"
OUT="sample_sheet_hilo_mexicana.tsv"
for f in "$RUNINFO" "$COORDS"; do
  [ -s "$f" ] || { echo "missing $MET/$f; run 01_find_data.sh first"; exit 1; }
done

# COORDS: $1 = BioSample, $3 = lat_lon ("LAT N LON W"), $4 = geo_loc_name ("Mexico:Locality")
# RUNINFO (csv): $1 = Run, $5 = bases, $22 = BioProject, $26 = BioSample, $29 = ScientificName
awk -F'\t' '
  FNR==NR { pop[$1]=$4; ll[$1]=$3; next }
  {
    n=split($0,a,",")
    run=a[1]; bases=a[5]; bp=a[22]; bs=a[26]; sci=a[29]
    if (bp=="PRJNA657016" && sci ~ /mexicana/ && bases ~ /^[0-9]+$/) {
      loc=ll[bs]; p=pop[bs]
      m=split(loc,c," ")
      if (m>=4 && p!="") {
        lat=c[1]; if(c[2]=="S") lat=-lat
        lon=c[3]; if(c[4]=="W") lon=-lon
        cov=bases/2.3e9
        sub(/^Mexico:/,"",p)
        printf "%s\t%s\t%s\t%s\t%s\t%.1f\n", p, run, bs, lat, lon, cov
      }
    }
  }
' "$COORDS" "$RUNINFO" | sort -k1,1 -k2,2 > "$OUT.tmp"
{ printf "population\trun\tbiosample\tlat\tlon\tcoverageX\n"; cat "$OUT.tmp"; } > "$OUT"
rm -f "$OUT.tmp"
echo "written $MET/$OUT: $(( $(wc -l < "$OUT") - 1 )) runs, $(tail -n +2 "$OUT" | cut -f1 | sort -u | wc -l) populations"
tail -n +2 "$OUT" | awk -F'\t' '{n[$1]++; c[$1]+=$6} END{for(p in n) printf "  %-22s runs=%2d  pooled coverage=%5.1fx\n", p, n[p], c[p]}' | sort
