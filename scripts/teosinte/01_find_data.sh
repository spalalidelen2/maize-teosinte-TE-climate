#!/usr/bin/env bash
# Retrieve NCBI metadata for georeferenced teosinte whole-genome resequencing
# runs (SRA runinfo and BioSample coordinates) with Entrez Direct.
# Requires internet access (login node):
#   conda activate maize_te
#   bash scripts/teosinte/01_find_data.sh [metadata_dir]
# Output files are written to the metadata directory (default: ./metadata):
#   teosinte_wgs_runinfo.csv          all teosinte WGS runs (SRA runinfo)
#   teosinte_biosample_geo.tsv        BioSample accession, organism, lat_lon, locality, altitude
#   teosinte_biosample_WITHCOORDS.tsv the subset with coordinates
#   teosinte_wgs_georef_JOINED.tsv    WGS runs whose BioSample has coordinates
set -euo pipefail
OUT="${1:-metadata}"
mkdir -p "$OUT"
command -v esearch &>/dev/null || { echo "entrez-direct (esearch) not found; activate the maize_te environment"; exit 1; }

PARVI='"Zea mays subsp. parviglumis"[Organism]'
MEXI='"Zea mays subsp. mexicana"[Organism]'

# all teosinte WGS runs (genomic libraries only)
esearch -db sra -query "(${PARVI} OR ${MEXI}) AND \"wgs\"[Strategy] AND \"genomic\"[Source]" \
  | efetch -format runinfo > "$OUT/teosinte_wgs_runinfo.csv"
echo "teosinte WGS runs: $(( $(wc -l < "$OUT/teosinte_wgs_runinfo.csv") - 1 ))"

# BioSample coordinates (lat_lon), locality and altitude attributes
esearch -db biosample -query "(${PARVI} OR ${MEXI})" \
  | efetch -format docsum \
  | xtract -pattern DocumentSummary \
      -element Accession Organism \
      -block Attribute -if Attribute@harmonized_name -equals lat_lon -element Attribute \
      -block Attribute -if Attribute@harmonized_name -equals geo_loc_name -element Attribute \
      -block Attribute -if Attribute@harmonized_name -equals altitude -element Attribute \
  > "$OUT/teosinte_biosample_geo.tsv"
awk -F'\t' 'NF>=3 && $3 ~ /[0-9]/' "$OUT/teosinte_biosample_geo.tsv" > "$OUT/teosinte_biosample_WITHCOORDS.tsv"
echo "teosinte BioSamples: $(wc -l < "$OUT/teosinte_biosample_geo.tsv") | with coordinates: $(wc -l < "$OUT/teosinte_biosample_WITHCOORDS.tsv")"

# join: WGS runs whose BioSample has coordinates
cut -f1 "$OUT/teosinte_biosample_WITHCOORDS.tsv" | sort -u > "$OUT/_geo_biosamples.txt"
head -1 "$OUT/teosinte_wgs_runinfo.csv" > "$OUT/teosinte_wgs_georef_JOINED.tsv"
grep -F -f "$OUT/_geo_biosamples.txt" "$OUT/teosinte_wgs_runinfo.csv" >> "$OUT/teosinte_wgs_georef_JOINED.tsv" || true
rm -f "$OUT/_geo_biosamples.txt"
echo "georeferenced WGS runs: $(( $(wc -l < "$OUT/teosinte_wgs_georef_JOINED.tsv") - 1 ))"
