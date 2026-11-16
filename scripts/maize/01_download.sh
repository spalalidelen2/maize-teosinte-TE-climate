#!/bin/bash
# Download the FASTQ files of the 129 sympatric-maize runs (PRJNA657016) from
# ENA, using the fastq_ftp links in the maize sample sheet.  Optional: the
# alignment script (02_align.sh) downloads any run that is not yet present.
#   bash scripts/maize/01_download.sh
# Reads : metadata/sample_sheet_hilo_maize.tsv
# Writes: maize/fastq/<run>_1.fastq.gz, <run>_2.fastq.gz
set -uo pipefail
PROJ="${PROJ:-$PWD}"
cd "$PROJ"
SHEET="metadata/sample_sheet_hilo_maize.tsv"
OUT="maize/fastq"
mkdir -p "$OUT"
[[ -f "$SHEET" ]] || { echo "missing $SHEET"; exit 1; }
n=0; ok=0; fail=0
while IFS=$'\t' read -r pop popN run hilo lat lon elev cov zea ftp; do
  [[ "$pop" == "population" || -z "$run" ]] && continue
  n=$((n+1))
  IFS=';' read -ra urls <<< "$ftp"
  for u in "${urls[@]}"; do
    [[ -z "$u" ]] && continue
    fn="$OUT/$(basename "$u")"
    [[ -s "$fn" ]] && continue
    if wget -q -O "$fn" "https://$u"; then ok=$((ok+1)); else echo "failed: $u"; rm -f "$fn"; fail=$((fail+1)); fi
  done
done < "$SHEET"
echo "runs: $n | files downloaded: $ok | failed: $fail | fastq files present: $(ls "$OUT"/*.fastq.gz 2>/dev/null | wc -l)"
