#!/usr/bin/env bash
#SBATCH --job-name=maize_align
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=8:00:00
#SBATCH --array=1-129%25
#SBATCH --output=logs/maize_align_%A_%a.log
# SLURM array: one task per maize sequencing run (same fastp and bwa mem -M
# settings as the teosinte pipeline, scripts/teosinte/06_download_align.sh).
# FASTQ files are taken from maize/fastq/ when present (01_download.sh) or
# downloaded from ENA; finished BAMs are skipped on resubmission.
#   mkdir -p logs; sbatch scripts/maize/02_align.sh
# Output: $PROJ/maize/bam/<population>/<run>.bam
set -uo pipefail
module load mamba 2>/dev/null || true
eval "$(conda shell.bash hook)"; conda activate maize_te
PROJ="${PROJ:-$PWD}"
REF=$PROJ/reference/B73v5_plus_TE.fa
SHEET=$PROJ/metadata/sample_sheet_hilo_maize.tsv
THREADS=8
LINE=$(( SLURM_ARRAY_TASK_ID + 1 ))
POP=$(sed -n "${LINE}p" "$SHEET" | cut -f1)
SRR=$(sed -n "${LINE}p" "$SHEET" | cut -f3)
FTP=$(sed -n "${LINE}p" "$SHEET" | cut -f10)
[ -n "${SRR:-}" ] || { echo "no run on line $LINE"; exit 1; }
POPSAFE=${POP// /_}
BAMDIR=$PROJ/maize/bam/$POPSAFE; mkdir -p "$BAMDIR"; OUT=$BAMDIR/${SRR}.bam
if [ -s "$OUT" ] && samtools quickcheck "$OUT" 2>/dev/null; then echo "$SRR already aligned"; exit 0; fi
TMP=$PROJ/maize/tmp/$SRR; rm -rf "$TMP"; mkdir -p "$TMP"; cd "$TMP"
FQ=$PROJ/maize/fastq
IFS=';' read -r U1 U2 <<< "$FTP"
get() {   # $1 = ENA url, $2 = destination
  if [ -s "$FQ/$(basename "$2")" ]; then cp "$FQ/$(basename "$2")" "$2"; return 0; fi
  wget -q -O "$2" "https://$1"
}
get "$U1" "${SRR}_1.fastq.gz" || { echo "R1 download failed"; exit 1; }
get "$U2" "${SRR}_2.fastq.gz" || { echo "R2 download failed"; exit 1; }
fastp -i "${SRR}_1.fastq.gz" -I "${SRR}_2.fastq.gz" -o "${SRR}_1.tr.fastq" -O "${SRR}_2.tr.fastq" \
      -w $THREADS -j "${SRR}.fastp.json" -h "${SRR}.fastp.html" || { echo "fastp failed"; exit 1; }
bwa mem -M -t $THREADS -R "@RG\tID:${SRR}\tSM:${SRR}\tLB:${POPSAFE}\tPL:ILLUMINA" \
    "$REF" "${SRR}_1.tr.fastq" "${SRR}_2.tr.fastq" | samtools sort -@ $THREADS -o "$OUT" - || { echo "bwa/sort failed"; exit 1; }
samtools index "$OUT"
cd "$PROJ"; rm -rf "$TMP"
echo "done $SRR"; ls -lh "$OUT"
