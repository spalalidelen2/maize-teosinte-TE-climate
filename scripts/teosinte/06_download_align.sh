#!/usr/bin/env bash
#SBATCH --job-name=dl_align
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=8:00:00
#SBATCH --array=1-245%25
#SBATCH --output=logs/dl_align_%A_%a.log
# SLURM array: one task per sequencing run of the teosinte sample sheet.
#   prefetch -> fasterq-dump -> fastp -> bwa mem -M (paired) -> sorted, indexed BAM
# FASTQ files are deleted after alignment.  Finished BAMs are skipped, so the
# array can be resubmitted to complete failed tasks.
#   mkdir -p logs; sbatch scripts/teosinte/06_download_align.sh
# Output: $PROJ/bam/<population>/<run>.bam
set -euo pipefail
module load mamba 2>/dev/null || true
eval "$(conda shell.bash hook)"; conda activate maize_te
PROJ="${PROJ:-$PWD}"
REF=$PROJ/reference/B73v5_plus_TE.fa
SHEET=$PROJ/metadata/sample_sheet_hilo_mexicana.tsv
THREADS=8

LINE=$(( SLURM_ARRAY_TASK_ID + 1 ))              # skip header
POP=$(sed -n "${LINE}p" "$SHEET" | cut -f1)
SRR=$(sed -n "${LINE}p" "$SHEET" | cut -f2)
[ -n "${SRR:-}" ] || { echo "no run on line $LINE"; exit 1; }
POPSAFE=${POP// /_}
BAMDIR=$PROJ/bam/$POPSAFE
mkdir -p "$BAMDIR"
OUT=$BAMDIR/${SRR}.bam
if [ -s "$OUT" ] && samtools quickcheck "$OUT" 2>/dev/null; then echo "$SRR already aligned"; exit 0; fi

TMP=$PROJ/fastq_tmp/$SRR
rm -rf "$TMP"; mkdir -p "$TMP"; cd "$TMP"
prefetch -O "$TMP" "$SRR"
fasterq-dump "$TMP/$SRR/$SRR.sra" --split-files -e $THREADS -t "$TMP" -O "$TMP"
fastp -i "${SRR}_1.fastq" -I "${SRR}_2.fastq" -o "${SRR}_1.tr.fastq" -O "${SRR}_2.tr.fastq" \
      -w $THREADS -j "${SRR}.fastp.json" -h "${SRR}.fastp.html"
bwa mem -M -t $THREADS -R "@RG\tID:${SRR}\tSM:${SRR}\tLB:${POPSAFE}\tPL:ILLUMINA" \
    "$REF" "${SRR}_1.tr.fastq" "${SRR}_2.tr.fastq" | samtools sort -@ $THREADS -o "$OUT" -
samtools index "$OUT"
cd "$PROJ"; rm -rf "$TMP"
echo "done $SRR"; ls -lh "$OUT"
