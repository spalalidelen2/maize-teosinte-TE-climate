#!/usr/bin/env bash
#SBATCH --job-name=merge_pools
#SBATCH --array=1-17
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=6:00:00
#SBATCH --output=logs/merge_%A_%a.log
# Merge the per-individual BAMs of each teosinte population into one
# population pool (17 pools).
#   sbatch scripts/teosinte/07_merge_pools.sh
# Output: $PROJ/pool/<population>.bam
set -euo pipefail
module load mamba 2>/dev/null || true
eval "$(conda shell.bash hook)"; conda activate maize_te
PROJ="${PROJ:-$PWD}"
cd "$PROJ"; mkdir -p pool
POP=$(tail -n +2 metadata/sample_sheet_hilo_mexicana.tsv | cut -f1 | sort -u | sed -n "${SLURM_ARRAY_TASK_ID}p")
POPSAFE=${POP// /_}
N=$(ls bam/"$POPSAFE"/*.bam 2>/dev/null | wc -l)
[ "$N" -gt 0 ] || { echo "no BAM files for $POP"; exit 1; }
echo "merging $N BAMs of $POP"
samtools merge -@ 4 -f "pool/${POPSAFE}.bam" bam/"$POPSAFE"/*.bam
samtools index "pool/${POPSAFE}.bam"
echo "done: $(samtools view -c pool/${POPSAFE}.bam) reads"
