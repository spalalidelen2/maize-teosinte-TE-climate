#!/usr/bin/env bash
#SBATCH --job-name=maize_merge
#SBATCH --array=1-14
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=6:00:00
#SBATCH --output=logs/maize_merge_%A_%a.log
# Merge the per-individual maize BAMs into one pool per population (14 pools).
#   sbatch scripts/maize/03_merge_pools.sh
# Output: $PROJ/maize/pool/<population>.bam
set -uo pipefail
module load mamba 2>/dev/null || true
eval "$(conda shell.bash hook)"; conda activate maize_te
PROJ="${PROJ:-$PWD}"
cd "$PROJ"; mkdir -p maize/pool
POP=$(tail -n +2 metadata/sample_sheet_hilo_maize.tsv | cut -f1 | sort -u | sed -n "${SLURM_ARRAY_TASK_ID}p")
POPSAFE=${POP// /_}
N=$(ls maize/bam/"$POPSAFE"/*.bam 2>/dev/null | wc -l)
[ "$N" -gt 0 ] || { echo "no BAM files for $POP"; exit 1; }
echo "merging $N BAMs of $POP"
samtools merge -@ 4 -f "maize/pool/${POPSAFE}.bam" maize/bam/"$POPSAFE"/*.bam || { echo "merge failed"; exit 1; }
samtools index "maize/pool/${POPSAFE}.bam"
echo "done: $(samtools view -c maize/pool/${POPSAFE}.bam) reads"
