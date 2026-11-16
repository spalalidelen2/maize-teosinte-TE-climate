#!/usr/bin/env bash
#SBATCH --job-name=maize_ppileup
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=48:00:00
#SBATCH --output=logs/maize_ppileup_%j.log
# PoPoolationTE2 ppileup across the 14 maize pools (minimum mapping quality
# 15, as for teosinte).  Writes maize/ppileup/pop_order.txt (sample index to
# population BAM).
#   sbatch scripts/maize/04_ppileup.sh
set -uo pipefail
module load mamba 2>/dev/null || true
eval "$(conda shell.bash hook)"; conda activate maize_te
PROJ="${PROJ:-$PWD}"
cd "$PROJ"; mkdir -p maize/ppileup
JAR=tools/popte2.jar; HIER=reference/te_hierarchy.txt
mapfile -t BAMS < <(ls maize/pool/*.bam | sort)
printf '%s\n' "${BAMS[@]}" | nl > maize/ppileup/pop_order.txt
ARGS=(); for b in "${BAMS[@]}"; do ARGS+=(--bam "$b"); done
java -Xmx16g -jar "$JAR" ppileup "${ARGS[@]}" --hier "$HIER" --map-qual 15 --output maize/ppileup/maize.ppileup.gz || { echo "ppileup failed"; exit 1; }
ls -lh maize/ppileup/maize.ppileup.gz
