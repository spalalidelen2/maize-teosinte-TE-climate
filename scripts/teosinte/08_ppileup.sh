#!/usr/bin/env bash
#SBATCH --job-name=ppileup
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=48:00:00
#SBATCH --output=logs/ppileup_%j.log
# PoPoolationTE2 ppileup across the 17 teosinte population pools
# (minimum mapping quality 15).  The pool order (ppileup/pop_order.txt) maps
# the sample index used by PoPoolationTE2 to the population BAM.
#   sbatch scripts/teosinte/08_ppileup.sh
set -euo pipefail
module load mamba 2>/dev/null || true
eval "$(conda shell.bash hook)"; conda activate maize_te
PROJ="${PROJ:-$PWD}"
cd "$PROJ"; mkdir -p ppileup
JAR=tools/popte2.jar
HIER=reference/te_hierarchy.txt
mapfile -t BAMS < <(ls pool/*.bam | sort)
printf '%s\n' "${BAMS[@]}" | nl > ppileup/pop_order.txt
ARGS=(); for b in "${BAMS[@]}"; do ARGS+=(--bam "$b"); done
java -Xmx16g -jar "$JAR" ppileup "${ARGS[@]}" --hier "$HIER" --map-qual 15 --output ppileup/hilo.ppileup.gz
ls -lh ppileup/hilo.ppileup.gz
