#!/usr/bin/env bash
#SBATCH --job-name=te_freq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=12:00:00
#SBATCH --output=logs/te_freq_%j.log
# PoPoolationTE2 signature and frequency chain for the teosinte pools:
#   subsamplePpileup (target coverage 10) -> identifySignatures (mode separate,
#   min-count 2) -> frequency -> filterSignatures (max-otherte-count 2,
#   max-structvar-count 2) -> pairupSignatures (min-distance -200, max-distance 300)
#   sbatch scripts/teosinte/09_te_frequency.sh
# Output: $PROJ/te_out/hilo.teinsertions
set -euo pipefail
module load mamba 2>/dev/null || true
eval "$(conda shell.bash hook)"; conda activate maize_te
PROJ="${PROJ:-$PWD}"
cd "$PROJ"; mkdir -p te_out
JAR=tools/popte2.jar
REF=reference/B73v5_plus_TE.fa
HIER=reference/te_hierarchy.txt
PP=ppileup/hilo.ppileup.gz
TC=10
java -Xmx16g -jar "$JAR" subsamplePpileup --ppileup "$PP" --target-coverage $TC --output te_out/hilo.ss.ppileup.gz
java -Xmx16g -jar "$JAR" identifySignatures --ppileup te_out/hilo.ss.ppileup.gz --mode separate --min-count 2 --output te_out/hilo.signatures
java -Xmx16g -jar "$JAR" frequency --ppileup te_out/hilo.ss.ppileup.gz --signature te_out/hilo.signatures --output te_out/hilo.freqsig
java -Xmx16g -jar "$JAR" filterSignatures --input te_out/hilo.freqsig --output te_out/hilo.filtered.freqsig --max-otherte-count 2 --max-structvar-count 2
java -Xmx16g -jar "$JAR" pairupSignatures --signature te_out/hilo.filtered.freqsig --ref-genome "$REF" --hier "$HIER" --min-distance -200 --max-distance 300 --output te_out/hilo.teinsertions
wc -l te_out/hilo.teinsertions
