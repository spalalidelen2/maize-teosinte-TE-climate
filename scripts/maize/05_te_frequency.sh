#!/usr/bin/env bash
#SBATCH --job-name=maize_te_freq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=12:00:00
#SBATCH --output=logs/maize_te_freq_%j.log
# PoPoolationTE2 signature and frequency chain for the maize pools.  Same
# parameters as the teosinte chain (scripts/teosinte/09_te_frequency.sh)
# except the subsampling target coverage of 8 (maize pools are 7.9-16.1x
# versus at least 14.8x for teosinte).
#   sbatch scripts/maize/05_te_frequency.sh
# Output: $PROJ/maize/te_out/maize.teinsertions
set -uo pipefail
module load mamba 2>/dev/null || true
eval "$(conda shell.bash hook)"; conda activate maize_te
PROJ="${PROJ:-$PWD}"
cd "$PROJ"; mkdir -p maize/te_out
JAR=tools/popte2.jar; REF=reference/B73v5_plus_TE.fa; HIER=reference/te_hierarchy.txt
PP=maize/ppileup/maize.ppileup.gz; TC=8
java -Xmx16g -jar "$JAR" subsamplePpileup --ppileup "$PP" --target-coverage $TC --output maize/te_out/maize.ss.ppileup.gz || { echo "subsamplePpileup failed"; exit 1; }
java -Xmx16g -jar "$JAR" identifySignatures --ppileup maize/te_out/maize.ss.ppileup.gz --mode separate --min-count 2 --output maize/te_out/maize.signatures || { echo "identifySignatures failed"; exit 1; }
java -Xmx16g -jar "$JAR" frequency --ppileup maize/te_out/maize.ss.ppileup.gz --signature maize/te_out/maize.signatures --output maize/te_out/maize.freqsig || { echo "frequency failed"; exit 1; }
java -Xmx16g -jar "$JAR" filterSignatures --input maize/te_out/maize.freqsig --output maize/te_out/maize.filtered.freqsig --max-otherte-count 2 --max-structvar-count 2 || { echo "filterSignatures failed"; exit 1; }
java -Xmx16g -jar "$JAR" pairupSignatures --signature maize/te_out/maize.filtered.freqsig --ref-genome "$REF" --hier "$HIER" --min-distance -200 --max-distance 300 --output maize/te_out/maize.teinsertions || { echo "pairupSignatures failed"; exit 1; }
wc -l maize/te_out/maize.teinsertions
