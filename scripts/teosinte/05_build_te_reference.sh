#!/usr/bin/env bash
#SBATCH --job-name=te_reference
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=6:00:00
#SBATCH --output=logs/te_reference_%j.log
# Build the TE-merged reference used by PoPoolationTE2:
#   1. mask the genomic TE copies of B73 RefGen_V5 with the structural TE annotation
#      (bedtools maskfasta);
#   2. append the curated Maize TE Consortium (MTEC) consensus library
#      distributed with EDTA (maizeTE02052020), after replacing non-ACGTN
#      characters by N (bwa cannot index other symbols);
#   3. derive the TE hierarchy file (id, family, order) from the library
#      sequence names, with family = order_superfamily so that a superfamily
#      shared by two orders (e.g. DTA under DNA and MITE) is not ambiguous;
#      the resulting file is distributed as data/te_hierarchy.txt;
#   4. index the merged reference with bwa.
# Inputs (in $PROJ/reference): B73v5.fa, B73 RefGen_V5 structural TE annotation
# (Zm-B73-REFERENCE-NAM-5.0.TE.gff3 from MaizeGDB), maizeTE02052020 (EDTA/database).
#   sbatch scripts/teosinte/05_build_te_reference.sh
set -euo pipefail
module load mamba 2>/dev/null || true
eval "$(conda shell.bash hook)"; conda activate maize_te
PROJ="${PROJ:-$PWD}"
REF="$PROJ/reference"
cd "$REF"
TEGFF="${TEGFF:-Zm-B73-REFERENCE-NAM-5.0.TE.gff3}"
TELIB="${TELIB:-maizeTE02052020}"

# 1. mask annotated TE copies in the genome
awk -F'\t' '!/^#/ {print $1"\t"($4-1)"\t"$5}' "$TEGFF" > B73v5_TE.bed
bedtools maskfasta -fi B73v5.fa -bed B73v5_TE.bed -fo B73v5.masked.fa

# 2. clean the TE consensus library (only ACGTN allowed by bwa)
awk '/^>/{print; next}{gsub(/[^ACGTNacgtn]/,"N"); print}' "$TELIB" > maizeTE02052020.clean.fa

# 3. TE hierarchy: id, family (order_superfamily), order.  Library headers are
#    of the form ">name#ORDER/SUPERFAMILY" (EDTA convention).
awk '/^>/ {
  id=substr($1,2); split(id,a,"#"); split(a[2],b,"/");
  ord=b[1]; sf=(b[2]==""?b[1]:b[2]);
  print id"\t"ord"_"sf"\t"ord
}' maizeTE02052020.clean.fa | sort -u > te_hierarchy.body
{ printf "id\tfamily\torder\n"; cat te_hierarchy.body; } > te_hierarchy.txt
rm -f te_hierarchy.body
echo "TE hierarchy: $(( $(wc -l < te_hierarchy.txt) - 1 )) sequences, $(tail -n +2 te_hierarchy.txt | cut -f2 | sort -u | wc -l) families, $(tail -n +2 te_hierarchy.txt | cut -f3 | sort -u | wc -l) orders"

# 4. merged reference and bwa index
cat B73v5.masked.fa maizeTE02052020.clean.fa > B73v5_plus_TE.fa
samtools faidx B73v5_plus_TE.fa
bwa index B73v5_plus_TE.fa
echo "done: $REF/B73v5_plus_TE.fa"
