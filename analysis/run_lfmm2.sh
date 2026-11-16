#!/bin/bash
# Run the LFMM2 analyses (R package LEA) for teosinte and maize.
# Intended for the HPC cluster where the project conda environment (R + LEA)
# is installed; run from the repository root on a compute node, e.g.
#   srun --nodes=1 --ntasks=1 --mem=8g --time=1:00:00 --pty bash
#   bash analysis/run_lfmm2.sh maize_te
# The SNP matrices (data/snp_matrix_*.tsv) must be present; see
# scripts/shared/build_snp_matrix.py.
set -euo pipefail
ENV="${1:-maize_te}"
module load mamba 2>/dev/null || true
ENVPATH="$(conda env list | awk -v e="$ENV" '$1==e{print $NF}')"
RS="${ENVPATH:+$ENVPATH/bin/}Rscript"
for f in data/snp_matrix_teosinte.tsv data/snp_matrix_maize.tsv \
         data/te_teosinte/te_freq_matrix.tsv data/te_maize/te_freq_matrix.tsv; do
  [[ -f "$f" ]] || { echo "missing input: $f"; exit 1; }
done
"$RS" -e 'if(!requireNamespace("LEA",quietly=TRUE)) stop("LEA is not installed in this environment"); cat("LEA", as.character(packageVersion("LEA")), "\n")'
"$RS" analysis/12_lfmm2_teosinte.R
"$RS" analysis/16_lfmm2_maize.R
