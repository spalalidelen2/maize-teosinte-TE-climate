#!/usr/bin/env bash
# Create the project conda environment (maize_te) on the HPC cluster from
# environment.yml.  Run once on a login node:
#   bash scripts/teosinte/00_setup_env.sh
set -euo pipefail

if command -v module &>/dev/null; then
    module load mamba 2>/dev/null || module load anaconda 2>/dev/null || true
fi
if ! command -v conda &>/dev/null; then
    BASE_DIR="${WORK:-$HOME}/miniforge3"
    if [ ! -d "$BASE_DIR" ]; then
        wget -q https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh -O /tmp/miniforge.sh
        bash /tmp/miniforge.sh -b -p "$BASE_DIR"
    fi
    # shellcheck disable=SC1091
    source "$BASE_DIR/etc/profile.d/conda.sh"
fi
SOLVER=conda
command -v mamba &>/dev/null && SOLVER=mamba

# remove partially downloaded packages and limit parallelism (login-node limits)
conda clean --all -y 2>/dev/null || true
export CONDA_FETCH_THREADS=2 CONDA_VERIFY_THREADS=1 CONDA_EXECUTE_THREADS=1

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if conda env list | grep -q "maize_te"; then
    $SOLVER env update -n maize_te -f "$SCRIPT_DIR/environment.yml"
else
    $SOLVER env create -f "$SCRIPT_DIR/environment.yml"
fi

eval "$(conda shell.bash hook)"
conda activate maize_te
fastp --version 2>&1 | head -1 || true
bwa 2>&1 | grep -i version || true
samtools --version 2>&1 | head -1 || true
prefetch --version 2>&1 | head -1 || true
Rscript -e 'cat("R:", as.character(getRversion()), "| LEA:", requireNamespace("LEA", quietly=TRUE), "\n")' || true
