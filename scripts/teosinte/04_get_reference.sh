#!/usr/bin/env bash
# Download the B73 RefGen_V5 assembly (Zm-B73-REFERENCE-NAM-5.0) from MaizeGDB
# and index it with samtools faidx.
#   bash scripts/teosinte/04_get_reference.sh
# PROJ: project directory on the cluster (override with the PROJ variable).
set -euo pipefail
PROJ="${PROJ:-$PWD}"
REF="$PROJ/reference"
mkdir -p "$REF" && cd "$REF"
URL="https://download.maizegdb.org/Zm-B73-REFERENCE-NAM-5.0/Zm-B73-REFERENCE-NAM-5.0.fa.gz"
curl -L -C - -o B73v5.fa.gz "$URL"
gunzip -c B73v5.fa.gz > B73v5.fa
samtools faidx B73v5.fa
echo "sequences: $(grep -c '^>' B73v5.fa)"
grep '^>' B73v5.fa | head -12
