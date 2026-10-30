#!/usr/bin/env bash
# 03b_species_refs.sh - extract the per-gene sequences of the other complex-species
# reference genomes (config/species_refs.tsv) so 04_align.sh can put them under the kansasii REF.
#
# Usage: ./03b_species_refs.sh [-g GTDB_REPRESENTATIVES_DIR]   (genome of a species: <dir>/<accession>_query.fasta)
# Output: $AMR/work/species_refs/REFSP_<species>/genes/<tag>.fa  (header <tag>__REFSP_<species>)
# Already-extracted species are skipped; delete the directory to redo one.
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/env.sh"

GENOMES=$GTDB_REPRESENTATIVES
while getopts g: o; do case $o in g) GENOMES=$OPTARG ;; *) exit 1 ;; esac; done

while IFS=$'\t' read -r sp acc _; do
  [[ -z "$sp" || "$sp" == \#* ]] && continue
  out="$AMR/work/species_refs/REFSP_$sp"
  if [[ -d "$out/genes" ]]; then echo "[skip] REFSP_$sp"; continue; fi
  fa="$GENOMES/${acc}_query.fasta"
  [[ -s "$fa" ]] || { echo "missing genome for $sp ($acc): $fa" >&2; exit 1; }
  OUT_DIR="$out" "$HERE/03_sample_seqs.sh" asm "REFSP_$sp" "$fa"
done < "$CONFIG/species_refs.tsv"
