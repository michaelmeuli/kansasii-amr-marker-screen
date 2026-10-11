#!/usr/bin/env bash
# 03b_species_refs.sh - extract the per-gene sequences of the other complex-species
# reference genomes (species_refs.tsv) so 04_align.sh can put them under the kansasii REF.
#
# Usage: ./03b_species_refs.sh
# Output: $ARM/work/species_refs/REFSP_<species>/genes/<tag>.fa  (header <tag>__REFSP_<species>)
# Already-extracted species are skipped; delete the directory to redo one.
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/env.sh"

while IFS=$'\t' read -r sp acc path; do
  [[ -z "$sp" || "$sp" == \#* ]] && continue
  out="$ARM/work/species_refs/REFSP_$sp"
  if [[ -d "$out/genes" ]]; then echo "[skip] REFSP_$sp"; continue; fi
  fa="$ROOT/$path"
  [[ -s "$fa" ]] || { echo "missing genome for $sp ($acc): $fa" >&2; exit 1; }
  OUT_DIR="$out" "$HERE/03_sample_seqs.sh" asm "REFSP_$sp" "$fa"
done < "$HERE/species_refs.tsv"
