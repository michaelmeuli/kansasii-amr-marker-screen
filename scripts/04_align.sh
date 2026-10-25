#!/usr/bin/env bash
# 04_align.sh - build one MAFFT alignment per gene: reference + the other complex-species
# references (03b_species_refs.sh) + every sample.
#
# Usage: ./04_align.sh [THREADS]
#
# For coding genes this also produces a codon-aware alignment (translate -> align
# protein -> back-translate), which is what you want before calling codon numbers:
# a plain nucleotide alignment can open gaps out of frame and shift every
# downstream codon by one or two.
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/env.sh"

THREADS=${1:-4}
ALN=$ARM/work/aln
mkdir -p "$ALN"

# Which genes are protein-coding, taken from targets.tsv feature type
mapfile -t CODING < <(awk -F'\t' '!/^#/ && $5 ~ /CDS/ {print $1}' "$CONFIG/targets.tsv")
is_coding () { for g in "${CODING[@]}"; do [[ "${1%%_c[0-9]*}" == "$g" ]] && return 0; done; return 1; }

for refgene in "$ARM"/work/ref_genes/*.fa; do
  tag=$(basename "$refgene" .fa)
  in="$ALN/$tag.input.fa"
  # kansasii REF: >REF_<tag>|locus|NC_022663.1:s-e(+)  ->  ><tag>_REF_kansasii_NC022663
  sed -E "1s/^>REF_[^|]*\|[^|]*\|([A-Z]+)_?([0-9]+)\.[0-9]+:.*/>${tag}_REF_kansasii_\1\2/" "$refgene" > "$in"
  # other complex-species references directly under the kansasii REF, in species_refs.tsv order:
  # ><tag>_REF_<species>_<GCF accession without version and underscore>
  while IFS=$'\t' read -r sp acc _path; do
    [[ -z "$sp" || "$sp" == \#* ]] && continue
    r="$ARM/work/species_refs/REFSP_$sp/genes/$tag.fa"
    acc_id=${acc%.*}; acc_id=${acc_id//_/}   # GCF_002086675.1 -> GCF002086675
    [[ -e "$r" ]] && sed "1s/.*/>${tag}_REF_${sp}_${acc_id}/" "$r" >> "$in"
  done < "$CONFIG/species_refs.tsv"
  # control strains (LABEL set, e.g. Mkan329-126..133) first, then samples by species then id;
  # rows named <gene>_<id>_<species|abbreviated label>
  python3 "$HERE/sort_samples.py" "$SCREENING_MAP_RESULTS" "$ARM/work/samples" "$tag" >> "$in"
  found=$(grep "^>" "$in" | grep -vc "^>${tag}_REF_" || true)
  if (( found == 0 )); then
    echo "[skip] $tag: reference only, no sample sequences"; continue
  fi

  # --adjustdirectionaccurately is intentionally NOT used: 03 already
  # strand-corrects, and direction flipping would silently rename sequences.
  mafft --auto --thread "$THREADS" --preservecase "$in" > "$ALN/$tag.nt.aln.fasta" 2> "$ALN/$tag.mafft.log"

  if is_coding "$tag"; then
    if command -v macse >/dev/null 2>&1; then
      macse -prog alignSequences -seq "$in" \
        -out_NT "$ALN/$tag.codon.aln.fasta" -out_AA "$ALN/$tag.aa.aln.fasta" \
        >> "$ALN/$tag.mafft.log" 2>&1 || echo "[warn] macse failed for $tag, using nt alignment"
    else
      # translatorx-style fallback with MAFFT on the protein level
      python3 "$HERE/backtranslate.py" --nt "$in" --threads "$THREADS" \
        --out-nt "$ALN/$tag.codon.aln.fasta" --out-aa "$ALN/$tag.aa.aln.fasta"
    fi
  fi
  echo "[ok] $tag  ($((found)) samples + reference)"
done

echo
echo "Alignments in $ALN/"
echo "  *.nt.aln.fasta     plain nucleotide - open these in Jalview"
echo "  *.codon.aln.fasta  codon-aware nucleotide (coding genes) - use for codon calls"
echo "  *.aa.aln.fasta     protein alignment"
