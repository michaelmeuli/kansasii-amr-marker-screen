# Sourced by the pipeline shell scripts: shared project paths (kansasii_paths.sh) plus
#   AMR_DIR   where refs/ and work/ live (default: $OUT_AMR/results/assembly; run_amr_bam.sh uses $OUT_AMR/results/bam)
#   HERE      this scripts/ directory;  CONFIG  ../config (markers.tsv, targets.tsv, species_refs.tsv)
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "${KANSASII_ROOT:-/shares/sander.imm.uzh/MM/kansasii}/repos/kansasii-docs/templates/kansasii_paths.sh"
CONFIG=$(cd "$HERE/../config" && pwd)
AMR=${AMR_DIR:-$OUT_AMR/results/assembly}

# targets_done DIR: true if DIR/genes exists and DIR/targets.list (written by 03_sample_seqs.sh)
# covers every current target in $AMR/work/ref_genes, so adding a target re-extracts old samples.
targets_done () {
  [[ -d "$1/genes" && -s "$1/targets.list" ]] || return 1
  local t
  for t in "$AMR"/work/ref_genes/*.fa; do
    grep -qxF "$(basename "$t" .fa)" "$1/targets.list" || return 1
  done
}
