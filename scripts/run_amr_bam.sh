#!/usr/bin/env bash
# run_amr_bam.sh - read-mapping variant of run_amr.sh (all Mkan329 samples).
#
# Usage: ./run_amr_bam.sh prepare     # once: shared refs, bwa index, targets
#        ./run_amr_bam.sh finish      # after all samples: 04_align + 05_annotate
# Per-sample work is 06_bam_sample.sh (see submit_amr_bam.sbatch for the Slurm array).
# Output dir: $AMR_DIR (default OUTPUT_ROOT/amr_bam); refs/ is shared with output/amr.
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/env.sh"
export AMR_DIR=${AMR_DIR:-$OUTPUT_ROOT/amr_bam}
AMR=$AMR_DIR
command -v bwa >/dev/null 2>&1 || export PATH=$PATH:${CONDA_ENVS_DIR:-$HOME/data/conda/envs}/snippy/bin   # bwa lives in the snippy env

case "${1:?prepare|finish}" in
prepare)
  mkdir -p "$AMR"
  [[ -e "$AMR/refs" ]] || ln -s "$OUTPUT_ROOT/amr/refs" "$AMR/refs"
  python3 "$HERE/02_make_targets.py"
  [[ -s "$AMR/refs/reference.fna.bwt" ]] || bwa index "$AMR/refs/reference.fna"
  ;;
finish)
  "$HERE/03b_species_refs.sh"
  "$HERE/04_align.sh" "${SLURM_CPUS_PER_TASK:-4}"
  python3 "$HERE/05_annotate.py"
  bash "$MANIFEST_SCRIPT" "$AMR" "$REPOS/kansasii-amr-marker-screen" run_amr_bam.sh   # provenance
  ;;
*) echo "usage: $0 prepare|finish" >&2; exit 1 ;;
esac
