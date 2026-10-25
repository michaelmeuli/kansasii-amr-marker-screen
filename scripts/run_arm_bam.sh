#!/usr/bin/env bash
# run_arm_bam.sh - read-mapping variant of run_arm.sh (all Mkan329 samples).
#
# Usage: ./run_arm_bam.sh prepare     # once: shared refs, bwa index, targets
#        ./run_arm_bam.sh finish      # after all samples: 04_align + 05_annotate
# Per-sample work is 06_bam_sample.sh (see submit_arm_bam.sbatch for the Slurm array).
# Output dir: $ARM_DIR (default OUTPUT_ROOT/amr_bam); refs/ is shared with output/arm.
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/env.sh"
export ARM_DIR=${ARM_DIR:-$OUTPUT_ROOT/amr_bam}
ARM=$ARM_DIR
command -v bwa >/dev/null 2>&1 || export PATH=$PATH:${CONDA_ENVS_DIR:-$HOME/data/conda/envs}/snippy/bin   # bwa lives in the snippy env

case "${1:?prepare|finish}" in
prepare)
  mkdir -p "$ARM"
  [[ -e "$ARM/refs" ]] || ln -s "$OUTPUT_ROOT/arm/refs" "$ARM/refs"
  python3 "$HERE/02_make_targets.py"
  [[ -s "$ARM/refs/reference.fna.bwt" ]] || bwa index "$ARM/refs/reference.fna"
  ;;
finish)
  "$HERE/03b_species_refs.sh"
  "$HERE/04_align.sh" "${SLURM_CPUS_PER_TASK:-4}"
  python3 "$HERE/05_annotate.py"
  bash "$MANIFEST_SCRIPT" "$ARM" "$REPOS/kansasii-amr-marker-screen" run_arm_bam.sh   # provenance
  ;;
*) echo "usage: $0 prepare|finish" >&2; exit 1 ;;
esac
