#!/usr/bin/env bash
# run_arm_bam.sh - read-mapping variant of run_arm.sh (all Mkan329 samples).
#
# Usage: ./run_arm_bam.sh prepare     # once: shared refs, bwa index, targets
#        ./run_arm_bam.sh finish      # after all samples: 04_align + 05_annotate
# Per-sample work is 06_bam_sample.sh (see submit_arm_bam.sbatch for the Slurm array).
# Output dir: $ARM_DIR (default $KANSASII_ROOT/output/amr_bam); refs/ is shared with output/arm.
set -euo pipefail
export ARM_DIR=${ARM_DIR:-${KANSASII_ROOT:-/shares/sander.imm.uzh/MM/kansasii}/output/amr_bam}
source "$(dirname "${BASH_SOURCE[0]}")/env.sh"
command -v bwa >/dev/null 2>&1 || export PATH=$PATH:${CONDA_ENVS_DIR:-$HOME/data/conda/envs}/snippy/bin   # bwa lives in the snippy env

case "${1:?prepare|finish}" in
prepare)
  mkdir -p "$ARM"
  [[ -e "$ARM/refs" ]] || ln -s "$ROOT/output/arm/refs" "$ARM/refs"
  python3 "$HERE/02_make_targets.py"
  [[ -s "$ARM/refs/reference.fna.bwt" ]] || bwa index "$ARM/refs/reference.fna"
  ;;
finish)
  "$HERE/03b_species_refs.sh"
  "$HERE/04_align.sh" "${SLURM_CPUS_PER_TASK:-4}"
  python3 "$HERE/05_annotate.py"
  bash "$ROOT/repos/kansasii-docs/bin/write_manifest.sh" "$ARM" "$HERE" run_arm_bam.sh   # provenance
  ;;
*) echo "usage: $0 prepare|finish" >&2; exit 1 ;;
esac
