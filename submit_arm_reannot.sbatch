#!/usr/bin/env bash
# After the MIC analysis and the results table: 05 again (arm, amr_bam), then 07 (asm vs bam) and 08 (MUT calls with MIC).
#SBATCH --job-name=arm-reannot
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --output=/shares/sander.imm.uzh/MM/kansasii/output/arm/logs/%x-%j.log

set -euo pipefail
source "${KANSASII_ROOT:-/shares/sander.imm.uzh/MM/kansasii}/repos/kansasii-docs/templates/kansasii_paths.sh"
REPO=$REPOS/kansasii-amr-marker-screen
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate env_amr
ARM_DIR=$OUTPUT_ROOT/arm python3 "$REPO/scripts/05_annotate.py"
ARM_DIR=$OUTPUT_ROOT/amr_bam python3 "$REPO/scripts/05_annotate.py"
python3 "$REPO/scripts/07_compare_modes.py"
python3 "$REPO/scripts/08_mut_with_mic.py"

# provenance (see kansasii-docs/docs/PROVENANCE.md)
bash "$MANIFEST_SCRIPT" "$OUTPUT_ROOT/arm" "$REPO" submit_arm_reannot.sbatch "$SCREENING_MAP_RESULTS" "$OUTPUT_ROOT/mic"
