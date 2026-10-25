#!/usr/bin/env bash
# AMR-marker screen on the Mkan329 assemblies (asm mode): 01 -> 02 -> 03 per new sample -> 03b -> 04 -> 05.
# Output: output/arm/{refs,work}. Run submit_arm_reannot.sbatch afterwards (needs MIC + results table).
#SBATCH --job-name=kansasii-amr-marker-screen
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=/shares/sander.imm.uzh/MM/kansasii/output/arm/logs/%x-%j.log

set -euo pipefail
source "${KANSASII_ROOT:-/shares/sander.imm.uzh/MM/kansasii}/repos/kansasii-docs/templates/kansasii_paths.sh"
REPO=$REPOS/kansasii-amr-marker-screen
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate env_amr
"$REPO/scripts/run_arm.sh"

# provenance (see kansasii-docs/docs/PROVENANCE.md)
bash "$MANIFEST_SCRIPT" "$OUTPUT_ROOT/arm" "$REPO" submit_arm.sbatch "$ASSEMBLIES" "$GTDB_REPRESENTATIVES" "$SCREENING_MAP_RESULTS"
