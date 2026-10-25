#!/usr/bin/env bash
# Read-mapping variant, one array task per sample. Prepare once (scripts/run_arm_bam.sh prepare), then
#   sbatch --array=0-<N-1>%20 submit_arm_bam.sbatch
# and afterwards scripts/run_arm_bam.sh finish (inside env_amr); scripts/../../immensekansasii/scripts/run_downstream.sh --with-bam does all three.
#SBATCH --job-name=arm-bam
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --array=0-185%20   # default: 186 samples; pass --array=0-<N-1>%20 when the sample count changes
#SBATCH --output=/shares/sander.imm.uzh/MM/kansasii/output/arm/logs/%x-%A_%a.log

set -euo pipefail
source "${KANSASII_ROOT:-/shares/sander.imm.uzh/MM/kansasii}/repos/kansasii-docs/templates/kansasii_paths.sh"
REPO=$REPOS/kansasii-amr-marker-screen
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate env_amr
export ARM_DIR=$OUTPUT_ROOT/amr_bam
mapfile -t IDS < <(ls "$ASSEMBLIES")
"$REPO/scripts/06_bam_sample.sh" "${IDS[$SLURM_ARRAY_TASK_ID]}"
