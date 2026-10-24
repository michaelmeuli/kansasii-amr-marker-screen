#!/usr/bin/env bash
#SBATCH --job-name=arm-bam
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --array=0-185%20   # one task per sample; pass --array=0-<N-1>%20 when the sample count changes
#SBATCH --output=/shares/sander.imm.uzh/MM/kansasii/output/arm/logs/%x-%A_%a.log
# Prepare once (./run_arm_bam.sh prepare), then sbatch this array; afterwards run
#   sbatch --dependency=afterok:<arrayjob> --wrap './run_arm_bam.sh finish'  (inside env_amr)
set -euo pipefail
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate env_amr
H=/shares/sander.imm.uzh/MM/kansasii/repos/kansasii-amr-marker-screen
ROOT=${KANSASII_ROOT:-/shares/sander.imm.uzh/MM/kansasii}
export ARM_DIR=$ROOT/output/amr_bam
mapfile -t IDS < <(ls "$ROOT/runs/mkan329/assembly/results")
"$H/06_bam_sample.sh" "${IDS[$SLURM_ARRAY_TASK_ID]}"
