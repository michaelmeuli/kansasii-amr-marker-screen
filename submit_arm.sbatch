#!/usr/bin/env bash
#SBATCH --job-name=kansasii-arm
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=/shares/sander.imm.uzh/MM/kansasii/output/arm/logs/%x-%j.log

set -euo pipefail
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate env_arm
/shares/sander.imm.uzh/MM/kansasii/repos/kansasii-arm/run_arm.sh

# provenance (see kansasii-docs/docs/PROVENANCE.md)
bash /shares/sander.imm.uzh/MM/kansasii/repos/kansasii-docs/bin/write_manifest.sh \
  "/shares/sander.imm.uzh/MM/kansasii/output/arm" \
  /shares/sander.imm.uzh/MM/kansasii/repos/kansasii-arm submit_arm.sbatch 
