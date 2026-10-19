#!/usr/bin/env bash
#SBATCH --job-name=arm-reannot
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --output=%x-%j.log

set -euo pipefail
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate env_arm
HERE=/shares/sander.imm.uzh/MM/kansasii/repos/kansasii-arm
ROOT=/shares/sander.imm.uzh/MM/kansasii
ARM_DIR=$ROOT/output/arm python3 "$HERE/05_annotate.py"
ARM_DIR=$ROOT/output/arm_bam python3 "$HERE/05_annotate.py"
python3 "$HERE/07_compare_modes.py"
python3 "$HERE/08_mut_with_mic.py"

# provenance (see kansasii-docs/docs/PROVENANCE.md)
bash /shares/sander.imm.uzh/MM/kansasii/repos/kansasii-docs/bin/write_manifest.sh \
  "/shares/sander.imm.uzh/MM/kansasii/output/arm" \
  /shares/sander.imm.uzh/MM/kansasii/repos/kansasii-arm submit_arm_reannot.sbatch /shares/sander.imm.uzh/MM/kansasii/data/imm/screening_map_results.csv /shares/sander.imm.uzh/MM/kansasii/output/mic
