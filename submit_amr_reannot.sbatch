#!/usr/bin/env bash
# After the MIC analysis and the results table: 05 again (amr, amr_bam), then 07 (asm vs bam) and 08 (MUT calls with MIC).
#SBATCH --job-name=amr-reannot
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --output=/shares/sander.imm.uzh/MM/kansasii/output/kansasii-amr-marker-screen/logs/%x-%j.log

set -euo pipefail
source "${KANSASII_ROOT:-/shares/sander.imm.uzh/MM/kansasii}/repos/kansasii-docs/templates/kansasii_paths.sh"
REPO=$REPOS/kansasii-amr-marker-screen
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate env_amr
AMR_DIR=$OUT_AMR/results/assembly python3 "$REPO/scripts/05_annotate.py"
AMR_DIR=$OUT_AMR/results/bam python3 "$REPO/scripts/05_annotate.py"
python3 "$REPO/scripts/07_compare_modes.py"
python3 "$REPO/scripts/08_mut_with_mic.py"

# provenance (see kansasii-docs/docs/PROVENANCE.md)
bash "$MANIFEST_SCRIPT" "$OUT_AMR/results/assembly" "$REPO" submit_amr_reannot.sbatch "$SCREENING_MAP_RESULTS" "$OUT_MIC/results"
