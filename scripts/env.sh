# Sourced by the pipeline shell scripts: shared project paths (kansasii_paths.sh) plus
#   ARM_DIR   where refs/ and work/ live (default: $OUTPUT_ROOT/arm; run_arm_bam.sh uses $OUTPUT_ROOT/amr_bam)
#   HERE      this scripts/ directory;  CONFIG  ../config (markers.tsv, targets.tsv, species_refs.tsv)
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "${KANSASII_ROOT:-/shares/sander.imm.uzh/MM/kansasii}/repos/kansasii-docs/templates/kansasii_paths.sh"
CONFIG=$(cd "$HERE/../config" && pwd)
ARM=${ARM_DIR:-$OUTPUT_ROOT/arm}
