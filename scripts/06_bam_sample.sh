#!/usr/bin/env bash
# 06_bam_sample.sh - read-mapping genotyping of one sample: map reads to refs/reference.fna,
# extract the target genes (03_sample_seqs.sh bam), keep only a target-region BAM.
#
# Usage: AMR_DIR=<bam-mode dir> ./06_bam_sample.sh SAMPLE_ID [THREADS]
# Inputs : ASSEMBLIES/<id>/0_trimming/<id>_r{1,2}.fastq.gz
# Outputs: $AMR/work/samples/<id>/ (as 03), $AMR/bam/<id>.targets.bam (+ .bai)
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/env.sh"

ID=${1:?sample id}
THREADS=${2:-${SLURM_CPUS_PER_TASK:-4}}
command -v bwa >/dev/null 2>&1 || export PATH=$PATH:${CONDA_ENVS_DIR:-$HOME/data/conda/envs}/snippy/bin   # bwa lives in the snippy env
R=$ASSEMBLIES/$ID/0_trimming
REF=$AMR/refs/reference.fna
mkdir -p "$AMR/bam"

[[ -s "$REF.bwt" ]] || { echo "index missing: run run_amr_bam.sh first" >&2; exit 1; }
TMP=$(mktemp -d "${TMPDIR:-/tmp}/amr_bam.XXXXXX"); trap 'rm -rf "$TMP"' EXIT

bwa mem -t "$THREADS" "$REF" "$R/${ID}_r1.fastq.gz" "$R/${ID}_r2.fastq.gz" 2>/dev/null \
  | samtools sort -@ 2 -o "$TMP/$ID.bam" -
samtools index "$TMP/$ID.bam"

"$HERE/03_sample_seqs.sh" bam "$ID" "$TMP/$ID.bam"

samtools view -b -L "$AMR/work/targets.padded.bed" -o "$AMR/bam/$ID.targets.bam" "$TMP/$ID.bam"
samtools index "$AMR/bam/$ID.targets.bam"
