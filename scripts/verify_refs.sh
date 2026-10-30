#!/usr/bin/env bash
# verify_refs.sh - check the numbering references in $AMR/refs against NCBI, and the
# resolved target coordinates against the ATCC 12478 GFF. Prints a PASS/FAIL table.
# Reads config/seq_refs.tsv (file, db, accession, range, ...). Needs network (NCBI efetch).
# Usage: verify_refs.sh [OUTDIR]     (default $AMR/verify_refs; writes sequences.tsv, MANIFEST.txt)
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/env.sh"

OUT=${1:-$AMR/verify_refs}
mkdir -p "$OUT"
seqonly () { grep -v '>' "$1" | tr -d '\n\r ' ; }   # sequence without header/newlines

fail=0
printf "file\taccession\tlocal_len\tncbi_len\tlocal_sha256\tstatus\n" > "$OUT/sequences.tsv"
while IFS=$'\t' read -r file db acc range _; do
  [[ "$file" == \#* ]] && continue
  url="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=$db&id=$acc&rettype=fasta&retmode=text"
  [[ "$range" != "-" ]] && url+="&seq_start=${range%-*}&seq_stop=${range#*-}"
  curl -fsSL "$url" -o "$OUT/$file.ncbi"
  l=$(seqonly "$AMR/refs/$file"); n=$(seqonly "$OUT/$file.ncbi")
  if [[ "$l" == "$n" ]] && head -1 "$AMR/refs/$file" | grep -q "^>$acc"; then st=PASS; else st=FAIL; fail=1; fi
  printf "%s\t%s\t%d\t%d\t%s\t%s\n" "$file" "$acc" "${#l}" "${#n}" \
    "$(printf %s "$l" | sha256sum | cut -c1-16)" "$st" >> "$OUT/sequences.tsv"
  sleep 0.4   # NCBI: <= 3 requests/s without an API key
done < "$CONFIG/seq_refs.tsv"
column -t -s$'\t' "$OUT/sequences.tsv"

# every resolved target must sit on a GFF gene/rRNA feature with the same locus tag and coordinates
python3 -I - "$AMR/work/targets.resolved.tsv" "$AMR/refs/reference.gff" "$OUT/targets.tsv" <<'PY' || fail=1
import csv, re, sys
res, gff, out = sys.argv[1:]
feat = {}
for line in open(gff):
    if line.startswith('#'): continue
    f = line.rstrip('\n').split('\t')
    if len(f) < 9 or f[2] not in ('gene', 'rRNA', 'CDS'): continue
    m = re.search(r'locus_tag=([^;]+)', f[8])
    if m: feat.setdefault(m.group(1), set()).add((f[0], int(f[3]), int(f[4]), f[6]))
bad = 0
with open(out, 'w') as o:
    o.write('gene\tlocus_tag\tcontig\tstart\tend\tstrand\tstatus\n')
    for r in csv.DictReader(open(res), delimiter='\t'):
        key = (r['contig'], int(r['start_1based']), int(r['end_1based']), r['strand'])
        ok = key in feat.get(r['locus_tag'], ())
        bad += not ok
        o.write('\t'.join([r['gene'], r['locus_tag'], *map(str, key), 'PASS' if ok else 'FAIL']) + '\n')
print(f'targets: {bad} of the resolved rows do not match the GFF' if bad else 'targets: all resolved rows match the GFF')
sys.exit(1 if bad else 0)
PY
"${KANSASII_ROOT:-/shares/sander.imm.uzh/MM/kansasii}/repos/kansasii-docs/bin/write_manifest.sh" "$OUT" "$(dirname "$HERE")" scripts/verify_refs.sh \
  "$CONFIG/seq_refs.tsv" "$AMR/refs" "$AMR/work/targets.resolved.tsv"
[[ $fail -eq 0 ]] && echo "ALL PASS" || { echo "SOME CHECKS FAILED (see $OUT)"; exit 1; }
