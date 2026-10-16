#!/usr/bin/env python3
"""
05_annotate.py - map published resistance markers onto your MAFFT alignment and
write a Jalview features file plus a machine-readable call table.

The hard part this solves: almost no marker in markers.tsv is expressed in
M. kansasii coordinates. rrl/rrs positions are E. coli numbering; rpoB, gyrA,
rpsL, gid codons are M. tuberculosis or M. avium numbering. Writing those
numbers straight onto an M. kansasii alignment is wrong, sometimes by tens of
residues. So for each marker this script:

  1. aligns the M. kansasii reference gene to the external numbering reference
     (protein-level for CDS, nucleotide-level for rRNA),
  2. walks the pairwise alignment to convert external position -> M. kansasii
     position,
  3. converts M. kansasii position -> column in your multi-sample alignment,
  4. reads off what each sample actually has at that column.

Markers with numbering_ref = NONE are reported but not placed; they need manual
curation against the original paper's amplicon.

Sample/gene QC: a call is only trusted when the sample is a non-contaminated member of the
M. kansasii complex (species + contamination_flag from screening_map_results.csv) AND the
gene was extracted with >= --min-ident % BLAST identity and >= --min-cov coverage of the
reference gene. Otherwise the `qc` column says why, `differs` is NOCALL and aa_call is empty.
marker_calls.pass.tsv holds only the rows with qc == ok.

<gene>.aa.calls.jalview_features.txt lists only mutated positions: features typed
"<call>:<aa>" (e.g. MUT:R) for call MUT (red), OTHER or DIFF (orange). MUT = sample has the
marker's mutant residue; OTHER = a residue that is neither the marker's wt nor mut; DIFF =
marker without wt/mut letters whose amino acid differs from the reference. wt, SYN, NOCALL and
UNVERIFIED get no feature. The reference row is flagged MUT when it carries the mutant allele
(ref_aa_check == mut) and OTHER when it carries neither; SAMPLE MUT at such a marker is then
non-discriminating (see below).

`discriminating` is "no" when the reference strain carries the marker's mutant allele
(ref_aa_check == mut; e.g. rpsL K43R, rrs A514C: ATCC 12478 has the "mutant" residue), so an
aa_call of MUT only means "same as the reference" and is not evidence of resistance; "yes" when the
reference carries the wild-type allele; empty when the reference check is unavailable/suspect.

Usage:
  python3 05_annotate.py            # defaults: $KANSASII_ROOT/output/arm/{refs,work}
  python3 05_annotate.py --aln-dir DIR --refs DIR --markers markers.tsv \
      --resolved targets.resolved.tsv --outdir DIR
"""
from __future__ import annotations

import argparse
import csv
import os
import re
import subprocess
import sys
import tempfile
from collections import defaultdict
from pathlib import Path
from typing import Any

ROOT = Path(os.environ.get("KANSASII_ROOT", "/shares/sander.imm.uzh/MM/kansasii"))
ARM = Path(os.environ.get("ARM_DIR", ROOT / "output" / "arm"))
HERE = Path(__file__).resolve().parent

SPECIES_REF_PREFIX = "REF_"  # other complex-species references (<gene>_REF_<species>_<GCF>); never called
COMPLEX = {"kansasii", "persicum", "pseudokansasii", "innocens", "attenuatum", "ostraviense", "gastri"}

CODONS: dict[str, str] = {}
_B, _AA = "TCAG", "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
for i, b1 in enumerate(_B):
    for j, b2 in enumerate(_B):
        for k, b3 in enumerate(_B):
            CODONS[b1 + b2 + b3] = _AA[i * 16 + j * 4 + k]

# Jalview marker features are coloured by the nucleotide the sequence carries
BASE_COLOURS = {
    "A": "64f73f",  # green
    "C": "ffb340",  # yellow-orange
    "G": "e41a1c",  # red
    "T": "3c88ee",  # blue
    "N": "999999",  # N / other IUPAC
}

# protein-coordinate Jalview features: by amino acid (Zappo-like physico-chemical classes)
_AA_CLASSES = {
    "ffafaf": "ILVAM",  # aliphatic / hydrophobic
    "ffc800": "FWY",    # aromatic
    "6464ff": "KRH",    # positive
    "ff0000": "DE",     # negative
    "00ff00": "STNQ",   # hydrophilic
    "ff00ff": "PG",     # conformationally special
    "ffff00": "C",      # cysteine
}
AA_COLOURS = {aa: col for col, aas in _AA_CLASSES.items() for aa in aas}
AA_COLOURS["STOP"] = "000000"
AA_COLOURS["X"] = "999999"

# ... or by call state, for mutated positions only: MUT red, OTHER / DIFF orange.
# Feature type is "<call>:<aa>" so the residue is visible; wt / SYN / NOCALL / UNVERIFIED
# positions get no feature at all.
STATE_COLOURS = {
    "MUT": "e41a1c",
    "OTHER": "ff7f00",
    "DIFF": "ff7f00",
}


def read_fasta(path: str) -> list[tuple[str, str]]:
    out: list[tuple[str, str]] = []
    name: str | None = None
    buf: list[str] = []
    for line in open(path, encoding="utf-8"):
        if line.startswith(">"):
            if name:
                out.append((name, "".join(buf)))
            name, buf = line[1:].strip(), []
        else:
            buf.append(line.strip())
    if name:
        out.append((name, "".join(buf)))
    return out


def translate(nt: str) -> str:
    nt = nt.upper().replace("-", "")
    return "".join(CODONS.get(nt[i:i + 3], "X") for i in range(0, len(nt) - len(nt) % 3, 3))


def mafft_pair(a_name: str, a_seq: str, b_name: str, b_seq: str) -> tuple[str, str]:
    """Align two sequences, return the aligned pair."""
    with tempfile.NamedTemporaryFile("w", suffix=".fa", delete=False, encoding="utf-8") as tf:
        tf.write(f">{a_name}\n{a_seq}\n>{b_name}\n{b_seq}\n")
        p = tf.name
    res = subprocess.run(["mafft", "--auto", "--quiet", "--preservecase", p],
                         capture_output=True, text=True, check=True)
    with tempfile.NamedTemporaryFile("w", suffix=".fa", delete=False, encoding="utf-8") as tf:
        tf.write(res.stdout)
        q = tf.name
    recs = dict(read_fasta(q))
    os.unlink(p)
    os.unlink(q)
    return recs[a_name], recs[b_name]


def transfer(ext_seq: str, mkn_seq: str, ext_pos: int) -> tuple[int | None, tuple[str | None, str | None]]:
    """External 1-based position -> M. kansasii 1-based position, or None."""
    ea, ma = mafft_pair("EXT", ext_seq, "MKN", mkn_seq)
    ei = mi = 0
    for ec, mc in zip(ea, ma):
        if ec != "-":
            ei += 1
        if mc != "-":
            mi += 1
        if ei == ext_pos:
            return (mi if mc != "-" else None), (ec, mc)
    return None, (None, None)


def ungapped_to_column(aligned_seq: str) -> dict[int, int]:
    """Map 1-based ungapped position -> 0-based alignment column."""
    m: dict[int, int] = {}
    u = 0
    for col, ch in enumerate(aligned_seq):
        if ch != "-":
            u += 1
            m[u] = col
    return m


def row_suffix(sname: str, tag: str) -> str:
    """Alignment row name `<gene>_<rest>` -> `<rest>` (gene tags such as folP_c1 contain `_`)."""
    return sname[len(tag) + 1:] if sname.startswith(tag + "_") else sname


def sample_id(sname: str, tag: str) -> str:
    """Plain sample id from an alignment row name `<gene>_<id>_<species>` (ids contain no `_`)."""
    return row_suffix(sname, tag).split("_", 1)[0]


def load_meta(path: str) -> dict[str, tuple[str, str]]:
    """sample id -> (species, contamination_flag) from screening_map_results.csv."""
    out: dict[str, tuple[str, str]] = {}
    if not os.path.exists(path):
        print(f"[warn] no sample metadata at {path}: species/contamination QC skipped", file=sys.stderr)
        return out
    with open(path, encoding="utf-8", newline="") as fh:
        for row in csv.DictReader(fh):
            out[row["PROBENNUMMER"]] = (row.get("species", "").strip(), row.get("contamination_flag", "").strip())
    return out


def load_mixed(samples_dir: str, sample: str) -> dict[int, tuple[str, str, float]]:
    """Genomic position -> (ref, alt, alt fraction) of the 10-90 % alt-allele sites that
    03_sample_seqs.sh (bam mode) wrote to mixed_sites.tsv; empty in asm mode."""
    path = os.path.join(samples_dir, sample, "mixed_sites.tsv")
    out: dict[int, tuple[str, str, float]] = {}
    if not os.path.exists(path):
        return out
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 5:
                continue
            ref_n, alt_n = (int(x) for x in f[4].split(",")[:2])
            if ref_n + alt_n:
                out[int(f[1])] = (f[2], f[3], alt_n / (ref_n + alt_n))
    return out


def load_hit(samples_dir: str, sample: str, tag: str) -> tuple[float, float] | None:
    """(percent identity, fraction of the reference gene covered) of the extracted gene, or None."""
    path = os.path.join(samples_dir, sample, "blast", f"{tag}.hit")
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return None
    with open(path, encoding="utf-8") as fh:
        f = fh.readline().split()
    return float(f[3]), float(f[4]) / float(f[5])


def sample_qc(sample: str, tag: str, meta: dict[str, tuple[str, str]], samples_dir: str,
              min_ident: float, min_cov: float) -> str:
    """'ok' or a ';'-joined list of reasons this sample/gene must not be called."""
    why: list[str] = []
    species, contam = meta.get(sample, ("", ""))
    if contam:
        why.append("contaminated")
    if species and species not in COMPLEX:
        why.append(f"non_complex({species})")
    hit = load_hit(samples_dir, sample, tag)
    if hit is None:
        if os.path.isdir(os.path.join(samples_dir, sample, "blast")):   # BLAST QC exists only in asm mode
            why.append("no_hit")
    else:
        pid, cov = hit
        if pid < min_ident:
            why.append(f"low_identity({pid:.1f})")
        if cov < min_cov:
            why.append(f"partial({cov:.0%})")
    return ";".join(why) or "ok"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--aln-dir", default=str(ARM / "work" / "aln"))
    ap.add_argument("--refs", default=str(ARM / "refs"))
    ap.add_argument("--markers", default=str(HERE / "markers.tsv"))
    ap.add_argument("--resolved", default=str(ARM / "work" / "targets.resolved.tsv"))
    ap.add_argument("--outdir", default=str(ARM / "work" / "report"))
    ap.add_argument("--samples", default=str(ARM / "work" / "samples"))
    ap.add_argument("--sample-meta", default=str(ROOT / "output" / "screening_map_results.csv"))
    ap.add_argument("--min-ident", type=float, default=80.0,
                    help="minimum BLAST %% identity of the extracted gene to the reference gene "
                         "(a floor against wrong-gene hits; species gating does the real work, and "
                         "complex species differ a lot at e.g. gid: 88%% persicum, 80%% attenuatum)")
    ap.add_argument("--min-cov", type=float, default=0.9,
                    help="minimum fraction of the reference gene covered by the extracted gene")
    ap.add_argument("--min-cov-gene", action="append", default=[], metavar="GENE=FRAC",
                    help="per-gene override of --min-cov (repeatable); built-in: gyrA=0.6, because "
                         "the reference gyrA carries a ~1.26 kb insertion that most samples lack")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    meta = load_meta(args.sample_meta)
    cov_gene = {"gyrA": 0.6}
    for item in args.min_cov_gene:
        g, _, v = item.partition("=")
        cov_gene[g] = float(v)
    qc_cache: dict[tuple[str, str], str] = {}
    numbering: dict[str, str] = {}
    coords: dict[str, tuple[int, int, str]] = {}   # gene tag -> (start, end, strand)
    mixed_cache: dict[str, dict[int, tuple[str, str, float]]] = {}
    with open(args.resolved, encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            numbering[r["gene"]] = r["numbering_ref"]
            coords[r["gene"]] = (int(r["start_1based"]), int(r["end_1based"]), r["strand"])

    markers: defaultdict[str, list[dict[str, str]]] = defaultdict(list)
    with open(args.markers, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            markers[f[1]].append(dict(zip(
                ["drug", "gene", "kind", "marker", "position", "numbering_ref",
                 "evidence_species", "evidence", "notes", "doi"], f)))

    calls_path = os.path.join(args.outdir, "marker_calls.tsv")
    unplaced: list[dict[str, str]] = []
    pass_path = os.path.join(args.outdir, "marker_calls.pass.tsv")
    with open(calls_path, "w", newline="", encoding="utf-8") as cf, \
         open(pass_path, "w", newline="", encoding="utf-8") as pf:
        cw = csv.writer(cf, delimiter="\t", lineterminator="\n")
        pw = csv.writer(pf, delimiter="\t", lineterminator="\n")
        header = ["gene", "drug", "marker", "evidence", "numbering_ref",
                  "ext_position", "mkn_position", "aln_column",
                  "ref_state", "sample", "sample_state", "differs",
                  "ref_aa", "sample_aa", "marker_wt_aa", "marker_mut_aa",
                  "ref_aa_check", "aa_call", "discriminating", "species", "contamination_flag", "hit_pident",
                  "hit_cov", "qc", "mixed", "doi"]
        cw.writerow(header)
        pw.writerow(header)

        for aln_file in sorted(os.listdir(args.aln_dir)):
            if not aln_file.endswith((".codon.aln.fasta", ".nt.aln.fasta")):
                continue
            tag = aln_file.split(".")[0]
            gene = re.sub(r"_c\d+$", "", tag)
            if gene not in markers:
                continue
            # prefer the codon-aware alignment when both exist
            codon = os.path.join(args.aln_dir, f"{tag}.codon.aln.fasta")
            path = codon if os.path.exists(codon) else os.path.join(args.aln_dir, aln_file)
            if path != os.path.join(args.aln_dir, aln_file) and aln_file.endswith(".nt.aln.fasta"):
                continue

            recs = read_fasta(path)
            ref_name, ref_aln = recs[0]
            ref_ungapped = ref_aln.replace("-", "")
            colmap = ungapped_to_column(ref_aln)
            is_cds = os.path.exists(codon)

            feats: list[tuple[dict[str, str], int, int, str, int]] = []
            feat_ref_check: dict[int, str] = {}  # feature index -> ref_check (wt/mut/neither/"")
            feat_wt_mut: dict[int, tuple[str, str]] = {}  # feature index -> (wt_aa, mut_aa)
            calls: dict[tuple[int, str], str] = {}  # (feature index, sequence) -> state
            for m in markers[gene]:
                nref, pos = m["numbering_ref"], m["position"]
                if nref == "NONE" or pos == "NA":
                    unplaced.append(m)
                    continue

                ext_path = os.path.join(args.refs, f"{nref}.fa")
                if not os.path.exists(ext_path):
                    unplaced.append({**m, "notes": m["notes"] + f" [missing {ext_path}]"})
                    continue
                ext_seq = read_fasta(ext_path)[0][1]
                mkn_cmp = translate(ref_ungapped) if nref.endswith("_PROT") else ref_ungapped

                lo, hi = (pos.split("-") + [pos])[:2] if "-" in pos else (pos, pos)
                mkn_lo: int | None = None
                mkn_hi: int | None = None
                for p in (int(lo), int(hi)):
                    mkn_p, (ec, mc) = transfer(ext_seq, mkn_cmp, p)
                    if mkn_p is None:
                        continue
                    if mkn_lo is None:
                        mkn_lo, ref_state_ext = mkn_p, (ec, mc)
                    mkn_hi = mkn_p
                if mkn_lo is None or mkn_hi is None:
                    unplaced.append({**m, "notes": m["notes"] + " [no alignment anchor]"})
                    continue

                # protein position -> first nt of that codon
                def to_nt(p: int) -> int:
                    return (p - 1) * 3 + 1 if nref.endswith("_PROT") else p

                nt_lo, nt_hi = to_nt(mkn_lo), to_nt(mkn_hi) + (2 if nref.endswith("_PROT") else 0)
                c_lo, c_hi = colmap.get(nt_lo), colmap.get(nt_hi)
                if c_lo is None:
                    unplaced.append({**m, "notes": m["notes"] + " [outside alignment]"})
                    continue
                c_hi = c_hi if c_hi is not None else c_lo

                ref_state = ref_aln[c_lo:c_hi + 1]
                feats.append((m, c_lo, c_hi, ref_state, mkn_lo))

                # Coding markers are compared at the amino-acid level so synonymous
                # changes (SYN) are not reported as differences.
                use_aa = nref.endswith("_PROT") and "-" not in ref_state and len(ref_state) % 3 == 0
                ref_aa = translate(ref_state) if use_aa else ""
                mm = re.match(r"^([A-Z])\d+([A-Z*])(?![A-Za-z0-9])", m["marker"])
                wt_aa, mut_aa = (mm.group(1), mm.group(2)) if mm and use_aa else ("", "")
                # single-nucleotide markers (rRNA, e.g. A2058G): same wt/MUT/OTHER logic on the base
                nm = re.match(r"^([ACGTU])\d+([ACGTU])(?![A-Za-z0-9])", m["marker"])
                is_nt_marker = bool(nm) and not use_aa and len(ref_state) == 1
                if nm and is_nt_marker:
                    wt_aa, mut_aa = nm.group(1).replace("U", "T"), nm.group(2).replace("U", "T")
                # does the reference itself carry the wild-type or the mutant residue?
                # "neither" means the position transfer is probably off.
                ref_check = ""
                if wt_aa and len(ref_aa) == 1:
                    ref_check = "wt" if ref_aa == wt_aa else ("mut" if ref_aa == mut_aa else "neither")
                elif is_nt_marker:
                    rb = ref_state.upper()
                    ref_check = "wt" if rb == wt_aa else ("mut" if rb == mut_aa else "neither")
                feat_ref_check[len(feats) - 1] = ref_check
                feat_wt_mut[len(feats) - 1] = (wt_aa, mut_aa)

                # genomic window of this marker in the reference gene, for the mixed-site lookup
                g_start, g_end, g_strand = coords.get(tag, (0, 0, "+"))
                g_a, g_b = ((g_start + nt_lo - 1, g_start + nt_hi - 1) if g_strand == "+"
                            else (g_end - nt_hi + 1, g_end - nt_lo + 1))

                for sname, saln in recs[1:]:
                    if row_suffix(sname, tag).startswith(SPECIES_REF_PREFIX):
                        continue
                    s_state = saln[c_lo:c_hi + 1]
                    # missing data (depth-masked N or alignment gap) must never read as
                    # a result: decide NOCALL before comparing to the reference
                    s_up = s_state.upper()
                    s_aa = ""
                    aa_call = ""
                    sid = sample_id(sname, tag)
                    if (sid, tag) not in qc_cache:
                        qc_cache[(sid, tag)] = sample_qc(
                            sid, tag, meta, args.samples, args.min_ident,
                            cov_gene.get(re.sub(r"_c\d+$", "", tag), args.min_cov))
                    qc = qc_cache[(sid, tag)]
                    species, contam = meta.get(sid, ("", ""))
                    hit = load_hit(args.samples, sid, tag)
                    if qc != "ok":
                        # wrong species / contaminated / poorly extracted gene: whatever the
                        # alignment shows at this column is not a trustworthy genotype
                        differs = "NOCALL"
                    elif "N" in s_up or "-" in s_up:
                        differs = "NOCALL"
                    elif use_aa:
                        s_aa = translate(s_up)
                        if s_aa != ref_aa:
                            differs = "YES"
                        else:
                            differs = "SYN" if s_up != ref_state.upper() else "no"
                        if ref_check == "neither":
                            # the reference has neither the marker's wt nor mut residue, so
                            # the numbering does not fit this gene: do not make a call
                            aa_call = "UNVERIFIED"
                        elif wt_aa and wt_aa != mut_aa and len(s_aa) == 1:  # synonymous markers: no aa call
                            aa_call = "MUT" if s_aa == mut_aa else ("wt" if s_aa == wt_aa else "OTHER")
                    else:
                        differs = "YES" if s_up != ref_state.upper() else "no"
                        if is_nt_marker:
                            if ref_check == "neither":
                                aa_call = "UNVERIFIED"
                            elif len(s_up) == 1:
                                aa_call = "MUT" if s_up == mut_aa else ("wt" if s_up == wt_aa else "OTHER")
                    if sid not in mixed_cache:
                        mixed_cache[sid] = load_mixed(args.samples, sid)
                    mixed = ";".join(f"{pos}:{r}>{a}@{fr:.2f}" for pos, (r, a, fr)
                                     in sorted(mixed_cache[sid].items()) if g_a <= pos <= g_b)
                    # the reference strain carries the marker's mutant allele: MUT is then the
                    # reference/species baseline, so the marker cannot be read as resistance
                    discriminating = {"wt": "yes", "mut": "no"}.get(ref_check, "")
                    out_row = [tag, m["drug"], m["marker"], m["evidence"], nref,
                               pos, mkn_lo, c_lo + 1, ref_state, f"{tag}__{sid}", s_state,
                               differs, ref_aa, s_aa, wt_aa, mut_aa, ref_check, aa_call,
                               discriminating,
                               species, contam, f"{hit[0]:.2f}" if hit else "",
                               f"{hit[1]:.3f}" if hit else "", qc, mixed, m["doi"]]
                    cw.writerow(out_row)
                    if qc == "ok":
                        pw.writerow(out_row)
                    if differs == "NOCALL":
                        calls[(len(feats) - 1, sname)] = "NOCALL"
                    elif aa_call:
                        calls[(len(feats) - 1, sname)] = aa_call
                    else:
                        calls[(len(feats) - 1, sname)] = {"YES": "DIFF", "SYN": "SYN"}.get(differs, "wt")

            # ---- Jalview features file -------------------------------------
            jf = os.path.join(args.outdir, f"{tag}.jalview_features.txt")
            # Jalview feature coords are 1-based on the named sequence, in UNGAPPED
            # residue numbering, so build the column -> residue index map per sequence.
            col_idx: dict[str, dict[int, int]] = {}
            for sname, saln in recs:
                smap: dict[int, int] = {}
                u = 0
                for col, ch in enumerate(saln):
                    if ch != "-":
                        u += 1
                        smap[col] = u
                col_idx[sname] = smap
            lines: list[str] = []
            used: set[str] = set()
            for m, c_lo, c_hi, ref_state, mkn_p in feats:
                desc = f"{m['drug']} | {m['marker']} | {m['evidence']} | {m['evidence_species']} | doi:{m['doi']}"
                # one single-residue feature per marker nucleotide, typed by the base that
                # sequence carries there, so the colour shows the allele
                for sname, saln in recs:
                    for col in range(c_lo, c_hi + 1):
                        idx = col_idx[sname].get(col)
                        if idx is None:
                            continue
                        base = saln[col].upper().replace("U", "T")
                        btype = base if base in "ACGT" else "N"
                        used.add(btype)
                        lines.append(f"{desc}\t{sname}\t-1\t{idx}\t{idx}\t{btype}\n")
            with open(jf, "w", encoding="utf-8") as fh:
                for btype in sorted(used):
                    fh.write(f"{btype}\t{BASE_COLOURS[btype]}\n")
                fh.write("\nSTARTGROUP\tAMR\n")
                fh.writelines(lines)
                fh.write("ENDGROUP\tAMR\n")
            print(f"[ok] {tag}: {len(feats)} markers placed -> {jf}")

            # ---- protein-coordinate features (for <gene>.aa.aln.fasta) -------
            aa_path = os.path.join(args.aln_dir, f"{tag}.aa.aln.fasta")
            if is_cds and os.path.exists(aa_path):
                aa_recs = read_fasta(aa_path)
                aa_by_name = dict(aa_recs)
                aa_ref = aa_recs[0][1]
                aa_idx: dict[str, dict[int, int]] = {}
                for sname, saa in aa_recs:
                    amap: dict[int, int] = {}
                    u = 0
                    for col, ch in enumerate(saa):
                        if ch != "-":
                            u += 1
                            amap[col] = u
                    aa_idx[sname] = amap
                by_aa: list[str] = []
                by_state: list[str] = []
                used_aa: set[str] = set()
                used_state: set[str] = set()
                for fi, (m, c_lo, c_hi, _rs, _mp) in enumerate(feats):
                    # only codon-aligned windows map to protein columns
                    if c_lo % 3 or (c_hi - c_lo + 1) % 3:
                        continue
                    desc = f"{m['drug']} | {m['marker']} | {m['evidence']} | {m['evidence_species']} | doi:{m['doi']}"
                    for sname, saa in aa_recs:
                        for acol in range(c_lo // 3, c_hi // 3 + 1):
                            idx = aa_idx[sname].get(acol)
                            if idx is None or acol >= len(saa):
                                continue
                            aa = saa[acol].upper()
                            atype = "STOP" if aa == "*" else (aa if aa in AA_COLOURS else "X")
                            used_aa.add(atype)
                            by_aa.append(f"{desc}\t{sname}\t-1\t{idx}\t{idx}\t{atype}\n")
                            if sname == ref_name:
                                stype = "REF"  # replaced below by the reference's own call
                            elif (fi, sname) in calls:
                                stype = calls[(fi, sname)]
                            elif "X" in aa or aa == "-":  # species refs (no per-sample call)
                                stype = "NOCALL"
                            else:
                                # species references: no per-sample call. Judge against the
                                # marker's own wt/mut residues when it has them (the reference
                                # strain may itself carry the mutant), else against the reference.
                                wt_r, mut_r = feat_wt_mut.get(fi, ("", ""))
                                if len(wt_r) == 1 and len(aa) == 1 and wt_r != mut_r:
                                    stype = ("MUT" if aa == mut_r
                                             else "wt" if aa == wt_r else "OTHER")
                                else:
                                    stype = "wt" if aa == aa_ref[acol].upper() else "DIFF"
                            if sname == ref_name:
                                # the reference is only flagged where it carries the mutant
                                # allele (or neither wt nor mut); it is never "DIFF" from itself
                                rc = feat_ref_check.get(fi, "")
                                stype = {"mut": "MUT", "neither": "OTHER"}.get(rc, "")
                            if stype not in STATE_COLOURS:  # wt / SYN / NOCALL / UNVERIFIED: no feature
                                continue
                            sdesc = desc
                            if sname == ref_name and feat_ref_check.get(fi) == "mut":
                                sdesc += " | ref carries mutant (non-discriminating)"
                            ttype = f"{stype}:{atype}"
                            used_state.add(ttype)
                            by_state.append(f"{sdesc}\t{sname}\t-1\t{idx}\t{idx}\t{ttype}\n")
                for suffix, body, used_t, cols in (
                        ("aa.jalview_features.txt", by_aa, used_aa, AA_COLOURS),
                        ("aa.calls.jalview_features.txt", by_state, used_state, STATE_COLOURS)):
                    aa_out = os.path.join(args.outdir, f"{tag}.{suffix}")
                    if not body:
                        if os.path.exists(aa_out):
                            os.remove(aa_out)  # stale file from an earlier run
                        continue
                    with open(aa_out, "w", encoding="utf-8") as fh:
                        for t in sorted(used_t):
                            fh.write(f"{t}\t{cols[t.split(':')[0]]}\n")
                        fh.write("\nSTARTGROUP\tAMR\n")
                        fh.writelines(body)
                        fh.write("ENDGROUP\tAMR\n")
                    print(f"[ok] {tag}: protein features -> {aa_out}")

    up = os.path.join(args.outdir, "unplaced_markers.tsv")
    with open(up, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["drug", "gene", "marker", "evidence", "reason_or_notes", "doi"])
        for m in unplaced:
            w.writerow([m["drug"], m["gene"], m["marker"], m["evidence"], m["notes"], m["doi"]])

    print(f"\ncalls:    {calls_path}  (qc == ok only: {pass_path})")
    print(f"unplaced: {up}  <- review these by hand, they are not failures of the data")


if __name__ == "__main__":
    main()
