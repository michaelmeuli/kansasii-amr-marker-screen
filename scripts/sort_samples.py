#!/usr/bin/env python3
"""Print the per-sample records of one gene: control strains first, then samples by species and id.

Usage: sort_samples.py screening_map_results.csv samples_dir gene_tag
Sample headers  >gene__ID  become  >gene_ID_species  (species from the `species` column; blank -> unknown).
Control strains (non-empty LABEL, e.g. Mkan329-126..133) come first as  >gene_ID_<abbreviated LABEL>,
e.g. rpsL_Mkan329-126_tuberculosis_ATCC27294.
"""
from __future__ import annotations

import csv
import re
import sys
from pathlib import Path


def abbrev_label(label: str) -> str:
    """'M. tuberculosis H37Rv (ATCC 27294)' -> 'tuberculosis_ATCC27294' (genus and strain name dropped)."""
    m = re.match(r"M\.\s+(\S+)\s*(?:.*?)\(\s*([A-Za-z]+)\s*([\w.-]+)\s*\)", label)
    if m:
        return f"{m.group(1)}_{m.group(2)}{m.group(3)}"
    return re.sub(r"[^\w-]+", "_", label.replace("M. ", "")).strip("_")


def main() -> None:
    meta_path, samples_dir, tag = sys.argv[1:4]
    species: dict[str, str] = {}
    controls: dict[str, str] = {}
    with open(meta_path, encoding="utf-8", newline="") as fh:
        for row in csv.DictReader(fh):
            species[row["PROBENNUMMER"]] = row.get("species", "").strip() or "unknown"
            if row.get("LABEL", "").strip():
                controls[row["PROBENNUMMER"]] = abbrev_label(row["LABEL"].strip())
    items = []
    for fa in Path(samples_dir).glob(f"*/genes/{tag}.fa"):
        sid = fa.parent.parent.name
        if sid in controls:
            suffix, grp = controls[sid], 0
        else:
            suffix, grp = species.get(sid, "unknown").replace(" ", "-").replace("_", "-"), 1
        items.append((grp, suffix if grp else "", sid, suffix, fa))
    for _grp, _k, sid, suffix, fa in sorted(items):
        for line in fa.read_text().splitlines():
            if line.startswith(">"):
                line = f">{tag}_{sid}_{suffix}"
            print(line)


if __name__ == "__main__":
    main()
