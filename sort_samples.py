#!/usr/bin/env python3
"""Print the per-sample records of one gene sorted by species then sample id.

Usage: sort_samples.py screening_map_results.csv samples_dir gene_tag
Headers  >gene__ID  become  >gene__ID_species  (species from the `species` column; blank -> unknown).
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path


def main() -> None:
    meta_path, samples_dir, tag = sys.argv[1:4]
    species: dict[str, str] = {}
    with open(meta_path, encoding="utf-8", newline="") as fh:
        for row in csv.DictReader(fh):
            species[row["PROBENNUMMER"]] = row.get("species", "").strip() or "unknown"
    items = []
    for fa in Path(samples_dir).glob(f"*/genes/{tag}.fa"):
        sid = fa.parent.parent.name
        sp = species.get(sid, "unknown").replace(" ", "-").replace("_", "-")
        items.append((sp, sid, fa))
    for sp, sid, fa in sorted(items):
        for line in fa.read_text().splitlines():
            if line.startswith(">"):
                line = f"{line}_{sp}"
            print(line)


if __name__ == "__main__":
    main()
