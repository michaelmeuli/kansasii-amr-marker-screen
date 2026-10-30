# CLAUDE.md

See README.md. Easy to get wrong:

- **Environment:** `env_amr` (`environment.yml`); system python has no samtools/mafft/blast.
  Type check with `~/data/conda/envs/env_typecheck/bin/mypy` (strict, py3.11).
- **Paths:** `env.sh` (shell) and the `AMR`/`ROOT` constants (python) derive
  `output/amr` (`OUTPUT_ROOT/amr`; override `AMR_DIR`). refs/ and work/ are NOT in the repo.
  `config/targets.tsv` / `markers.tsv` / `species_refs.tsv` are resolved relative to the script, not the cwd.
- **Paths:** `scripts/kansasii_paths.py` is a vendored copy of `kansasii-docs/templates/kansasii_paths.py` (never edit it here; `bash repos/check_sync.sh --fix`); the shell steps use `scripts/env.sh`, which sources the shell twin `kansasii_paths.sh`.
- **Real runs go through `sbatch submit_amr.sbatch`**, not the login node. Input assemblies: `ASSEMBLIES/*/1_unicycler/*.fasta` (`runs/mkan329/assembly/results/`).
- **No hardcoded coordinates:** everything comes from the GFF in `02_make_targets.py`.
- **Marker numbering is external** (E. coli rRNA, Mtb/M. avium protein); `05_annotate.py`
  transfers it by alignment. Spot-check `mkn_position` after any reference change.
- **`NOCALL` is not wild-type** (depth-masked N or gap). Research use only, not clinical.
- `config/markers.tsv` / `targets.tsv` / `species_refs.tsv` are tab-separated; keep them that way.
