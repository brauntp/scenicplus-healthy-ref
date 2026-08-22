# Cluster runbook (ARC)

Copy-paste order of operations. Every command is one of: **login node** (cheap,
minutes) or **sbatch** (queued). Nothing here needs the analysis sandbox.

Placeholders to fill once, then reuse: `<ARC_ACCOUNT>`, `<ARC_PARTITION>`,
`<N_CPU>`, `<MEM>`, `<HH:MM:SS>`, `<CONDA_BASE>`, `<NODE_SCRATCH_ROOT>`.
Find the first two with `sacctmgr show assoc user=$USER format=account,partition`
and `sinfo -s`; `<CONDA_BASE>` is the output of `conda info --base`.

---

## 0. Clone and make the log directory

```bash
git clone https://github.com/brauntp/scenicplus-healthy-ref.git
cd scenicplus-healthy-ref
mkdir -p slurm/logs          # SLURM discards output silently if this is absent
```

---

## STAGE 1 — Inspect (login node, minutes, READ-ONLY)

**Do this first and stop.** Its output sets every parameter below. Both scripts
open files read-only; `inspect_archr.R` never calls `saveArchRProject`.

```bash
# ArchR side: peaks, matrices, fragments; also writes the consensus peak BED
Rscript 00_inspect/inspect_archr.R \
    --proj /path/to/ArchRProject \
    --out  report_archr

# RNA + ATAC AnnData side
python 00_inspect/inspect_anndata.py \
    --rna  /path/to/RNA.h5ad \
    --atac /path/to/ATAC.h5ad \
    --out  report_py
# single integrated object instead:
#   python 00_inspect/inspect_anndata.py --mudata integrated.h5mu --out report_py
```

`inspect_archr.R` needs only ArchR (jsonlite optional). `inspect_anndata.py`
needs anndata + mudata — the scGLUE env from stage 3 is fine, or any env with
those two.

**Send back `report_archr.md`, `report_py.md` and both `.json` files.** They
answer: genome build; whether peaks are fixed-width 501 bp; which `obs` column
works as `key_to_group_by` and whether its levels agree across modalities;
whether the GLUE latent is in `obsm`; which cell types are too small or
single-modality. Those set the pairing parameters and the memory plan.

---

## STAGE 2 — Database decision (login node, free)

Cheap and worth doing before anything expensive. `stage 1` wrote
`report_archr_consensus_peaks.bed`.

```bash
# Pull the 1.8M SCREEN region IDs out of the remote feather's Arrow footer
# (~121 MB range request, not the 33 GiB file). Needs outbound HTTPS.
python 04_db/fetch_db_regions.py --out screen_db_regions.parquet

# How much of YOUR peak set can that database represent?
python 04_db/peak_overlap_audit.py \
    --peaks       report_archr_consensus_peaks.bed \
    --db-regions  screen_db_regions.parquet \
    --out-prefix  audit_screen
```

Read the VERDICT block. Then either:

```bash
# (a) precomputed is good enough -> download it (~33 GiB + ~13 GiB, resumable)
bash 04_db/download_precomputed_db.sh --dest /path/to/db
#   --dry-run to preview, --skip-scores if you only need cisTarget (not DEM)

# (b) dropout too high -> build a custom DB on your own peaks
sbatch --account=<ARC_ACCOUNT> --partition=<ARC_PARTITION> \
       04_db/slurm_build_db.sbatch
```

Note: the ~150–350 bp SCREEN regions being narrower than your 501 bp peaks is
**not** a penalty. pycistarget's filter is `Overlap_query > frac OR
Overlap_target > frac`, and `Overlap_query` is normalized by the *database*
region's width — so a 501 bp peak containing a 272 bp cCRE scores 1.0. The
audit measures what actually matters: peaks with no qualifying cCRE at all.

---

## STAGE 3 — Environments (login node, once, slow)

**Two environments, deliberately.** scGLUE needs torch and a newer
pandas/numpy than `pycisTopic` tolerates; forcing them together breaks the
`pandas==1.5.0` pin. They hand off through the `.h5mu` file, and that
cross-version handoff is tested — a file written by a modern
mudata/anndata stack reads correctly under the pipeline's pinned
`mudata==0.2.3` / `anndata==0.10.5.post1`, including the `.to_df()` calls
`infer_region_to_gene` makes.

```bash
# Pipeline env. Git deps are pinned to explicit commits, so this is reproducible.
mamba env create -f 03_pipeline/environment.yml
conda activate scenicplus
scenicplus --help                                     # must print CLI usage
python -c "import pycisTopic, pycistarget, scenicplus; print('ok')"
pip check                                             # worth the 5 seconds
```

Put the env on a filesystem visible to compute nodes, and set `CONDA_PKGS_DIRS`
somewhere with several GB free — not a quota'd `$HOME`.

Your existing scGLUE env covers stage 5. It needs anndata, mudata,
scikit-learn, numpy, pandas — nothing else.

---

## STAGE 4 — Export + topic modelling (sbatch; the long ATAC step)

```bash
# 4a. Export the PeakMatrix from ArchR (chunked over cells).
# Writes peak_matrix.mtx, barcodes.tsv, peaks.bed, regions.txt,
# cell_metadata.tsv, export_manifest.tsv into --out-dir.
Rscript 01_cistopic/export_from_archr.R \
    --archr-project /path/to/ArchRProject \
    --group-col     <CELLTYPE_OBS_KEY> \
    --out-dir       cistopic_input

# 4b. cisTopic object + LDA + imputed accessibility + region sets
export MALLET_MEMORY=64g       # NOT a Python argument -- see below
python 01_cistopic/run_cistopic.py \
    --matrix            cistopic_input/peak_matrix.mtx \
    --barcodes          cistopic_input/barcodes.tsv \
    --regions           cistopic_input/regions.txt \
    --cell-metadata     cistopic_input/cell_metadata.tsv \
    --group-col         <CELLTYPE_OBS_KEY> \
    --n-topics          "5 10 15 20 30 40 50" \
    --mallet-path       "$(which mallet)" \
    --mallet-memory     64g \
    --n-cpu             <N_CPU> \
    --tmp-path          "$TMPDIR" \
    --region-set-folder region_sets \
    --out-dir           cistopic_out
```

Useful memory flags on `run_cistopic.py`, straight from stage-1 numbers:
`--impute-on-variable-regions` (impute only variable regions),
`--max-impute-gb` (refuse to allocate beyond a ceiling rather than OOM),
`--impute-chunk-size`, `--impute-scale-factor`.

**Two traps here, both verified from source:**

1. `run_cgs_models_mallet` has **no memory keyword**. Heap comes only from
   `MALLET_MEMORY`, exported before the call. A Java `OutOfMemoryError` cannot
   be fixed by any Python argument.
2. `impute_accessibility` returns a **dense** `regions × cells` matrix — about
   93 GB at 100k cells × 250k peaks, ~186 GB with the transpose live. **Read
   `docs/MEMORY.md` before submitting.** The cheapest fix is to impute only the
   regions in your region sets (`selected_regions`); the scalable one is to
   impute and pair per cell type.

Run this under `sbatch` with `--mem` from the `docs/MEMORY.md` table, and point
`--tmp-path` at node-local scratch.

---

## STAGE 5 — Pair (scGLUE env, minutes to an hour)

This is the custom step — the reason this repo exists.

```bash
conda activate <your scGLUE env>
python 02_pair/glue_metacells.py \
    --rna                 /path/to/RNA.h5ad \
    --atac                cistopic_out/imputed_accessibility.h5ad \
    --latent-key          X_glue \
    --group-key           <CELLTYPE_OBS_KEY> \
    --cells-per-metacell  25 \
    --min-cells-per-group 50 \
    --out                 ACC_GEX.h5mu \
    --diagnostics         pairing_diagnostics
```

**Read `pairing_diagnostics.csv` before going further.** The
`median_crossmodal_gap` column is the one that matters: a large value for a cell
type means RNA and ATAC do not co-occupy that region of the latent space, so
pairing there is extrapolation and any eRegulon resting on it is weak. Set
`--cells-per-metacell` from the benchmark in `docs/pairing_benchmark.png` —
k=50 gives the largest effect size, k=25 is more robust when integration is
noisy.

Then gate:

```bash
conda activate scenicplus
python 03_pipeline/validate_h5mu.py ACC_GEX.h5mu   # must exit 0
```

---

## STAGE 6 — Run SCENIC+ (sbatch, 12–24 h)

```bash
cp 03_pipeline/config.template.yaml 03_pipeline/config.yaml
# fill in: cisTopic_obj_fname, GEX_anndata_fname, region_set_folder,
#          ctx_db_fname, dem_db_fname, path_to_motif_annotations,
#          combined_GEX_ACC_mudata <- ACC_GEX.h5mu  (CRITICAL, see below)
#          temp_dir, key_to_group_by

# Always dry-run first: confirms the DAG SKIPS prepare_GEX_ACC
03_pipeline/run_pipeline.sh --config 03_pipeline/config.yaml --cores 8 -- -n

sbatch --account=<ARC_ACCOUNT> --partition=<ARC_PARTITION> \
       slurm/scenicplus.sbatch 03_pipeline/config.yaml
```

**The critical line is `combined_GEX_ACC_mudata: ACC_GEX.h5mu`.** That is the
output slot of `prepare_GEX_ACC`; putting our paired file there makes Snakemake
treat the rule as satisfied. In the dry-run plan, `prepare_GEX_ACC` must NOT
appear. If it does, Snakemake will overwrite our GLUE-paired metacells with
randomly-paired ones — silently, with no error — and every downstream link
collapses to cell-type resolution. This is why `run_pipeline.sh` pins
`--rerun-triggers mtime`; do not call `snakemake` directly.

**If compute nodes have no outbound internet**, `download_genome_annotations`
calls pybiomart and will fail. Run that one rule on a login node first:

```bash
03_pipeline/run_pipeline.sh --config 03_pipeline/config.yaml --cores 4 \
    -- --until download_genome_annotations
```

---

## Outputs

| file | what it is |
|---|---|
| `eRegulon_direct.tsv` / `eRegulons_extended.tsv` | the eGRNs — TF → region → gene |
| `region_to_gene_adj.tsv` | **the peak-to-gene linkages you asked for** |
| `tf_to_gene_adj.tsv` | TF → gene adjacencies |
| `AUCell_direct.h5mu` | per-metacell eRegulon activity |
| `scplusmdata.h5mu` | everything assembled |

Report `pairing_diagnostics.csv` and the audit's dropped-peak fraction in
methods alongside any of these. Both bound what the eGRNs can support.

---

## If it dies

| symptom | cause | fix |
|---|---|---|
| `Killed`, no traceback | OOM killer | `sacct -j <id> -o MaxRSS`; raise `--mem`; see `docs/MEMORY.md` |
| Java `OutOfMemoryError` | MALLET heap | raise `MALLET_MEMORY`; no Python arg will help |
| `prepare_GEX_ACC` in the dry-run plan | `combined_GEX_ACC_mudata` not pointed at `ACC_GEX.h5mu` | fix the config; do not proceed |
| `bedtools not found` | runtime dep missing | it is in `environment.yml`; confirm the env is active |
| pybiomart connection error | compute node has no internet | run `--until download_genome_annotations` on a login node |
| walltime kill mid-GBM | `--time` too low | raise it; completed outputs are reused on resubmit |
