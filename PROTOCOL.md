# Protocol — reproducing this analysis by hand

The linear path, start to finish. Every command, in order, with the check that
tells you it worked before you spend the next block of cluster time.

This document deliberately contains **no troubleshooting**. When something
breaks, `RUNBOOK.md` has the diagnosis and the fix for every failure this
project actually hit. Come back here once it is running again.

**What it produces:** peak-to-gene regulatory linkages and eGRNs (TF → region →
gene) from an unpaired scRNA + scATAC healthy hematopoietic reference,
co-embedded with scGLUE.

**What it costs:** ~2 days wall-clock, most of it in two jobs (topic modelling,
then the SCENIC+ DAG at 12–24 h).

---

## Before you start

Settle these seven values. Step 1 tells you all of them; nothing downstream
works until they are decided.

| value | where it comes from | used in |
|---|---|---|
| `$REF` | path to the reference directory | all steps |
| ArchRProject path | your ArchR output | steps 1, 4 |
| `<CELLTYPE_OBS_KEY>` | an `obs` column present in **both** modalities | steps 4, 5, 6 |
| latent source | `X_glue` in the object, **or** the embeddings TSV | step 5 |
| label source | `obs` column, **or** the transferred-labels TSV | step 5 |
| database choice | measured dropout from step 2 | steps 2, 6 |
| `<N_CPU>`, `--mem` | your partition limits, sized by `docs/MEMORY.md` | steps 4, 6 |

Two environments, built once in step 3: `scenicplus` (the full stack) and a
light pairing env. They hand off through a single `.h5mu` file.

**Run `source setenv.sh` first in every new shell.** It sets `$REF` and the
pairing parameters and verifies the inputs exist. Batch jobs do not inherit
your login shell's variables; `slurm/pairing.sbatch` sources it too.

---

## Step 0 — Get the code

```bash
cd /home/groups/MaxsonLab/braun/analysis
git clone https://github.com/brauntp/scenicplus-healthy-ref.git
cd scenicplus-healthy-ref
mkdir -p slurm/logs          # SLURM discards output silently if this is absent
```

**✓ Check:** `ls slurm/logs` succeeds, and `test -x 03_pipeline/run_pipeline.sh`
returns 0. (A ZIP or `rsync` copy loses executable bits — clone instead.)

---

## Step 1 — Inspect the inputs

**Read-only, login node, minutes.** Do this first and stop. Its output sets
every parameter above.

```bash
source setenv.sh

# ArchR side: peaks, matrices, fragments; also writes the consensus peak BED
Rscript 00_inspect/inspect_archr.R \
    --archr-project /path/to/ArchRProject \
    --out           report_archr

# RNA + ATAC AnnData side
python 00_inspect/inspect_anndata.py \
    --rna  "$REF/rna.h5ad" \
    --atac "$REF/atac.h5ad" \
    --out  report_py
```

**✓ Check:** `report_archr_consensus_peaks.bed` exists (step 2 needs it), and
the reports name a cell-type column present in both modalities.

**Decide from the reports:** whether the GLUE latent lives in the objects
(`obsm["X_glue"]`) or in `combined_glue_embeddings.tsv`; same for the labels.
Step 5 accepts either. Also note which cell types are too small or appear in
only one modality.

---

## Step 2 — Choose the motif database

**Login node, free, minutes.** Region IDs must match your peaks, so this is a
measurement, not a preference.

```bash
# Read the 1.8M SCREEN region IDs from the remote feather's Arrow footer
# (~121 MB range request, not the 33 GiB file)
python 04_db/fetch_db_regions.py --out screen_db_regions.parquet

# How much of YOUR peak set can that database represent?
python 04_db/peak_overlap_audit.py \
    --peaks       report_archr_consensus_peaks.bed \
    --db-regions  screen_db_regions.parquet \
    --out-prefix  audit_screen
```

Read the VERDICT block, then take one branch:

```bash
# (a) dropout acceptable -> download precomputed (~33 GiB + ~13 GiB, resumable)
bash 04_db/download_precomputed_db.sh --dest /path/to/db

# (b) dropout too high -> build a custom database on your own peaks
sbatch --account=<ARC_ACCOUNT> --partition=<ARC_PARTITION> \
       04_db/slurm_build_db.sbatch
```

**✓ Check:** both feathers present (rankings **and** scores). The DEM rule is an
unconditional target in the Snakefile, so a scores-only-skipped database fails
the run hours in, after cisTarget has already spent its time.

**Record the dropped-peak fraction.** It bounds what motif enrichment can find
and belongs in your methods.

---

## Step 3 — Build the environments

**Login node, once, slow.** Put the env on a filesystem compute nodes can see,
and point `CONDA_PKGS_DIRS` somewhere with several GB free — not a quota'd
`$HOME`.

```bash
bash 03_pipeline/preflight_env.sh          # 10 s; must exit 0 before you invest

mamba env create -f 03_pipeline/environment.yml
conda activate scenicplus
scenicplus --help                                      # must print CLI usage
python -c "import pycisTopic, pycistarget, scenicplus; print('ok')"
pip check

# Light pairing env (~1 min, no compilation)
mamba env create -f 03_pipeline/pairing_env.yml
conda activate scplus-pairing
bash 03_pipeline/preflight_pairing.sh
```

**✓ Check:** all four commands above exit 0.

---

## Step 4 — Export peaks and run topic modelling

**sbatch. The long ATAC step.** Size `--mem` from `docs/MEMORY.md` before
submitting and point `--tmp-path` at node-local scratch.

```bash
conda activate scenicplus

# 4a. Export the PeakMatrix from ArchR (chunked over cells)
Rscript 01_cistopic/export_from_archr.R \
    --archr-project /path/to/ArchRProject \
    --group-col     <CELLTYPE_OBS_KEY> \
    --out-dir       cistopic_input

# 4b. cisTopic object + LDA + imputed accessibility + region sets
export MALLET_MEMORY=64g          # heap comes ONLY from this variable
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

**✓ Check:** `cistopic_out/imputed_accessibility.h5ad` exists and `region_sets/`
is non-empty.

`impute_accessibility` returns a **dense** `regions × cells` matrix (~93 GB at
100k cells × 250k peaks). If that does not fit, use
`--impute-on-variable-regions` and `--max-impute-gb` rather than raising `--mem`
indefinitely.

---

## Step 5 — Pair the modalities

**The custom step — the reason this repo exists.** Minutes to an hour.

Stock SCENIC+ pairs unpaired data by drawing RNA and ATAC cells *independently*
within a cell-type label, which destroys every within-cell-type program. This
step instead anchors metacells on the shared GLUE latent space so the paired
cells come from the same point in the manifold. See README for the benchmark.

```bash
conda activate scplus-pairing

python 02_pair/glue_metacells.py \
    --rna                 "$REF/rna.h5ad" \
    --atac                cistopic_out/imputed_accessibility.h5ad \
    --latent-key          X_glue \
    --group-key           <CELLTYPE_OBS_KEY> \
    --cells-per-metacell  50 \
    --min-cells-per-group 50 \
    --out                 ACC_GEX.h5mu \
    --diagnostics         pairing_diagnostics
```

If step 1 showed the latent and labels live in the **TSVs**, swap those two
flags (everything else is identical):

```bash
    --latent-tsv          "$REF/combined_glue_embeddings.tsv" \
    --obs-tsv             "$REF/atac_metadata_with_transferred_labels.tsv" \
```

Then gate before committing to a multi-hour job:

```bash
conda activate scenicplus
python 03_pipeline/validate_h5mu.py ACC_GEX.h5mu       # must exit 0
```

**✓ Check:** `validate_h5mu.py` exits 0, **and** you have read
`pairing_diagnostics.csv`. The `median_crossmodal_gap` column is the one that
matters: a large value for a cell type means RNA and ATAC do not co-occupy that
region of the latent space, so pairing there is extrapolation and any eRegulon
resting on it is weak.

**On `--cells-per-metacell`:** use 50. It gives the highest median ρ on true
links at every noise level in `docs/pairing_sensitivity.csv`. The cost is
observation count, not robustness — k=10/25/50 yield 798/320/160 metacells for
the region-to-gene GBM to regress on. Drop to 25 only when a cell type is too
small for a stable fit.

---

## Step 6 — Run SCENIC+

**sbatch, 12–24 h.**

```bash
cp 03_pipeline/config.template.yaml 03_pipeline/config.yaml
```

Fill in: `cisTopic_obj_fname`, `GEX_anndata_fname`, `region_set_folder`,
`ctx_db_fname`, `dem_db_fname`, `path_to_motif_annotations`, `temp_dir`,
`key_to_group_by`, and — critically —
`combined_GEX_ACC_mudata: /absolute/path/to/ACC_GEX.h5mu`.

```bash
# ALWAYS dry-run first
03_pipeline/run_pipeline.sh --config 03_pipeline/config.yaml --cores 8 -- -n

sbatch slurm/scenicplus.sbatch
```

**✓ Check — this is the one that matters most:** `prepare_GEX_ACC` must **NOT**
appear in the dry-run plan. That rule's output slot is
`combined_GEX_ACC_mudata`; pointing it at your paired file is what makes
Snakemake treat the rule as already satisfied. If the rule appears, Snakemake
will overwrite your GLUE-paired metacells with randomly-paired ones — silently,
no error — and every downstream link collapses to cell-type resolution. Fix the
config before submitting; do not proceed.

This is also why you launch through `run_pipeline.sh`, which pins
`--rerun-triggers mtime`. Do not call `snakemake` directly.

If compute nodes have no outbound internet, run the annotation rule on a login
node first:

```bash
03_pipeline/run_pipeline.sh --config 03_pipeline/config.yaml --cores 4 \
    -- --until download_genome_annotations
```

---

## Step 7 — Read the results

The final object is ~41 GB and ~99% of it is the input matrices carried
through. The biology is in the two eRegulon TSVs. Do not transfer the h5mu.

```bash
conda activate scplus-pairing

# What did the pipeline find? Reads the two eRegulon TSVs, not the 41 GB object.
python 05_report/summarize_eregulons.py \
    --direct    03_pipeline/eRegulon_direct.tsv \
    --extended  03_pipeline/eRegulons_extended.tsv \
    --out-prefix docs/eregulon_summary

# The peak-to-gene links, filtered and written as TSV + BEDPE
python 05_report/export_peak_gene_links.py \
    --adj     03_pipeline/region_to_gene_adj.tsv \
    --out-dir 05_report/plot_bundle

# Assemble everything a local plotting session needs, then print the rsync line
bash 05_report/make_plot_bundle.sh
```

(Every path above is that script's default, so the flags can be dropped if you
kept the standard layout.)

**✓ Check:** `05_report/plot_bundle/` exists and contains the peak-gene links.
Expect a few hundred MB — if it is tens of GB, something copied the h5mu.

### Outputs

| file | what it is |
|---|---|
| `region_to_gene_adj.tsv` | **the peak-to-gene linkages** |
| `eRegulon_direct.tsv` / `eRegulons_extended.tsv` | the eGRNs — TF → region → gene |
| `tf_to_gene_adj.tsv` | TF → gene adjacencies |
| `AUCell_direct.h5mu` | per-metacell eRegulon activity |
| `scplusmdata.h5mu` | everything assembled (41 GB; leave it on the cluster) |

---

## The checkpoints, in one place

Each gate is cheap and sits in front of something expensive. In order:

| after step | check | why it is here |
|---|---|---|
| 1 | consensus peak BED written; group key present in both modalities | everything downstream is parameterized from this |
| 2 | both feathers present; dropped-peak fraction recorded | DEM is unconditional — a partial DB fails hours in |
| 3 | `preflight_env.sh`, `scenicplus --help`, `pip check` all exit 0 | catches env breakage in 10 s instead of mid-solve |
| 4 | `imputed_accessibility.h5ad` + non-empty `region_sets/` | the long ATAC step actually produced its outputs |
| 5 | `validate_h5mu.py` exits 0; `pairing_diagnostics.csv` read | bounds what any eRegulon can support |
| 6 | **`prepare_GEX_ACC` absent from the dry-run plan** | otherwise the custom pairing is silently discarded |

---

## Reporting this analysis

Report alongside any result, because both bound what the eGRNs can support:

- the cross-modal gap per cell type from `pairing_diagnostics.csv`
- the dropped-peak fraction from the step-2 audit

And state the limits plainly. A region-to-gene link here means accessibility and
expression covary across the co-embedded manifold — **not** that they were
observed in the same nucleus. If the integration mislocates a population,
pairing mixes the wrong cells and nothing downstream flags it. A healthy
reference is a baseline, not a contrast; it cannot support "absent in disease"
claims in either direction. README has the full list.

---

## When something breaks

`RUNBOOK.md` — every failure this project hit, with the diagnosis and the fix:
OOM kills and the measured memory figures, MALLET heap, the cisTarget and AUCell
fork failures, env and pin drift, BioMart, resuming a partial DAG.

To check the documentation itself still matches the code:

```bash
bash docs/fetch_snakefile.sh    # pinned Snakefile, checksum-verified
python docs/verify_claims.py    # re-derives every number in these docs
```

`docs/PIPELINE_DAG.md` maps all 13 Snakefile rules — useful when a run fails and
you need to know what the failing rule consumes.
