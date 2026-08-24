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
cd /home/groups/MaxsonLab/braun/analysis
git clone https://github.com/brauntp/scenicplus-healthy-ref.git
cd scenicplus-healthy-ref
mkdir -p slurm/logs          # SLURM discards output silently if this is absent
```

The repo is public, so this needs no credentials, and a clone preserves the
executable bits on all 15 scripts. Pick up later fixes with `git pull`.

A ZIP download or `rsync` works too, but **both lose the executable bits** —
`sbatch` on a non-executable script fails. After either, run:

```bash
chmod +x 00_inspect/* 01_cistopic/* 02_pair/*.py 03_pipeline/*.sh \
         03_pipeline/*.py 04_db/*.sh 04_db/*.py docs/*.py slurm/*
```

---

## STAGE 1 — Inspect (login node, minutes, READ-ONLY)

**Do this first and stop.** Its output sets every parameter below. Both scripts
open files read-only; `inspect_archr.R` never calls `saveArchRProject`.

The reference lives here (confirmed on ARC):

```
$REF = /home/groups/MaxsonBraunScratch/worme/projects/scATAC/
       251112_hematopoiesis_ref/integration/output/02
    rna.h5ad   atac.h5ad   combined.h5ad
    combined_glue_embeddings.tsv
    atac_metadata_with_transferred_labels.tsv
    John_Dick_cell_type_palette.tsv          final.dill
```

```bash
REF=/home/groups/MaxsonBraunScratch/worme/projects/scATAC/251112_hematopoiesis_ref/integration/output/02

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

Two things about this layout to settle from the reports before going further.

**The GLUE embedding and the transferred labels may live in the TSVs rather
than inside the `.h5ad` objects.** `combined_glue_embeddings.tsv` and
`atac_metadata_with_transferred_labels.tsv` existing alongside the objects
suggests exactly that. The inspector reports which `obsm` keys and `obs`
columns are actually present; stage 5 accepts either source
(`--latent-tsv` / `--obs-tsv`).

**`atac.h5ad` is probably not a peak matrix.** It was the GLUE input, so it is
likely a gene-activity or tile matrix — SCENIC+ needs peaks. The inspector's
region-width report settles it. Either way stage 4a exports the PeakMatrix from
the ArchRProject, which is the authoritative peak source; `atac.h5ad` is used
only for its cell ids and latent coordinates.

## First thing in every new shell

```bash
source setenv.sh
```

This sets `$REF` and the pairing parameters, and verifies the three input files
exist. Do it on each login — a variable you exported at a prompt is gone when
you log out, and the next login node inherits nothing. An unset `$REF` makes
`"$REF/rna.h5ad"` expand to `/rna.h5ad`, which used to fail deep inside h5py;
the scripts now catch that and say so, but sourcing the file avoids it entirely.

`slurm/pairing.sbatch` sources it too, since a batch job does not inherit the
submitting shell's variables either. Edit `setenv.sh` rather than the sbatch.

## Which account and partition?

```bash
bash slurm/whoami_slurm.sh
```

Read-only. It reports what your own recent jobs used (`sacct` — the most
reliable answer, since those submissions demonstrably worked), which accounts
you're associated with (`sacctmgr`), the partition list with the default
marked, and which partitions have ≥80 GB per node and a walltime cap above
4 hours — both of which the pairing job needs.

**Settled — these follow `maxsonBraunLab/cutTag-pipeline`:**

```bash
conda activate scplus-pairing      # house style: activate BEFORE sbatch
mkdir -p jobs                      # once; SLURM will not create it
sbatch slurm/pairing.sbatch
```

| setting | value | source |
|---|---|---|
| partition | `batch` | `cluster.yaml` `__default__`, and both `run_pipeline_*.sh` wrappers |
| account | *not set* | no lab script specifies one; the default `maxsonlab` applies |
| logs | `jobs/<name>_%j.{out,err}` | `run_pipeline_conda.sh` (`--output=jobs/run_pipeline_%j.log`); `cluster.yaml` uses the nested `jobs/{rule}/{rule}_%j.*` for per-rule jobs |
| conda | activate before `sbatch` | `run_pipeline_conda.sh` header comments |

`slurm/pairing.sbatch` hardcodes `--partition=batch` accordingly, so no flags
are needed. `batch` has 142 nodes at 391 GB with a 36 h cap, against this job's
80 GB and ~4 h.

**Do not rely on the site default partition.** It is `interactive` — 19 nodes,
shared with people waiting at a prompt — which is why the lab's own scripts all
name `batch` explicitly. `sacct` confirms it: 1,534 of 1,535 recorded jobs on
this account ran on `batch`.

The lab's `callpeaks` rule requests 55 G and `diffbind` 40 G, so this job's 80 G
is large for the group's usual pattern but not out of scale for the partition.

Anything on the `sbatch` command line overrides the directives in the script,
which is why `pairing.sbatch` deliberately does not hardcode either one.

## Region sets: the gap I only found when writing the config

`input_data.region_set_folder` is a **required** config entry — SCENIC+ globs
`region_set_folder/<family>/*.bed` and runs one cisTarget/DEM enrichment per
BED — and nothing in this pipeline produced it. The canonical source is
`01_cistopic/run_cistopic.py`, which at this reference's scale means a ~40 GB
Matrix Market export of 393,832 × 163,969 followed by MALLET LDA: hours to days,
Java heap controllable only through `MALLET_MEMORY`.

The config's own comments name `DARs_cell_type` as a valid region-set family
alongside `topics_otsu`. Per-cell-type DARs are computable directly from
`ACC_GEX.h5mu` — 25,323 metacells already grouped — streaming in blocks, one
batch job at 32G (the first attempt at 16G was OOM-killed; see
`docs/MEMORY.md`):

```bash
sbatch slurm/region_sets.sbatch
```

**What this gives up, stated plainly.** Topics are unsupervised co-accessibility
programs; they can span cell types, split one label into sub-programs, or track a
gradient no label names. DARs are label-driven by construction, so any regulatory
program not aligned to `predicted_CellType_Broad` gets no region set — and an
eRegulon can only be found for a program some region set represents. A DAR-only
run finds cell-type-associated eRegulons well and shared or continuous ones
poorly. Region sets are also the runtime driver, so 24 DAR sets is much cheaper
than 100+ topics, which is part of why it finds less.

For a healthy hematopoietic reference where the question is *which TFs drive
which cell type*, DARs are the aligned choice. Topics can be added later as a
second subdirectory; the DAG runs both families.

Validated on fixtures with planted ground truth: perfect recovery (precision and
recall 1.000 on 1,800 planted regions across 6 groups) and zero regions with an
explicit warning under a null — including a tie-heavy null with 85% exact zeros,
since real accessibility is tie-heavy and the significance test has to survive
that.

**Pass `--diagnostics pairing_diagnostics.csv`.** The thin-group filter counts
*independent* observations, not metacells, and the two differ by the oversample
factor. Raw metacell counts here run from 22 (Pro-Monocyte) to 5,058 (Pro-B), so
every group clears a raw threshold of 5 — filtering on them would skip nothing.
The honest denominators are the `independent_metacell_equiv` values the pairing
job already computed:

| group | metacells | independent | `--min-independent 5` |
|---|---|---|---|
| Pro-Monocyte | 22 | 1 | skipped |
| cDC | 44 | 2 | skipped |
| Late GMP | 48 | 3 | skipped |
| Megakaryocyte Precursor | 60 | 5 | kept |
| Stromal | 78 | 9 | kept |

So **three** groups are skipped, not the two I first wrote — and only when the
diagnostics file is supplied. Without it the script falls back to
`metacells / oversample`, which is close but not the same number (Pro-Monocyte
`floor(22/8) = 2` against the diagnostics' 1) and skips only Pro-Monocyte; it
labels itself APPROXIMATE in the output. It refuses outright if the diagnostics
file does not name every group in the object.

Two of the three — Pro-Monocyte (1) and Late GMP (3) — are named in
`docs/QC_RESULT.md`'s thin-support list. cDC (2) is not: that list stops at the
five groups the QC report happened to enumerate, and cDC sits between them by
independent count. Its exclusion here follows the same rule, not a prior
finding.

**Runtime and memory, both measured — and both wrong on the first attempt.**
Measured on the real object (job 10714858): **42m43s, MaxRSS 2.80 GB** against
a 2 h walltime and 32G. Memory came in 14% over my 2.45 GB projection — the
first estimate in this project to land close, and the only one taken from
measuring the actual code path rather than reasoning about it. Runtime was 3×
under: I timed the per-block compute in isolation (2.14 ms per region) and
ignored reading each block from the backed object plus the BH/selection pass.

Two prior versions were killed, for reasons worth keeping:

1. **Walltime.** The first called `scipy.mannwhitneyu` once per group, which
   re-ranks the same block 24 times: 1.84 ms per region *per group*, i.e. 4.8 h
   against the 4 h I had set. Ranking once per block and deriving each group's U
   from its rank sum is 20× faster.
2. **Memory (job 10714844, exit 137).** I sized `--mem` from the block's own
   size (1.89 GB at 20,000 regions) and called it "flat", ignoring what gets
   computed *from* the block. Measured 9.59 GB: `rankdata` returns float64 — 2×
   a float32 block — argsorts int64 internally, and the tie terms sorted a
   second copy. Ranks and ties now come from one stable argsort in the block's
   dtype: 2.45 GB at a 4,000-region default.

Both rewrites keep the statistics exact — ranks identical to `rankdata` (max
|diff| 0.0, ties included) and p-values matching scipy to 5e-08, tie and
continuity corrections in place. Recall on the planted fixture stays 1.000 with
a 0.13% false-positive rate against a 5% FDR target.

## Which env does each step need?

Two environments, and the split is deliberate — the light steps must not be
blocked on the heavy env solving.

| step | needs | env |
|---|---|---|
| `00_inspect/inspect_h5ad_lite.py`, `compare_labels.py` | h5py | none — `run_lite.sh` finds a python that has it |
| `02_pair/attach_atac_labels.py` | h5py | none |
| `02_pair/aggregate_atac_sparse.py` (**including `--dry-run`**) | numpy, scipy, pandas, anndata, mudata | `scplus-pairing` (`03_pipeline/pairing_env.yml`) |
| `01_cistopic/*`, the SCENIC+ Snakemake pipeline | the full stack | `03_pipeline/environment.yml` |

**`--dry-run` still needs the pairing env.** It is cheap in time and memory —
it plans metacells and prints the footprint without aggregating — but it
imports the same libraries as the real run. Create the env first:

```bash
mamba env create -f 03_pipeline/pairing_env.yml   # ~1 min, no compilation
conda activate scplus-pairing
bash 03_pipeline/preflight_pairing.sh             # confirms the interpreter
```

The pairing env is four libraries and solves in seconds; the full SCENIC+ env
is ~900 packages with six source builds. Keep them separate so pairing can run
before the big one exists.

## Building the scenicplus env (the last long-running unknown)

Job 10715336 failed with `EnvironmentNameNotFound: scenicplus` -- the env was
never built. The config check passed, so this is the only thing missing.

```bash
screen -S spenv                           # NOT a bare login shell
bash 03_pipeline/create_env.sh            # runs the preflight, then builds
conda activate scenicplus
scenicplus --help                          # smoke test
```

**Use `create_env.sh`, not `mamba env create` directly.** A direct
`mamba env create` fails on this spec, and fails *after* the 234-package conda
solve has already succeeded:

```
pybedtools uses setuptools (...) for installation but setuptools was not found
CondaEnvException: Pip failed
```

The cause is neither pybedtools nor this env file. setuptools 84.0.0
(2026-08-08) removed the bundled `pkg_resources`; pybedtools 0.9.1's `setup.py`
opens with `import pkg_resources` and exits with that message when it fails.

**Which setuptools the build sees depends on which build path pip takes, and
which paths exist depends on the PIP VERSION.** Measured on pybedtools 0.9.1
with `setuptools==80.9.0` and `wheel` both installed in the env:

| pip | result | path taken |
|---|---|---|
| 24.0 | built | legacy `setup.py`, runs in the target env |
| 25.2 | built | legacy `setup.py`, runs in the target env |
| 26.2.1 | **FAILED** | legacy path removed; PEP 517 isolated env fetches setuptools 84 |

So two earlier fixes were both dead ends, and both looked correct when tested
against an older pip:

- `PIP_CONSTRAINT=setuptools<81` — not applied to build dependencies here.
- `setuptools<81` + `wheel` in the conda section — only helps where the legacy
  path still exists, i.e. pip < 26.

`PIP_NO_BUILD_ISOLATION` as an *environment variable* is inert on pip 26 at any
value (tested `1`, `true`, `yes`). Only the command-line flag works.

**The build runs in five stages**, which is why `create_env.sh` drives pip
itself instead of letting `mamba env create` run the yml's pip section:

1. Conda section only — the pip section is stripped from a temporary copy, because letting mamba run it would install everything under isolation.
2. The five sdists with no usable wheel, `--no-build-isolation --no-deps`: pybedtools 0.9.1, pyranges 0.0.111, pyrle 0.0.39, tspex 0.6.3, MACS2 2.2.9.1. Their build dependencies come from the conda section, already installed — which is why `setuptools<81` and `wheel` are conda dependencies, and why the script refuses to start without them.
3. The 34 pinned PyPI requirements, isolation **on**. This establishes the pins before anything else can override them.
4. The five git packages, `--no-deps`.
5. `pip check` names every unsatisfied requirement; those get installed, then the script verifies the pinned versions survived.

### Why stage 4 needs `--no-deps`

Resolving the git group normally fails outright:

```
ERROR: Cannot install ... because these package versions have conflicting dependencies.
    The user requested loomxpy 0.4.2 (from git+...LoomXpy@<commit>)
    scenicplus 1.0a2 depends on loomxpy 0.4.2 (from git+...LoomXpy@main)
ResolutionImpossible
```

scenicplus v1.0a2's own `requirements.txt` declares `loomxpy @
git+...LoomXpy@main`, and we pin a commit. **pip treats two different URLs for
one name as irreconcilable even when the versions match** — no version range
fixes it. Note which package is at fault: pycisTopic declares `loomxpy` with no
URL and resolves against our pin fine (verified locally); only scenicplus pins
the moving ref.

There is a second, independent reason. Resolving pycisTopic *with* its declared
dependencies silently overrides our pins — measured on pip 26.2.1:

| package | our pin | resolver picked |
|---|---|---|
| scanpy | 1.8.2 | 1.11.0 |
| anndata | 0.10.5.post1 | 0.11.4 |
| scikit-learn | 1.3.2 | 1.5.2 |
| scipy | 1.12.0 | 1.17.1 |
| numba | 0.59.0 | 0.67.0 |
| polars | 0.20.13 | 1.44.0 |
| pyranges | 0.0.111 | 0.0.127 |
| mudata | 0.2.3 | dropped |

Those pins are transcribed from scenicplus's own `requirements.txt` lock, so
losing them loses the thing this env file exists to reproduce.

Blanket `--no-deps` is **not** the answer either: the git packages declare 196
requirements our pinned list does not cover, and skipping them leaves an env
that imports but breaks at use. Hence stage 5 — `pip check` reports exactly
what is unsatisfied, so the repair is measured rather than guessed, and it runs
inside the pins stage 3 established.

Corrected along the way: `polars` 0.20.13 ships a `cp38-abi3` manylinux wheel
that works on 3.11, and `pyBigWig` 0.3.22 has a `cp311` manylinux wheel — both
were briefly miscounted as source builds because the check looked only for
`cp311` tags. The source-built set is exactly five.

Allow an hour or more: 234 conda packages (513 MB) plus 39 pip requirements.
Two distinct groups of five, easy to confuse: five *sdists* compile from source
(pybedtools, pyranges, pyrle, tspex, MACS2 — stage 2), and five *git*
dependencies are cloned and built at pinned refs (pycistarget, pycisTopic,
LoomXpy, pySCENIC, scenicplus — stage 4).

The preflight is worth the seconds. It checks the platform, solver, python pin,
pandas placement, git pin *reachability* and disk before the solve -- so a
deleted tag or an unreachable repo fails immediately rather than an hour in,
after the conda side has already succeeded. Verified 2026-08-24 for linux-64:
the conda specs solve to 234 packages with no conflicts, and all five git
dependencies resolve at their pinned refs.

`slurm/scenicplus.sbatch` now prints these instructions itself when activation
fails, and distinguishes "env does not exist" from "env exists but will not
activate".

### Fastest path: the h5py-only inspector

`.h5ad` is documented HDF5, so the facts stage 1 needs can be read without
anndata — and therefore without finding whichever env the integration ran in.
`run_lite.sh` locates a python that has h5py and re-executes with it:

```bash
REF=/home/groups/MaxsonBraunScratch/worme/projects/scATAC/251112_hematopoiesis_ref/integration/output/02
bash 00_inspect/run_lite.sh "$REF/rna.h5ad" "$REF/atac.h5ad" \
    --tsv "$REF/combined_glue_embeddings.tsv" \
          "$REF/atac_metadata_with_transferred_labels.tsv" \
    --out report_lite
```

It answers: cell and feature counts; whether ATAC `var_names` are peaks
(with width and `chr` prefix) or gene/tile names; whether the GLUE latent is in
`obsm` or only in the TSV; the `obs` columns and their levels; the cell-id
format — compared across the `.h5ad` files and the TSVs, which is how the
ArchR-vs-GLUE id mismatch surfaces now rather than at stage 5; and whether `X`
is raw counts or normalized.

If no python on the machine has h5py, it prints the one-line fix
(`pip install --user h5py` — binary wheel, seconds).

### Full inspectors (more detail, more dependencies)

**Neither runs in the base conda env** — `Rscript` is usually not on
PATH, and base has no `anndata`. You already have both somewhere: the GLUE
integration ran in a python env with anndata, and the ArchRProject was built in
an R with ArchR. Find them rather than building anything:

```bash
bash 00_inspect/find_inspect_env.sh          # read-only probe
# bash 00_inspect/find_inspect_env.sh --create   # only if nothing suitable exists
```

It walks every conda env for `anndata`/`mudata`, checks whether `Rscript` has
ArchR, and lists candidate `module load R/...` lines when R is absent. Then call
the interpreter it names directly:

```bash
/path/to/glue_env/bin/python 00_inspect/inspect_anndata.py --rna ... --atac ... --out report_py
```

`inspect_anndata.py` needs anndata (mudata only for `--mudata` input);
`inspect_archr.R` needs ArchR (jsonlite optional).

**If the R side is a hassle, do not block on it.** The python inspector answers
most of stage 1 — cell counts, where the latent lives, which label columns
exist, the cell-id format, and matrix provenance. Send that report first; the
ArchR one mainly adds peak geometry, which stage 2 also needs but stage 5 does
not.

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
#   --dry-run to preview. Do NOT pass --skip-scores on a stock checkout: DEM's
#   outputs are unconditional targets in the Snakefile, so the pipeline would
#   fail at the DEM rule after cisTarget had already spent its hours.

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
# Check first -- 10 seconds, and it catches the failures that otherwise appear
# minutes into a ~900-package install with a misleading solver trace.
bash 03_pipeline/preflight_env.sh          # must exit 0

# Pipeline env. Git deps are pinned to explicit commits, so this is reproducible.
mamba env create -f 03_pipeline/environment.yml
conda activate scenicplus
scenicplus --help                                     # must print CLI usage
python -c "import pycisTopic, pycistarget, scenicplus; print('ok')"
pip check                                             # worth the 5 seconds
```

Put the env on a filesystem visible to compute nodes, and set `CONDA_PKGS_DIRS`
somewhere with several GB free — not a quota'd `$HOME`.

**`mallet` is linux-64 only.** conda-forge ships exactly one build
(`mallet-2.0.8-ha770c72_0`, depends on `openjdk`) and no macOS build at all. ARC
is linux-64 so this resolves there; on a Mac the solve fails with *"mallet =\* \*
does not exist (perhaps a typo or a missing channel)"*, which looks like a
channel problem and isn't. The preflight checks this for your platform, and the
env file documents the upstream-tarball fallback.

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
    --rna                 "$REF/rna.h5ad" \
    --atac                cistopic_out/imputed_accessibility.h5ad \
    --latent-key          X_glue \
    --group-key           <CELLTYPE_OBS_KEY> \
    --cells-per-metacell  50 \
    --min-cells-per-group 50 \
    --out                 ACC_GEX.h5mu \
    --diagnostics         pairing_diagnostics
```

**If stage 1 showed the latent or the labels are NOT in the objects**, point at
the TSVs instead — a single combined embedding file covering both modalities is
handled (rows are matched by cell id, non-numeric columns ignored):

```bash
python 02_pair/glue_metacells.py \
    --rna                 "$REF/rna.h5ad" \
    --atac                cistopic_out/imputed_accessibility.h5ad \
    --latent-tsv          "$REF/combined_glue_embeddings.tsv" \
    --obs-tsv             "$REF/atac_metadata_with_transferred_labels.tsv" \
    --group-key           <CELLTYPE_OBS_KEY> \
    --cells-per-metacell  50 \
    --out                 ACC_GEX.h5mu \
    --diagnostics         pairing_diagnostics
```

Both paths refuse rather than silently subsetting: a cell present in the object
but absent from the embedding TSV is a hard error, as is zero cell-id overlap
with a metadata TSV. Every argument is recorded under `params` in
`pairing_diagnostics.json`, so which embedding was used stays traceable.

One caveat specific to this dataset: the ATAC cell ids in
`imputed_accessibility.h5ad` come from the ArchR export, while the latent
coordinates come from the GLUE run. If ArchR wrote `Sample#BARCODE` and GLUE
saw plain barcodes, the ids will not match and the script will say so. Fix it by
harmonising the ids, never by dropping the non-matching cells.

**Read `pairing_diagnostics.csv` before going further.** The
`median_crossmodal_gap` column is the one that matters: a large value for a cell
type means RNA and ATAC do not co-occupy that region of the latent space, so
pairing there is extrapolation and any eRegulon resting on it is weak.

**On `--cells-per-metacell`, from `docs/pairing_sensitivity.csv`:** k=50 gives
the highest median ρ on true links at *every* noise level tested (0.599, 0.508,
0.284, 0.144, 0.090, 0.041 as noise goes 0.15 → 2.5), and the highest AUROC at
four of six. Larger metacells average away more modality-specific noise, so
there is no noise regime in this sweep where a smaller k recovers more signal.
The AUROC ordering does invert at noise ≥ 1.5 (k=10 edges k=50 at 1.5), but the
spread across k there is only 0.058–0.069 while every k is already badly
degraded — that is run-to-run scatter in a regime you should not be operating
in, not a reason to pick a smaller k.

The real cost of large k is **observation count**, not robustness: on the same
input cells, k=10 → 798 metacells, k=25 → 320, k=50 → 160. Those are the rows
the region-to-gene GBM regresses on, so k=50 buys effect size by spending
degrees of freedom. Start at **k=50** and only drop to 25 if a cell type is too
small to yield enough metacells for a stable fit — which `--min-cells-per-group`
and the diagnostics CSV will tell you. If the diagnostics show a large
cross-modal gap, that is a signal to fix the integration or exclude the cell
type, not to change k.

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

# Config path is optional -- defaults to what make_config.sh writes.
# The lab convention is no --account, and the partition is set in the script.
sbatch slurm/scenicplus.sbatch
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
| `mallet =* * does not exist` at env-create | not a channel problem — mallet is linux-64 only | on linux-64 it resolves; elsewhere use the upstream tarball + `--mallet-path`. Run `preflight_env.sh` |
| pip: `requires a different Python` | `python=3.11` resolved past 3.11.8 | the patch pin is load-bearing; keep `python=3.11.8` |
| pybiomart connection error | compute node has no internet | run `--until download_genome_annotations` on a login node |
| walltime kill mid-GBM | `--time` too low | raise it; completed outputs are reused on resubmit |
