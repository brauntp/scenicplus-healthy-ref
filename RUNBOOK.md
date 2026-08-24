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

### The first real build drifted five pins — and why

Stage 5's repair ran **unconstrained**, so pip resolved the missing packages
freely and upgraded five pinned versions on the way:

| package | spec | installed |
|---|---|---|
| pandas | 1.5.0 | **3.0.5** |
| matplotlib | 3.6.3 | 3.11.1 |
| polars | 0.20.13 | 1.44.0 |
| statsmodels | 0.14.1 | 0.14.6 |
| tqdm | 4.66.2 | 4.70.0 |

A pandas *major* version jump against a spec that pins 1.5.0 is not a cosmetic
difference. Nothing forced it: the tightest declared floor among the git
packages is tmtoolkit's `pandas>=1.4.0`, which 1.5.0 satisfies — verified that
all five resolve at their pinned versions with pandas 1.5.0 held. "Already
installed" is not a constraint; pip will upgrade an installed package to satisfy
a new one.

The repair step now runs under `PIP_CONSTRAINT` built from the spec's own pins,
so pip either finds a solution within them or fails and says so. Drift is now a
build **failure** rather than a printed note — an env with pandas 3 where the
spec pins 1.5.0 would otherwise fail somewhere inside a multi-hour pipeline run.

Stage 5 also verifies what stage 4 installed `--no-deps`: it imports all five
git packages and checks the `scenicplus` console script responds. On the first
build `scenicplus --help` gave "command not found" while the build reported
success, because nothing checked.

**If an existing env has drifted**, repair it in place rather than rebuilding.
(The first version of this refused with "env 'scenicplus' does not exist" from
*inside an activated scenicplus prompt* — it parsed `micromamba env list`, which
fails silently to stdout when it cannot read `~/.condarc` while still exiting 0.
Env detection is now filesystem-based: `$CONDA_PREFIX` when active, then
`$MAMBA_ROOT_PREFIX/envs`, `<solver> info --base`/envs, and the usual install
roots. It also refuses if the interpreter version does not match the spec's
python pin, so a name collision cannot send a force-reinstall into the wrong
env.)

```bash
bash 03_pipeline/create_env.sh --repair-pins
```

That force-reinstalls the pinned versions under constraint, reinstalls the git
packages, repairs what `pip check` reports within the pins, and re-verifies.

### The cisTarget OOM: the recovery matrix, not the region sets

Equal-N sets at 5,000 regions did not fix it, and the staleness check proved why
that mattered: every output was newer than the newest `.bed`, so the sets were
**already** equal-N when the run died. A 0.6 GB per-worker database slice cannot
explain a 128 GB kill.

The term I had missed is in `cisTarget.run_ctx()` -> ctxcore `recovery()`:

```python
rccs = np.empty(shape=(n_features, rank_threshold))
```

`n_features` = every motif in the database (32,765); `rank_threshold` =
`int(ctx_rank_threshold x 1,837,304)`; dtype float64; then wrapped in `df_rccs`
as a second copy. At the configured `0.05` that is **44.9 GB per worker**, and
joblib runs `n_cpu` of them. Nothing about it scales with region-set size.

Fix it in the config, not the allocation:

```yaml
params_motif_enrichment:
  ctx_rank_threshold: 0.01     # was 0.05
```

9.0 GB per worker, 144 GB at 16 — or `--cpus-per-task 8` for 72 GB. Size it
first:

```bash
python 03_pipeline/size_cistarget_memory.py --region-set-folder 01_atac/region_sets
```

NES, AUC and the enriched-motif calls are **unchanged**: the AUC integrates only
the first `round(ctx_auc_threshold x total) - 1 = 9,186` ranks, which any
threshold above the floor still covers. See `docs/MEMORY.md` for what does change
(the leading-edge search) and why it does not affect motifs that clear NES.

Floor is **0.0055**, not 0.005: scenicplus truncates the curve length with
`int()` while ctxcore rounds the cutoff, so equality fails the assert by one.

### Rewriting the region sets in place does NOT invalidate the results

Both motif-enrichment rules declare the region-set **folder** as their input:

```python
input: region_set_folder=config["input_data"]["region_set_folder"]
```

Snakemake stats that directory. **A directory's mtime changes when an entry is
added or removed — not when a file inside is overwritten in place.** Verified
directly against snakemake:

| change to the region sets | snakemake |
|---|---|
| overwrite an existing `.bed` in place | `Nothing to be done` |
| add a new `.bed` | rerun |

`choose_dar_threshold.py --top-n N --write` rewrites the same filenames. So after
running it, snakemake happily resumes with `dem_results.hdf5` and everything
downstream still computed from the **previous** region sets.

This is worse than a crash. eGRNs would be assembled from cisTarget results on
one region-set definition and DEM results on another, and nothing would say so.

`pipeline_status.py` now compares every output at or downstream of motif
enrichment against the newest `.bed` mtime and prints the `rm` commands for
anything older. Run it after any `--write`:

```bash
python 03_pipeline/pipeline_status.py --config 03_pipeline/config.yaml
```

### A "login-node safe" tool that was not

`size_cistarget_memory.py` was documented as "runs in seconds on a login node."
Its first version counted motifs with

```python
pa.ipc.open_file(src).get_batch(i).num_rows          # primary
pf.read_table(db, columns=[]).num_rows               # "safe" fallback
```

Measured in a fresh process on a 3.04 GB / 400k-column fixture:

| approach | time | peak RSS |
|---|---|---|
| `get_batch(i).num_rows` | 1.43 s | **3341 MB** |
| `read_table(columns=[])` | 0.82 s | **3508 MB** |
| `read_table(columns=[one])` | 0.40 s | 251 MB |
| schema only (gives columns, not rows) | 0.20 s | 148 MB |

Both of the ones I shipped scale with **file size** — they map the whole thing
resident. On the real 33 GB rankings database that is tens of GB on a shared
login node. The fallback was *worse* than the primary path, so there was no safe
route through that function at all.

Fixed two ways. The count now comes from one bounded column read (~131 KB of
data regardless of file width): 0.40 s, 283 MB on the same fixture, a 12×
reduction. And `--db` is now opt-in rather than the documented default, because
the motif count is a fixed property of the motif collection (v10nr_clust =
32,765) and the script never needed the file. The default path — counting lines
in the `.bed` files — measures **0.16 s for 21 sets totalling 1.69 M lines**, so
that was never the slow part.

`slurm/size_cistarget.sbatch` runs it as a batch job for the case where you do
want `--db` against a database on a cold filesystem.

**The general rule this earned:** "read-only" is not "cheap." A single
`read_table` or `mudata.read()` on a login node can pull tens of GB resident, and
the tools in this repo that are genuinely login-node safe say *what* they open —
`peek_h5mu.py` reads HDF5 metadata only, `peak_overlap_audit.py` reads one
column, `pipeline_status.py` stats without opening. Any claim of login-node
safety in this repo should name the bytes it touches.

### Reading a `--keep-going` run: the log's tail lies

The run that got furthest ended with `region_to_gene` logging `Done!` and
`Finished job 9` — and still exited 1, with **two** `oom_kill` events. The tail
looked healthy because `run_pipeline.sh` passes `--keep-going`: a rule dies,
independent branches carry on, and a *later* success prints after the failure.

So `3 of 12 steps (25%) done` is not "25% of the pipeline exists". It counts
steps that **ran in that invocation** — the annotation builder had already
satisfied `download_genome_annotations` (13 rules → 12, confirming it worked),
and several rules from earlier runs were already satisfied and so not counted.
`pipeline_status.py` answers the actual question by stat-ing the outputs.

Three things it now reports, each of which took a measurement:

- **Which rule failed, and whether it was OOM.** The job's `.err` shows the traceback but not always the owning rule; the `.out` shows neither. Snakemake's own log has `Error in rule <name>:` with the shell exit code — verified against a real failing workflow. Exit **137** is SIGKILL, which under SLURM means the OOM killer.
- **Every failure, not the last.** Under `--keep-going` one run can fail several rules. Snakemake also repeats each failure in an end-of-run summary where no exit code precedes it, so a naive scan reports the same rule twice, once with `exit None` — deduplicated, preferring the entry that carries the code.
- **That the heavy rules serialise.** Verified against snakemake directly with two `threads: 16` rules at `--cores 16`: they run one at a time. So `--mem` sizes the largest single rule, not a sum — and the OOM'd rule named above is the one to size for.

Two rules can OOM where `region_to_gene` survived on the same allocation,
because `region_to_gene` slices to one gene's search space (~38 peaks) per task
while `tf_to_gene` regresses every gene against every TF over the full
expression matrix.

### The DAG ran: BioMart failed, and the job OOM-killed

First submission that got past path validation. Two independent failures.

**1. `download_genome_annotations` — not our bug.**

```
xml.etree.ElementTree.ParseError: mismatched tag: line 62, column 2
```

`pybiomart` asked `http://www.ensembl.org` for its dataset configuration and got
malformed XML. BioMart returns HTML error pages and truncated responses under
load, and `pybiomart` feeds the body straight to an XML parser. The rule then
makes a *second* network chain — NCBI esearch plus an assembly report — to derive
chromosome sizes. Two flaky services standing between a 24-hour job and its
first computation, for two small static files.

`03_pipeline/build_genome_annotation.py` builds them instead, from UCSC's API
(one call per chromosome, whole-chromosome responses in under a second; Ensembl
REST caps `/overlap/region` at 5 Mb, which would need ~620 paginated calls):

```bash
python 03_pipeline/build_genome_annotation.py --out-dir 03_pipeline \
    --genome hg38 --check-against ACC_GEX.h5mu
```

46 s, 278,455 protein-coding transcripts over **19,723 genes** — right for
GRCh38. Formats were read from the pinned source, not guessed, and verified:
column names exact, `Strand` as `+`/`-` (not 1/−1), `Transcript_type` filtered to
`protein_coding`, chromsizes `Start` all zero, and TSS = `Start` on `+` /
`End` on `−` for every row. Once both TSVs exist the rule drops out of the DAG.

`--check-against` exists because `get_search_space` matches genes with
`gene_annotation.query("Gene in @scplus_genes")` — a plain intersection that
returns **empty without erroring**. If the RNA object used Ensembl IDs the
search space would be silently empty. Tested both namespaces: 100% match on
symbols, 0% plus a warning on `ENSG…`.

**2. OOM at `motif_enrichment_cistarget` — my sizing error.**

`--mem=54G` came from the `region_to_gene` probe: 43.2 GB measured, +25%. That is
a *different rule*. Reading what the code loads:

| term | figure |
|---|---|
| 393,832 matched regions × 32,765 motifs × int32 | 48.1 GB |
| pyarrow Table + DataFrame both briefly live | **~96 GB peak** |

The good news is that `load_db()` calls `db.load(subset)`, not `load_full()`, and
`ctxcore` pushes the column list into `read_table(columns=...)` — so it reads only
database regions matching our peaks, not all 1.84M (~224 GB). But the audit found
near-total matching, so the subset is most of the file.

The two terms **do not add**: every DB-loading rule declares `threads: n_cpu`, so
at `--cores 16` snakemake serialises them. `--mem` covers the max, now **128G**.

Worth noting the sbatch header already said the database term "is usually the
binding constraint" and the triage text already named this exact OOM. Both were
right; the directive contradicted them because it was set from the rule that had
been *measured*, and a measured number for the wrong rule beat an unmeasured one
for the right rule.

### `combined_GEX_ACC_mudata` must be ABSOLUTE — the third base mismatch

Stages 1 and 2 passed; stage 3 failed with the paired object "not found" at
`03_pipeline/ACC_GEX.h5mu` while it sat at the repo root. `make_config.sh` had
reported the same key `[OK]`.

Two validators, two resolution bases. `make_config.sh` resolved relative
`output_data` paths against the **repo root**; snakemake resolves them against
`--directory`, which `run_pipeline.sh` sets to `dirname(config)` =
`03_pipeline/`. So the check passed against a base the run never uses.

Proven on the real Snakefile with the object at the repo root and
`--directory 03_pipeline`:

| config value | snakemake |
|---|---|
| `"ACC_GEX.h5mu"` (relative) | `MissingInputException` — and it tries to **rebuild** the paired object |
| `"<ABS_PATH>/ACC_GEX.h5mu"` | 13 jobs, `prepare_GEX_ACC` absent |

The relative form is not merely a failed lookup: it puts the run one step from
regenerating GLUE-paired metacells as randomly-paired ones. Three fixes:

1. the key is now absolute in the template — it is the **one** `output_data` entry naming a file that must already exist, so it must not float with the working directory, while the other seventeen are genuine outputs and stay relative by design;
2. `make_config.sh` now resolves `output_data` against snakemake's workdir, so it cannot certify a path the run resolves elsewhere — verified it now reports `MISS` for the relative form and `OK` for the absolute one, matching snakemake in both;
3. `run_pipeline.sh` distinguishes "absent" from "wrong base": if the value is relative and the file exists one level up, it says so and prints the fix. Verified it stays silent when the file is genuinely missing.

### The absent cisTopic object is a safety guard, not a gap

The second submission died on `run_pipeline.sh` demanding
`01_atac/cistopic_obj.pkl` exist, with my own message claiming "snakemake still
needs the path to resolve while building the DAG." **That claim was false**, and
two scripts in this repo disagreed about it: `make_config.sh` called the same key
branch-dead while `run_pipeline.sh` refused to start without it.

Measured against the real pinned Snakefile on a snakemake dry run:

| paired MuData | result |
|---|---|
| present | **13 jobs planned, `prepare_GEX_ACC` absent from the DAG** — its inputs are never evaluated |
| absent | `MissingInputException in rule prepare_GEX_ACC_multiome`, naming both files |

So snakemake does not evaluate a rule's inputs unless it needs the rule. And
the reason it doesn't need it is **not** `is_multiome`: both branches of the
Snakefile conditional (`prepare_GEX_ACC_multiome` / `_non_multiome`) declare the
*same* two inputs, and the flag only selects which variant is defined. What keeps
the rule out of the DAG is that its output already exists.

That inverts the design. **Those files being absent is load-bearing.** If
snakemake ever judged the paired object stale, a missing input aborts the run —
whereas with the files present it would regenerate the object through SCENIC+'s
own label-random pairing, replacing GLUE-paired metacells with randomly paired
ones and reporting success. `--rerun-triggers mtime` makes that unlikely; a
missing input makes it impossible.

Both scripts now say so, and the check is inverted: absent reports `[guard]`,
**present** is what gets flagged. Do not create placeholder files to satisfy a
checker — that disarms the guard.

The exemption is keyed to the **paired MuData's presence**, not to
`is_multiome`. An earlier version of `make_config.sh` reset the exempt set with
`if not multi`, which contradicted its own reasoning: on the non-multiome branch
it forced both keys back to required even with the paired object present.
Confirmed on the real Snakefile that this was wrong — with the mudata present,
both settings plan the same 13 jobs and omit `prepare_GEX_ACC`. Verified across
all four combinations of `is_multiome` × mudata presence; the two guard keys are
exempt whenever the mudata exists and required whenever it does not, and when
they *are* required the output says why (build the paired object rather than
supplying them).

### The first pipeline submission died in one second

`jobs/scenicplus-<id>.out` showed `started` and `finished` at the same
timestamp, `pipeline FAILED (exit 1)`, and no reason — because snakemake and
`run_pipeline.sh` write to stderr, which `#SBATCH --error` sends to the
**`.err`** file. **Read both logs; the `.out` file alone cannot explain a
failure.**

The `.err` file named it:

```
ERROR Snakefile not found: .../scenicplus-healthy-ref/src/scenicplus/src/scenicplus/snakemake/Snakefile
```

Two bugs in one line. The path has a doubled `src/scenicplus/`, and it points
inside this repository — which never contained a scenicplus checkout, so the
default could not resolve under any circumstance. The error's advice compounded
it by telling you to `git clone` the upstream repo, which was never needed: the
Snakefile **ships inside the installed package**, and scenicplus's own
`init_snakemake` command finds it with
`files("scenicplus.snakemake").joinpath("Snakefile")`.

`run_pipeline.sh` now asks the same question the upstream code does — resolving
from the installed package via `importlib.resources`, with `--snakefile` still
overriding. Verified across four cases: package present (resolves and proceeds),
package absent (actionable message naming the import check, not a clone),
explicit path (overrides), explicit but missing (still refuses).

The same doubled path was also asserted in `config.template.yaml`'s provenance
header; those are upstream-repository paths and are now labelled as such.

### `CXXABI_1.3.15 not found` — the smoke test's first real catch

Seven of eleven stages failed at import with:

```
ImportError: /lib64/libstdc++.so.6: version `CXXABI_1.3.15' not found
```

| stage | module | result |
|---|---|---|
| motif_enrichment_cistarget | `pycistarget.motif_enrichment_cistarget` | FAIL |
| motif_enrichment_dem | `pycistarget.motif_enrichment_dem` | FAIL |
| cistarget_io | `pycistarget.input_output` | FAIL |
| cistarget_result | `pycistarget.motif_enrichment_result` | FAIL |
| region_to_gene | `scenicplus.grn_builder.gsea_approach` | FAIL |
| eregulon_assembly | `scenicplus.grn_builder.modules` | FAIL |
| AUCell | `scenicplus.eregulon_enrichment` | FAIL |
| tf_to_gene | `arboreto.algo` | ok |
| cli_entrypoint | `scenicplus.cli.scenicplus` | ok |
| snakemake_driver | `snakemake` | ok |
| data_wrangling | `scenicplus.data_wrangling.adata_cistopic_wrangling` | ok |

Every functional check in the smoke test's section 3 also passed, including the
`pyarrow` feather round-trip — so pyarrow is fine and this is specific to the
`pyranges`-adjacent compiled chain. That split is why a per-stage import matrix
beats a single "does it import" check.

Not a Python version problem. PyPI wheels for the compiled packages in this
stack (`sorted-nearest`, `ncls`, `pyrle` and their kin, reached through
`pyranges`) are built against newer GCC than a RHEL-family login node ships.
The conda section listed `c-compiler` and `cxx-compiler` — the toolchain for
*building* — but not `libstdcxx-ng`, the runtime those built extensions *load*.
With no `libstdc++` under the env prefix, the dynamic linker falls back to the
system `/lib64` one, and that one is too old.

`environment.yml` now lists `libstdcxx-ng` and `libgcc-ng`; verified the conda
section still solves at the same 234 packages on linux-64, and that conda-forge
supplies libstdcxx 16.1.0 — well past the required ABI.

**Having the library is necessary but not sufficient.** The env already
shipped `libstdc++` at CXXABI 1.3.17 — newer than the 1.3.15 required — and the
error was byte-identical anyway. The measurement settled it:

```
as-is: 7/7 failed    with env lib dir: 0/7 failed
```

So the library was present, new enough, and simply not being used. The dynamic
linker resolves each extension's dependencies by `DT_RUNPATH`, then
`LD_LIBRARY_PATH`, then `ld.so.cache`, then the default directories — and the
offending extension's RPATH does not reach `$CONDA_PREFIX/lib`, so it got
`/lib64/libstdc++.so.6` while a newer copy sat in the same env.

The untruncated error path corrected an assumption I had made: the file is under
the env's own `lib/python3.11/lib-dynload`, one of **CPython's bundled extension
modules**, not a pip-installed wheel in `site-packages`. So this is not a
manylinux-wheel-versus-conda-build story, and `diagnose_cxxabi.sh` now scans
`lib-dynload` and the env `lib` directory alongside `site-packages` — its
original site-packages-only scan would have found nothing and said so
misleadingly.

**For an env that already exists**, don't rebuild:

```bash
conda activate scenicplus
bash 03_pipeline/fix_cxxabi.sh     # seconds
```

It installs the packages if absent, then imports the seven affected modules
**twice** — once as-is and once with `$CONDA_PREFIX/lib` prepended — and lets the
result pick the diagnosis rather than assuming one:

1. clean as-is: done, nothing further needed.
2. clean only with the lib dir: it writes `etc/conda/activate.d/zzz_libstdcxx.sh` in the env, so the prepend fires on every `conda activate` **including inside a batch job**. An `export` typed in a login shell would not reach `sbatch`; the hook does. It is guarded against duplicating on repeated activation, since `LD_LIBRARY_PATH` is inherited by children.
3. failing both ways: the search-path explanation is ruled out, so a wheel wants an ABI newer than conda-forge's libstdcxx provides — `03_pipeline/diagnose_cxxabi.sh` section 4 names the offending `.so` and its owning package, and section 5 shows it against scenicplus's lock. Pinning it to the lock version is then both the fix and a step toward the lock.

`slurm/scenicplus.sbatch` carries the same guard independently, for the
house-style path where the env is activated before `sbatch` and inherited: an
inherited env brings `LD_LIBRARY_PATH` only if it was already set at submit
time, so the job asserts it rather than assuming.

Re-running `fix_cxxabi.sh` after the hook is in place is safe and worth doing
once: it takes the "clean as-is" branch and *rewrites* the hook, so a hook
written by an earlier version of the script (carrying the superseded
wheel-RPATH explanation in its comment) is replaced with the correct one. The
hook's behaviour was always right; only its comment was wrong.

This is exactly the failure class flagged in the section below: `sorted-nearest`
is one of the 187 uncovered lock pins.

### The env builds, but our spec covers 35 of scenicplus's 222 pins

After a successful build or repair, pip prints a long list of
`scenicplus 1.0a2 requires X==a, but you have X==b`. That is not a failure and
not noise — it is worth understanding before trusting the env.

scenicplus v1.0a2's `requirements.txt` is a **full lock file: 222 pinned
packages.** `03_pipeline/environment.yml` carries 34 pip pins plus 21 conda
specs, so **187 of the upstream pins are uncovered** and pip resolves those
freely to current versions. In the alphabetical tail visible in one run, 9 of 26
were major-version jumps — `toolz` 0.12→1.1, `zope-interface` 6.2→8.6,
`xmltodict` 0.13→1.0, `url-normalize` 1.4→3.0 — and that slice is roughly a
seventh of the full list.

Two things are true at once: the warnings are pip comparing *declared* against
*installed* and are not evidence of breakage; and they are not evidence of
safety either. Reasoning about 187 packages is not the way to settle it —
exercise the code:

```bash
conda activate scenicplus
bash 03_pipeline/smoke_test.sh      # read-only, minutes
```

It checks the CLI responds, imports one module per Snakemake rule so a failure
*names the stage* that would break hours into a run, and runs the specific
operations SCENIC+ performs on the paired object: building a `MuData`, calling
`.to_df()` on a modality, the gradient-boosting regressor behind TF-to-gene, a
`pyarrow` feather round-trip as cisTarget does on the rankings database, and the
snakemake and ray entry points.

If a stage fails on a package pip warned about, add it at scenicplus's pinned
version to the pip section of `environment.yml` and re-run `--repair-pins`.
Adding all 187 up front would be the *reproducible* choice, but it also
re-introduces the resolver conflicts this file spent three rounds escaping, so
it is worth doing only for pins that demonstrably matter.

### `--repair-pins` v1 deleted packages — read this before using an old checkout

The first version of the repair ran `pip install --force-reinstall --no-deps -r
<all 34 pins>`. That was wrong twice:

- It omitted `--no-build-isolation`, and five of those 34 are the wheel-less sdists. On pip 26 they cannot build with isolation on — the very failure this wrapper exists to avoid — so the step could not succeed as written.
- `--force-reinstall` **uninstalls each package before reinstalling it.** A failure partway through the list leaves the already-processed packages *removed*.

That is what produced the contradictory report: 33 of 34 pins absent *and* the
git packages still importing. The diagnostic settled it — `mudata`, `sklearn`
and `matplotlib` appeared in both the ABSENT list and the import-FAIL list, so
they were genuinely gone, not victims of damaged metadata. The git packages
survived because they were never in the list being force-reinstalled.

**The repair now touches only what is actually wrong.** It reads the installed
versions, computes the drifted set, splits it (sdists get
`--no-build-isolation`), and installs with no `--force-reinstall` at all — a
plain pinned install downgrades an installed package and leaves everything else
alone. The git packages are reinstalled only if one of them fails to import.
Verified on a deliberately damaged venv: the drifted pin downgraded, the absent
sdist built, distribution count went **up**, residual drift empty.

### When the env's own reports contradict each other

`--repair-pins` once reported **"STILL DRIFTED: 33 of 34 ... installed
MISSING"** while, three lines later, all five git packages imported and the
`scenicplus` console script existed. Both cannot be true — scenicplus imports
pandas. So one report was wrong, and the answer is to ask the environment
directly rather than reason about it:

```bash
conda activate scenicplus
bash 03_pipeline/diagnose_env.sh      # read-only, seconds
```

It reports three things the flat "MISSING" list conflated:

1. **How many distributions each mechanism sees** — `importlib.metadata` versus `pip list`. Wildly different counts mean the query is broken, not the env. An interrupted `--force-reinstall` can leave partially-written `.dist-info` directories that one reader tolerates and another does not.
2. **Whether each package imports**, with its `__version__`.
3. **ABSENT versus DIFFERS per pin** — a distinction the repair's single "MISSING" label erased. A package that *imports* but reports ABSENT means damaged metadata: the code is present, the `.dist-info` is not.

The pin checker itself now queries `importlib.metadata` (the mechanism `import`
uses) rather than parsing `pip list`, and refuses to report a near-empty result
as "34 missing packages" — under 20 visible distributions in an env like this
one is a failed query, and it says so instead. Repair step 1's exit code and
error lines are also surfaced now; on the first run its outcome was invisible,
so a failure there was only inferable from the verification that followed.

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
