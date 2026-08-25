# The ~96G region-to-gene figure was 2.2× too high — measured on the cluster

The probe ran on the real object in the real environment:

| | |
|---|---|
| dense equivalent (25,323 × 34,112 + 25,323 × 393,832, float32) | 40.4 GB |
| loaded subset (scRNA whole + 20,000 of 393,832 ATAC columns) | 5.10 GB |
| peak RSS including both `.to_df()` calls | 5.46 GB |
| **measured overhead factor** | **1.07×** |
| extrapolated full-object peak | **43.2 GB** |
| suggested `--mem` (peak + 25%) | **54G** |

**Two independent measurements agree.** The QC job's `sacct` MaxRSS —
`mudata.read()` plus a full NaN scan — was 43.2 GB. The probe's extrapolation,
`mudata.read()` plus both `.to_df()` calls, is also 43.2 GB. Different code
paths, measured days apart, both landing at ~1.07× the dense size. That is the
object being resident and essentially nothing else: `.to_df()` and the NaN scan
each add ~3 GB, not 40.

The mechanism, confirmed directly: `.to_df()` on a dense float32 modality
returns a frame that **shares memory** with the AnnData's X
(`np.shares_memory(...)` is True) and adds +0.00 GB. There is no second copy.
The original 96G came from assuming there was — the same class of error as the
QC overshoot, and the third time an assumed access pattern was the mistake
rather than the arithmetic.

> **Superseded for the pipeline job.** `--mem=54G` sizes *one* rule —
> `region_to_gene` — and the pipeline OOM-killed at
> `motif_enrichment_cistarget` instead, whose footprint is set by the cisTarget
> database, not by the data. `slurm/scenicplus.sbatch` now requests **128G**.
> The derivation below is still the right one for `region_to_gene` in
> isolation; see "Two rules, two terms" at the end of this document for how they
> combine.

Set `--mem=54G`, not 96G. **One caveat that keeps this a floor:**
`region_to_gene` forks joblib workers that memory-map per-worker slices from
`temp_dir`. The parent's peak is what was measured; the workers add on top. If
the SCENIC+ job OOMs, that is where the extra went — raise `--mem` and lower
`--cpus-per-task` together, since worker count scales with cores.

```bash
# cheap, seconds, login-node safe — what produced the numbers above
python 03_pipeline/probe_region_to_gene_memory.py --h5mu ACC_GEX.h5mu

# authoritative, in the SCENIC+ env once it exists
python 03_pipeline/probe_region_to_gene_memory.py --h5mu ACC_GEX.h5mu --full
```

Check `sacct` MaxRSS afterwards and add the row to the table below.

---

# Measured memory on this reference (2026-08-23)

Predictions vs `sacct` MaxRSS, so future estimates can be calibrated:

| job | requested | predicted peak | **actual MaxRSS** | verdict |
|---|---|---|---|---|
| pairing (`aggregate_atac_sparse.py`) | 80G | ~55 GB | **63.5 GB** (79%) | good estimate |
| QC, first attempt | 32G | — | OOM-killed at 27.7 GB | too small |
| QC (`qc_paired.sbatch`) | 96G | ~80.7 GB | **43.2 GB** (45%) | **1.87x over** |
| region sets, first attempt | 16G | "flat ~2 GB" | **OOM-killed (exit 137)** | **~5x under** |
| region sets (`region_sets.sbatch`) | 32G | 2.45 GB measured | **2.80 GB** (9%) | **14% over — first close one** |

**The lesson from the region-set OOM.** I sized `--mem` from the streaming
block's own size (1.89 GB at 20,000 regions × 25,323 metacells) and called the
footprint "flat", ignoring everything computed *from* the block. Measured, that
block peaked at **9.59 GB**: `scipy.stats.rankdata` returns float64 (2× a
float32 block) and argsorts int64 internally, and the tie-term computation
sorted a second full copy.

Fixed by computing ranks and tie terms from one stable argsort, keeping ranks in
the block's dtype — 2.45 GB for a 4,000-region block, with ranks identical to
`rankdata` (max |diff| 0.0, ties included) and p-values matching scipy to
5e-08. The peak/block ratio is 5.1× at 20,000 regions and 6.5× at 4,000; it
falls with block size because the per-column loop's temporaries don't scale with
it.

Four resolved predictions, and the only one within 20% is the one taken from a
direct measurement of the actual code path rather than from arithmetic about it
(15% under, 1.87x over, ~5x under, 14% over). The
pattern is not arithmetic error — it is asserting an access pattern instead of
measuring one. `03_pipeline/probe_region_to_gene_memory.py` exists for exactly
this reason, and the same discipline should apply before every remaining
`--mem`.

**The lesson from the QC job.** The 80.7 GB figure assumed
`validate_h5mu.py` holds the object resident (40.4 GB) *and* materialises a full
second copy during the NaN scan. Actual overhead above resident was 2.8 GB (7%),
so the scan works in chunks. A "one full copy" assumption is the wrong default
for a *scan*.

**It is also the wrong default for `.to_df()`** — which is what I assumed next,
and it was wrong for the same reason. Superseded: the paragraph that once stood
here argued the region-to-gene estimate should keep a 2× assumption because
`.to_df()` "genuinely does build a new DataFrame". It does not, on a dense
float32 modality: the returned frame shares memory with the AnnData's X
(`np.shares_memory(...)` is True) and adds +0.00 GB. Measured overhead 1.07×,
not 2×, giving `--mem=54G` rather than 96G. See the section at the top of this
file for the numbers.

The generalisation worth keeping: **"one full copy" is the wrong default for any
access pattern you have not measured** — a chunked scan, a block-sharing frame
constructor. Measure it (`03_pipeline/probe_region_to_gene_memory.py`) instead of
choosing an assumption.

The pairing job's 63.5 GB against a 55 GB prediction is the opposite error and
the more dangerous one: 15% under. Its dense output arithmetic was exact
(40.4 GB); the sparse-input estimate was low, as recorded in `GROUPS.md`.

---

# Memory planning: the imputed-accessibility step

> **NOT ON THE PATH THIS PROJECT TOOK — kept for reference only.**
>
> This section describes `impute_accessibility`, the step SCENIC+'s stock
> unpaired path runs before building metacells. This pipeline does **not** call
> it: at 163,969 cells × 393,832 peaks a dense `int32` matrix is ~240 GB, which
> is unrunnable here. Because mean aggregation is linear, `02_pair/
> aggregate_atac_sparse.py` aggregates the **raw sparse counts** instead, which
> was verified numerically identical to densify-then-average. The pairing job
> then ran in 63.5 GB.
>
> The `~2×` transpose budget quoted below is an **unmeasured** projection of the
> same kind that overshot the QC job by 1.87× and the region-to-gene estimate by
> 2.2×. Treat every figure in this section as an order-of-magnitude sketch, and
> if you ever do take this path, measure it rather than trusting the table.

**If you take the stock path, this is the step most likely to kill a cluster
job.**

`pycisTopic.diff_features.impute_accessibility` returns a **dense**
`regions x cells` matrix (`int32` when `scale_factor` is an integer, which the
default `10**6` is). It is not sparse, and nothing downstream re-sparsifies it.
SCENIC+'s own `process_non_multiome_data` calls it on the full object before
building metacells, and `02_pair/glue_metacells.py` consumes the result.

### Footprint of a single dense int32 matrix

| cells | 150k regions | 250k regions | 400k regions |
|---|---|---|---|
| 20,000 | 11 GB | 19 GB | 30 GB |
| 50,000 | 28 GB | 47 GB | 75 GB |
| 100,000 | 56 GB | 93 GB | 149 GB |
| 200,000 | 112 GB | 186 GB | 298 GB |

The paired object is written **cells × regions**, so a transpose may be live at
the same time. The old advice here was to budget ~2× the table above — but that
factor is **assumed, never measured**, and assumed copy factors are what
overshot the QC job by 1.87× and the region-to-gene estimate by 2.2×. If you
take this path, measure the transpose rather than doubling: a 100k-cell
reference on a 250k-peak set is 93 GB for one matrix, and whether the peak is
93 or 186 GB decides whether it fits a standard compute node.

### Mitigations, in order of preference

1. **Impute only the regions you will actually use.** `impute_accessibility`
   takes `selected_regions`. Pass the union of the binarized-topic region sets
   and the per-cell-type DARs — the regions that reach motif enrichment anyway.
   Regions outside every region set cannot contribute an eRegulon, so imputing
   them is wasted memory. This typically cuts the matrix by 3–10×.
2. **Impute per cell type, pair per cell type, concatenate.** Pairing is already
   confined within `--group-key`, so it never needs two cell types resident at
   once. Peak memory becomes the largest single cell type rather than the whole
   reference. This is the recommended route for a >100k-cell reference and needs
   a per-group driver loop rather than one whole-object call.
3. **Set `scale_factor=1`** to get `float32` instead of `int32`. Same width, so
   this saves nothing on its own — it matters only because integer scaling
   quantizes low probabilities to zero, and `10**6` is chosen to make that
   quantization useful. Do not change it expecting a memory win.
4. **Request the memory.** If the reference is small enough (≤50k cells), a
   high-memory partition is the simplest answer. Size `--mem` from a measured
   run, not from the table doubled — see the note above on assumed copy
   factors.

### Unrelated but adjacent: MALLET

`run_cgs_models_mallet` has **no memory keyword** — the signature is
`(cistopic_obj, n_topics, n_cpu, n_iter, random_state, alpha, alpha_by_topic,
eta, eta_by_topic, top_topics_coh, tmp_path, save_path, reuse_corpus,
mallet_path)`. Heap size is controlled only by the `MALLET_MEMORY` environment
variable, which must be exported **before** the call. If LDA dies with a Java
`OutOfMemoryError`, no Python argument will fix it; raise `MALLET_MEMORY`.

Also set `tmp_path` to node-local scratch: Mallet writes a serialized corpus
that is large for a 100k-cell object, and a shared filesystem makes topic
modelling I/O-bound.


## AUCell: `ENOMEM at os.fork()` is not the OOM killer

`AUCell_direct` and `AUCell_extended` both failed at 128 GB with

```
  File ".../multiprocessing/popen_fork.py", line 66, in _launch
    self.pid = os.fork()
OSError: [Errno 12] Cannot allocate memory
```

**This is a different failure from the cisTarget OOM.** There is a full Python
traceback and no `oom_kill` event: the kernel *refused* to admit another address
space, rather than killing a process that had one. The shell exit code is an
ordinary 1, so anything keying on 137 reports it as a plain error.

### Where ~118 GB goes, in order

`calculate_auc` → `score_eRegulons` → `rank_data` → `aucell4r`:

| step | scATAC | scRNA |
|---|---|---|
| `mudata.read(ACC_GEX.h5mu)` | 37.2 GB | 3.2 GB |
| `.to_df()` on both | shares memory with a dense `.X` | — |
| `rank_data()` — **both computed before either enrichment runs** | 37.2 GB | 3.2 GB |
| `aucell4r` `RawArray(c_uint32)`, region step | 37.2 GB | — |
| **peak** | | **~118 GB** |

Two details that make this larger than it looks:

- `rank_data` does `ranking = np.zeros_like(mtx)` where `mtx = df.to_numpy()`.
  The `.astype('uint32')` applies to the row *values*, which are assigned back
  into a **float32** array — so the ranking is float32, the same width as the
  source, not the uint32 the values are cast to.
- `score_eRegulons` computes `gex_ranking` **and** `acc_ranking` up front, so
  the 37.2 GB region ranking stays resident throughout the gene AUC.

Then `os.fork()` must succeed. At ~118 GB inside a 128 GB cgroup there is no
headroom to admit the child, and fork returns ENOMEM.

### MEASURED: serial AUCell succeeded (job 10721352)

`--cpus-per-task 1`, both rules, no fork:

| | value |
|---|---|
| `AUCell_extended` enrichment | 12.9 min -> `AUCell_extended.h5mu` 87 MB |
| `AUCell_direct` enrichment | 14.1 min -> `AUCell_direct.h5mu` 183 MB |
| total wall (read + enrich + write, both) | **27.8 min** |
| read of the 40 GB paired object | ~18 s per rule |

So the entire fork obstacle cost 28 minutes of single-threaded work. There is no
parallel baseline to compare against — both parallel attempts died before a
single worker started — so the speedup forgone is unknown, not small.

Each rule is a separate process, so the paired object is read once per rule
rather than shared. At ~18 s that is not worth restructuring.

### CORRECTION: `--cpus-per-task` IS the lever, `--mem` is not

The 192 GB rerun failed with a **byte-identical** traceback. That retires the
sizing argument below, which I had reasoned out rather than tested. Four
observations settle it:

1. **It failed on `gex_AUC` — the gene matrix.** That ranking is 3.2 GB; the
   37 GB region `RawArray` had not been allocated yet. The failing allocation
   was small.
2. **128 GB and 192 GB produced the same traceback**, at the same call site.
3. **The node was `cnode-08-26` on `batch`**, whose nodes carry 400 GB+.
4. Resident at that point is the paired object plus both rank matrices, ~81 GB —
   comfortably inside 192 GB.

So this was never an allocation-size failure. `fork()` duplicates page tables,
and under strict overcommit accounting the kernel charges the child's worst case
instead of trusting copy-on-write. A parent holding ~81 GB forking 16 children
implies a commit charge no `--mem` on this cluster can satisfy — which is exactly
why raising it changed nothing.

**`aucell4r` takes a serial branch at `num_workers == 1`** — verified in the
pySCENIC source: the branch contains no `RawArray`, no `Process`, and no `fork`.
That is the fix.

```bash
sbatch --cpus-per-task=1 slurm/scenicplus.sbatch \
    --target AUCell_direct.h5mu --target AUCell_extended.h5mu
```

`--cpus-per-task` sets `n_cpu`, which reaches `aucell4r` unchanged as
`num_workers` (via `signature_enrichment`). Verified against snakemake that
naming two targets builds those two and leaves `scplus_mudata` alone.

Cost: the serial branch loops `enrichment4cells` over signatures instead of
chunking them across workers. Same total work, unmeasured wall time — neither
attempt reached a single worker, so there is no parallel baseline to compare
against.

### Superseded reasoning: why I thought `--cpus-per-task` would not help

This is the **opposite** of `motif_enrichment_cistarget`, and conflating them
wastes a run:

| rule | database/array | scales with `n_cpu`? |
|---|---|---|
| `motif_enrichment_cistarget` | each joblib worker loads its own | **yes** — throttle it |
| `AUCell_*` | `RawArray` allocated **once**, then forks | **no** — peak is flat above `n_cpu=1` |

| `n_cpu` | RawArray | peak |
|---|---|---|
| 1 | none (serial branch) | ~81 GB |
| 2–16 | 37.2 GB | ~118 GB |

`n_cpu=1` does take a serial branch that skips the `RawArray` entirely, but it
still needs ~81 GB and makes AUCell single-threaded over every eRegulon. **Raise
`--mem` instead** — now 192 GB, giving ~74 GB of headroom over the peak.

### What remains after this

`eGRN_direct` and `eGRN_extended` both completed before the failure, so
`eRegulon_direct.tsv` and `eRegulons_extended.tsv` are banked. Remaining:
`AUCell_direct`, `AUCell_extended`, `scplus_mudata`. The last reads the paired
object plus both AUC objects (~41 GB) and is not the binding rule.

## MEASURED: cisTarget at `--cpus-per-task 2` (job 10719633)

The subset run succeeded. Numbers from its log, replacing the projections:

| | value |
|---|---|
| `--mem` requested | 128 GB (131072 MB) |
| `--cpus-per-task` | 2 |
| wall time | **25.5 min** for 21 region sets |
| waves | 11 (2 sets each, then 1) |
| per wave | mean 2.3 min (range 2.1–2.6) |
| `ctx_results.hdf5` | 202 MB |

**The per-worker database load is confirmed directly by the log:** `Reading
cisTarget database` appears exactly **twice per wave** — once per joblib worker —
and takes ~60 s each before any region work begins. That is the mechanism this
document describes, observed rather than inferred.

Wall-time extrapolation from the measurement, for choosing a worker count:

| `--cpus-per-task` | waves | estimated wall |
|---|---|---|
| 1 | 21 | ~51 min |
| 2 | 11 | ~25 min (measured) |
| 4 | 6 | ~13 min |

### Regions actually used: 85% of what the BEDs contain

The sets hold 5,000 regions each, but cisTarget reported using 3,020–4,812
(mean 4,257) after mapping into the database at `fraction_overlap 0.4`:

| | regions used | retention |
|---|---|---|
| best (HSC_MPP) | 4,812 | 96% |
| median | ~4,300 | 86% |
| worst (Stromal) | 3,020 | **60%** |

The peak-overlap audit predicted near-total query-side matching. At 85% mean that
is broadly right but not exact, and **Stromal loses 40%** — worth remembering
when reading its eRegulons, since a stromal region set drawn from a
haematopoietic reference has the least support in a SCREEN cCRE database built
mostly on blood and standard cell lines.

This also means the memory tables above are conservative by ~15%: the db slice is
sized on regions requested, and fewer are actually loaded.

## The real cisTarget term: the recovery-curve matrix

Equal-N region sets did **not** fix the OOM, and that ruled out my previous
explanation. The staleness check reported every output newer than the newest
`.bed`, which means the sets were already equal-N at 5,000 regions when the run
died — so a per-worker database slice of 0.6 GB cannot account for 128 GB.

The term I had missed is in `cisTarget.run_ctx()`, which calls ctxcore's
`recovery()`:

```python
rccs = np.empty(shape=(n_features, rank_threshold))     # ctxcore/recovery.py
```

`n_features` is **every motif in the database** and `rank_threshold` is
`int(ctx_rank_threshold × total_regions_in_database)`. `np.empty`'s default dtype
is float64. It is then wrapped in `df_rccs`, a DataFrame with a MultiIndex over
every curve point, so two copies are live.

At the configured `ctx_rank_threshold: 0.05` against 1,837,304 database regions:

| term | figure |
|---|---|
| curve points per motif | 0.05 × 1,837,304 = 91,865 |
| `rccs`: 32,765 × 91,865 float64 | 22.4 GB |
| `+ df_rccs` | **44.9 GB per worker** |
| db slice at 5,000-region sets, for comparison | 0.6 GB |

**None of that depends on region-set size** — which is exactly why rewriting the
sets changed nothing. And joblib runs `n_cpu` workers, so at `--cores 16` the
recovery term alone is ~718 GB.

### The lever, and what it costs

`ctx_rank_threshold` is a config parameter, so this is directly tunable:

| `ctx_rank_threshold` | points | per worker | ×8 | ×16 |
|---|---|---|---|---|
| 0.0500 | 91,865 | 44.9 GB | 359 GB | 718 GB |
| 0.0200 | 36,746 | 17.9 GB | 144 GB | 287 GB |
| 0.0100 | 18,373 | 9.0 GB | 72 GB | 144 GB |
| 0.0055 | 10,105 | 4.9 GB | 39 GB | 79 GB |

**What it does not change:** NES, AUC, and which motifs are called enriched. The
AUC integrates `rccs[:, :rank_cutoff]` where `rank_cutoff = round(auc_threshold ×
total) - 1 = 9,186`, fixed by `ctx_auc_threshold` — any `rank_threshold` at or
above the floor still covers it.

**What it does change:** `leading_edge4row` compares each motif's curve against
`avgrcc + 2·std` over the full curve length, so a leading edge beyond the
truncation point cannot be found. Such a motif would not have cleared the NES
threshold anyway, since the AUC only integrates the first 9,186 ranks.

**The floor is 0.0055, not 0.005.** scenicplus passes `int(rank_threshold ×
total)` as the curve length while ctxcore computes `round(auc_threshold × total)`
and asserts `rank_cutoff <= curve_length`. At `rank_threshold == auc_threshold`
those are 9,186 and 9,187 — the assert fails on the truncation. Recommending
equality would have named a value that crashes.

## Correction: cisTarget loads the database ONCE PER WORKER

`motif_enrichment_cistarget` was OOM-killed at `--mem=128G` while
`motif_enrichment_dem` completed in the same run. The section below assumed one
database load per rule. That is wrong, and the code says so
(`scenicplus/cli/commands.py::run_motif_enrichment_cistarget`):

```python
cistarget_results = joblib.Parallel(n_jobs=n_cpu, temp_folder=temp_dir)(
    joblib.delayed(_run_cistarget_single_region_set)(
        ..., cistarget_db_fname=cistarget_db_fname, ...)
    for key in region_set_dict)
```

`_run_cistarget_single_region_set` receives the database **filename** and
constructs its own `cisTargetDatabase`. joblib's default backend is loky —
separate processes, no shared mapping — so **`n_cpu` region sets are in flight
at once, each holding its own slice.** At `--cores 16` that is up to sixteen
concurrent loads.

Two consequences:

1. **The peak is driven by the largest few region sets, not by the total and not
   by the largest single one.** It is the sum over the `n_cpu` largest sets that
   happen to run together.
2. **`--cpus-per-task` is as much a memory lever as `--mem`**, for this rule
   specifically. Halving it roughly halves the peak.

This also explains why DEM survived: its per-worker slice is a foreground set
plus a *capped* background, not a whole cell type's DAR set.

`03_pipeline/size_cistarget_memory.py` computes both from the real `.bed` files
and the real database, rather than from an assumed geometry:

```bash
python 03_pipeline/size_cistarget_memory.py \
    --region-set-folder 01_atac/region_sets \
    --db resources/cistarget_db/hg38_screen_v10_clust.regions_vs_motifs.rankings.feather
```

### Equal-N region sets make this predictable

With the DAR sets as written the sizes span ~27× (9,796 to 264,353 regions at
`min_log2fc 0.25`), so the peak is a tail-driven sum that moves whenever a cell
type gains regions. With `choose_dar_threshold.py --top-n N --write` every set
is the same size and the peak becomes a product:

| top-n | per worker (incl. transient 2×) | ×8 | ×16 | ×21 |
|---|---|---|---|---|
| 2,000 | 0.5 GB | 4 GB | 8 GB | 10 GB |
| 5,000 | 1.2 GB | 10 GB | 20 GB | 26 GB |
| 10,000 | 2.4 GB | 20 GB | 39 GB | 51 GB |
| 20,000 | 4.9 GB | 39 GB | 78 GB | 103 GB |

At `top-n 5,000` — comfortably above the ~2,000-region recovery-curve floor
derived in `docs/REGION_SETS.md` — all 21 sets fit in a 128 GB allocation with
every worker running. That is the recommended route, and the cost is the one
already documented: membership becomes rank-based.

## Two rules, two terms — and why the job needs the larger

The 43.2 GB figure above is `region_to_gene`. It is not the binding constraint.

**The database term.** `pycistarget`'s `cisTargetDatabase.load_db()` — with
`region_sets` supplied, which the Snakefile always does — calls `db.load(subset)`
rather than `load_full()`, and `ctxcore` pushes the column list down into
`pyarrow.feather.read_table(columns=...)`. So it reads only the database regions
that *match our peaks*, not all 1,837,304 of them. That is genuinely better than
loading the whole file, but at this project's geometry the subset is still most
of it:

| term | figure |
|---|---|
| 393,832 matched regions × 32,765 motifs × int32 | 48.1 GB |
| `subset_to_pandas` builds a pyarrow Table *then* a DataFrame, both briefly live | **~96 GB peak** |
| (`load_full()` on the whole database, for reference — path not taken) | ~224 GB |

The peak-overlap audit found near-total query-side matching, so assume most
peaks match rather than hoping for a small subset.

**They do not add.** Every database-loading rule declares
`threads: config["params_general"]["n_cpu"]`, so at `--cores 16` snakemake
serialises `motif_enrichment_cistarget`, `motif_enrichment_dem`, `eGRN_direct`
and `eGRN_extended` — two are never resident at once. `--mem` must cover the
**max** of the two terms, not the sum: ~96 GB, not ~140 GB.

`slurm/scenicplus.sbatch` requests **128G**, the peak plus a third. If it still
OOMs, the lever is `--cpus-per-task` — fewer joblib workers in `region_to_gene` —
rather than `--mem` alone.

**How this was missed:** the sbatch header already said the database term "is
usually the binding constraint," and the run_pipeline triage already said
"killed with no python traceback -> almost always SLURM OOM ... the cisTarget
rankings DB is memory-resident; raise --mem." Both were correct. The `--mem`
directive contradicted them because it was set from the one rule that had been
measured, and a measured number for the wrong rule beat an unmeasured one for
the right rule.
