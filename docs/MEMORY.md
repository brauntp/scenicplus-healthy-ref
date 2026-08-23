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
