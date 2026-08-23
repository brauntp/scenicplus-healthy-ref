# The ~96G region-to-gene figure is probably too high

That number came from assuming `.to_df()` builds a full copy on top of the
40.4 GB resident object. **Measured, it does not.** On a dense float32 `.h5mu`,
`.to_df()` returns a frame that shares memory with the AnnData's X
(`np.shares_memory(...)` is True) and adds +0.00 GB resident.

| assumption | implied peak | implied `--mem` |
|---|---|---|
| `.to_df()` copies (the old projection) | 80.8 GB | ~101G |
| measured overhead factor 1.34× | 54.1 GB | ~68G |

This is the third time an assumed access pattern was the error rather than the
arithmetic, so the figure is **not** being lowered on this observation alone.
Two things could invalidate it:

1. The no-copy behaviour depends on the pandas/anndata versions in the SCENIC+
   environment — a different env from the one measured. A version that
   consolidates blocks would copy.
2. `region_to_gene` forks joblib workers that memory-map per-worker slices from
   `temp_dir`. The parent's peak is a **floor**, not a ceiling.

So measure it there instead of projecting again:

```bash
# cheap, seconds, from a login node
python 03_pipeline/probe_region_to_gene_memory.py --h5mu ACC_GEX.h5mu

# authoritative: in the SCENIC+ env, in a job sized generously
python 03_pipeline/probe_region_to_gene_memory.py --h5mu ACC_GEX.h5mu --full
```

It reports peak RSS and a suggested `--mem`. Set the SCENIC+ job from that,
then check `sacct` MaxRSS against it afterwards and add the row to the table
below.

---

# Measured memory on this reference (2026-08-23)

Predictions vs `sacct` MaxRSS, so future estimates can be calibrated:

| job | requested | predicted peak | **actual MaxRSS** | verdict |
|---|---|---|---|---|
| pairing (`aggregate_atac_sparse.py`) | 80G | ~55 GB | **63.5 GB** (79%) | good estimate |
| QC, first attempt | 32G | — | OOM-killed at 27.7 GB | too small |
| QC (`qc_paired.sbatch`) | 96G | ~80.7 GB | **43.2 GB** (45%) | **1.87x over** |

**The lesson from the QC job.** The 80.7 GB figure assumed
`validate_h5mu.py` holds the object resident (40.4 GB) *and* materialises a full
second copy during the NaN scan. Actual overhead above resident was 2.8 GB (7%),
so the scan works in chunks. A "one full copy" assumption is the wrong default
for a *scan*.

It is not the wrong default for `.to_df()`, which genuinely does build a new
DataFrame — so the SCENIC+ region-to-gene estimate below keeps the 2x
assumption, but with this measurement as a reason to check MaxRSS on the first
run rather than trusting the projection.

The pairing job's 63.5 GB against a 55 GB prediction is the opposite error and
the more dangerous one: 15% under. Its dense output arithmetic was exact
(40.4 GB); the sparse-input estimate was low, as recorded in `GROUPS.md`.

---

# Memory planning: the imputed-accessibility step

**This is the step most likely to kill a cluster job.** Read before submitting.

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

The paired object is written **cells × regions**, so a transpose is live at the
same time: **budget ~2× the table above** for peak RSS. A 100k-cell reference on
a 250k-peak set needs roughly 190 GB — more than a standard compute node.

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
   high-memory partition is the simplest answer. Set `--mem` from the table
   above, doubled.

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
