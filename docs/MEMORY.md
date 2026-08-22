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
