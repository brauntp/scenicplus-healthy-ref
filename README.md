# SCENIC+ on the GLUE-integrated healthy hematopoietic reference

Peak-to-gene regulatory linkage from an **unpaired** scRNA + scATAC reference
co-embedded with scGLUE, run on SLURM (ARC).

Built to run on the cluster: the analysis sandbox this was authored in cannot
reach ARC, so every script validates its own inputs, fails loudly, and logs what
it did. Nothing here assumes a successful previous step.

---

## The one design decision that matters

**We do not use SCENIC+'s built-in unpaired path.** Understanding why is the
point of this repo.

SCENIC+ handles non-paired data in `process_non_multiome_data`
(`src/scenicplus/data_wrangling/adata_cistopic_wrangling.py`), which delegates
metacell construction to `generate_pseudocells_for_numpy`
(`src/scenicplus/utils.py`). That function does this, per cell-type label:

```python
for x in range(n_pseudobulk):
    random.seed(x)
    sample_idx = sample(list(idx), n_cell)   # uniform draw within the label
    ...
```

It is called **twice, independently** — once for RNA, once for ATAC. So the RNA
cells and the ATAC cells inside metacell `HSC_7` share nothing except the label
`HSC`. Region-to-gene regression downstream (`infer_region_to_gene`, GBM +
Spearman) therefore sees only **between-cell-type** variance plus sampling
noise. Every graded within-cell-type program — the continuous HSPC → progenitor
structure that a GLUE co-embedding exists to recover — is averaged away.

Two incidental bugs in the same function, worth knowing if you ever do use it:
`random.seed(x)` is set *inside* the metacell loop on the metacell index, so
metacell 0 of every cell type draws with seed 0 (draws are correlated across cell
types); and the pipeline's `params_general.seed` never reaches this code path.

### What we do instead

`02_pair/glue_metacells.py` builds the paired object directly:

1. Anchors are chosen per cell type by **farthest-point (maximin) sampling** over
   the pooled GLUE latent coordinates, so metacells tile the manifold the cell
   type occupies instead of oversampling its dense core.
2. Around each anchor, the *k* nearest RNA cells and *k* nearest ATAC cells are
   taken **from the same point in the shared latent space** and aggregated.
3. Output is a plain 2-modality MuData with matching `obs_names` — byte-for-byte
   the schema `prepare_GEX_ACC` emits.

Because `infer_region_to_gene` reads only `mdata["scRNA"].to_df()` and
`mdata["scATAC"].to_df()`, dropping this file at the config's
`combined_GEX_ACC_mudata` path makes Snakemake treat `prepare_GEX_ACC` as
satisfied. **No fork of scenicplus is required** — one upstream script, stock
pipeline thereafter.

### Does it actually help?

Simulation with a shared latent geometry (cell-type centroids + 40 independent
within-type programs), where peak *i* truly drives gene *i*, and decoy pairs
(peak *i*, gene *j*) are genuinely unlinked:

| pairing | median ρ, true links | median ρ, decoys | AUROC |
|---|---|---|---|
| SCENIC+ stock (label-random, k=25) | 0.009 | −0.001 | 0.529 |
| GLUE-anchored (this repo, k=25) | **0.490** | 0.009 | **1.000** |

Both rows are `noise=0.15, k=25` from `docs/pairing_sensitivity.csv`
(320 metacells each); regenerate with `python docs/benchmark_pairing.py`.

Stock pairing is at chance — it recovers none of the within-type links, exactly
as the code above predicts. See `docs/pairing_benchmark.png` for how this
degrades as integration quality falls: GLUE-anchored pairing holds AUROC ≥ 0.987
while modality-specific noise stays at or below 0.6× the program scale, falls to
0.89 at noise 1.0, and reaches 0.63 at 2.5 (best k at each level). **Pairing quality is bounded by
integration quality** — which is why the diagnostics below are not optional.

---

## What this cannot tell you

Read this before interpreting any eGRN out of this pipeline.

- **Pairing is inference, not measurement.** A region-to-gene link here means
  accessibility and expression covary across the co-embedded manifold. It does
  **not** mean they were observed together in one nucleus. With true multiome
  data that distinction disappears; with GLUE-paired data it never does.
- **Errors are inherited, not detected.** If the integration mislocates a
  population, pairing confidently mixes the wrong cells and nothing downstream
  will flag it. `pairing_diagnostics.csv` reports a cross-modal gap per cell
  type for exactly this reason — a large gap means the two modalities do not
  co-occupy that region of the latent space and pairing there is extrapolation.
- **Resolution is capped by metacell size.** With *k* = 25 and ~2× cell reuse,
  a cell type contributes roughly `n/12` independent observations. Rare
  populations produce too few metacells for a GBM fit; they are dropped by
  `--min-cells-per-group` rather than silently under-powered.
- **Region sets bound what motif enrichment can find.** A TF whose sites fall
  outside the ArchR consensus peak set cannot appear in any eRegulon. Absence of
  a TF is not evidence against it.
- **A healthy reference cannot support "absent in disease" claims** in either
  direction. It is a baseline, not a contrast.

---

## Peaks: we do not re-call them

Your ArchR reproducible peak set is used as-is. pycisTopic's
`get_consensus_peaks` extends each MACS summit by `peak_half_width` and
iteratively removes less-significant overlapping peaks, citing Corces et al.
2018 — the same iterative-overlap algorithm ArchR implements, at the same
default 501 bp fixed width. Re-calling would spend cluster hours reproducing
what you have.

What this *does* force is a database decision, because region IDs must match:

- **`04_db/download_precomputed_db.sh`** — precomputed hg38 SCREEN cCRE
  database. Fast, but SCREEN regions are ~150–350 bp and will not align to
  501 bp ArchR peaks; non-overlapping peaks drop out of motif enrichment
  silently.
- **`04_db/build_cistarget_db.sh`** — custom database on *your* peak set. Every
  peak scored, no dropout. Wraps `create_cistarget_motif_databases.py` (verified
  entry point) and needs the `cbust` binary on `$PATH` and executable, an hg38
  FASTA, and many-core hours.

Precomputed database sizes are substantial: the SCREEN rankings feather is
**~33 GiB** and the scores feather **~13 GiB**. Check free disk before starting;
`download_precomputed_db.sh` resumes and verifies checksums.

**Run `04_db/peak_overlap_audit.py` before choosing.** It reports what fraction
of your peaks the precomputed database can actually represent at SCENIC+'s
`fraction_overlap_w_ctx_database` (0.4), so the trade is a measured number
rather than a guess.

---

## Pipeline

```
00_inspect/   inspect_anndata.py     RNA/ATAC schema, peak geometry, shared group keys
              inspect_archr.R        ArchRProject peaks, matrices, fragments; writes consensus BED
01_cistopic/  export_from_archr.R    PeakMatrix + metadata (+ optional per-group fragments)
              run_cistopic.py        cisTopic object, LDA (Mallet), imputed accessibility, region sets
02_pair/      glue_metacells.py      GLUE-anchored paired MuData  <-- the custom step
03_pipeline/  validate_h5mu.py       gate before a multi-hour job
              config.template.yaml   annotated SCENIC+ config
              run_pipeline.sh        validate, then snakemake
04_db/        peak_overlap_audit.py  measure DB/peak-set compatibility FIRST
              build_cistarget_db.sh  custom DB on your peaks
slurm/        *.sbatch               submit wrappers (parameterize partition/account)
```

### Two things that will bite you on the cluster

1. **`impute_accessibility` returns a DENSE `regions × cells` matrix.** At 100k
   cells × 250k peaks that is ~93 GB, and ~186 GB with the transpose live.
   Read `docs/MEMORY.md` before submitting — it lists the mitigations
   (impute only region-set regions; or impute and pair per cell type).
2. **MALLET heap is not a Python argument.** `run_cgs_models_mallet` has no
   memory keyword; heap comes only from the `MALLET_MEMORY` environment
   variable, exported before the call. A Java `OutOfMemoryError` cannot be fixed
   from Python.

### Order of operations

1. **Inspect** (minutes, login node). Run both inspectors; send back the two
   reports. They resolve: genome build, peak width, which `obs` column can serve
   as `key_to_group_by`, whether the GLUE latent is where we expect, and which
   cell types are too small or present in only one modality.
2. **Audit the database choice** (`peak_overlap_audit.py`) and decide precomputed
   vs. custom on the measured dropout.
3. **Export + cisTopic/LDA** (the long ATAC step; Mallet needs `MALLET_MEMORY`).
4. **Pair** (`glue_metacells.py`) → **validate** (`validate_h5mu.py`).
5. **Run SCENIC+** on the paired object.

### Parameters worth thinking about, not defaulting

| parameter | default | why you might change it |
|---|---|---|
| `--cells-per-metacell` | 25 (script default; **use 50**) | k=50 has the highest median ρ on true links at every noise level in `docs/pairing_sensitivity.csv`, so larger is better for signal. The cost is observations, not robustness: k=10/25/50 → 798/320/160 metacells for the GBM to regress on. Drop below 50 only when a cell type is too small to yield a stable fit. |
| `--cells-per-metacell-atac` | = RNA | ATAC is sparser; a larger ATAC *k* is often justified. |
| `--min-cells-per-group` | 50 | Below this a cell type cannot support a metacell; raise it rather than trust a thin population. |
| `--anchor-mode` | `both` | `rna` or `atac` anchors on one modality's geometry when the other is much sparser. |
| `search_space_upstream/downstream` | 1 kb–150 kb | SCENIC+'s default. Anchor to insulation boundaries rather than a fixed span if you have Hi-C for these cell types. |

---

## Provenance

Written against **scenicplus v1.0a2** (commit `840dab8`) and the pycisTopic
`main` branch, with every function signature and config field verified against
the cloned source rather than the tutorials. `docs/pairing_sensitivity.csv` and
`docs/pairing_benchmark.png` are regenerated by `python docs/benchmark_pairing.py`
(add `--quick` for a 3-point sweep); it imports `02_pair/glue_metacells.py`
directly, so the benchmarked code is the shipped code.
