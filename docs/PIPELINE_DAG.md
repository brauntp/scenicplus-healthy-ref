# The SCENIC+ DAG, for an analyst who wants to verify it

**The Snakefile is not in this repository.** It ships inside the installed
`scenicplus` package at `<site-packages>/scenicplus/snakemake/Snakefile`, and
`03_pipeline/run_pipeline.sh:116-129` resolves it at runtime via
`importlib.resources` — the same lookup upstream's own `scenicplus
init_snakemake` uses.

That is a real interpretability problem: the computational core of this
pipeline is invisible at review time unless you first build a ~900-package
conda env. This document is the map, so the DAG can be checked without
running anything.

Get the file itself:

```bash
bash docs/fetch_snakefile.sh      # fetches the pinned tag, verifies SHA256
```

Everything below was derived from **scenicplus v1.0a2**,
`sha256 b0f2e1553838068bb10d151c93e2d950deb2dc72681bc727644fd2fe6a21517a`,
437 lines. Line numbers refer to that file. If `fetch_snakefile.sh` reports a
checksum mismatch, treat this document as stale until re-verified.

---

## Where this repo's custom step enters

This is the single most important thing to understand about the run.

```
  01_cistopic/run_cistopic.py ─┐
                               ├─► 02_pair/glue_metacells.py ──► paired MuData
  GLUE latent (RNA + ATAC)  ───┘         (THIS REPO)                   │
                                                                       │
  ════════════════════ everything below is STOCK scenicplus ═══════════╪════
                                                                       ▼
                                          config output_data.combined_GEX_ACC_mudata
```

The Snakefile's own `prepare_GEX_ACC` rule (line 11 multiome / line 28
non-multiome) **is never executed on this project**. `glue_metacells.py` writes
its output directly to the path the config names as
`output_data.combined_GEX_ACC_mudata`, so snakemake sees the rule's output
already present and treats it as satisfied.

That substitution only holds because of `--rerun-triggers mtime`
(`run_pipeline.sh:395-400`). Under snakemake's default triggers
(`mtime, params, input, code, software-env`) the `input`/`code` triggers would
fire — the hand-built file has no recorded provenance matching the rule — and
snakemake would decide the MuData is out of date and **overwrite it with the
stock label-random pairing**, silently discarding the entire reason this repo
exists. See README "The one design decision that matters".

**Verify it yourself before a long run:**

```bash
bash 03_pipeline/run_pipeline.sh --config <your.yaml> --validate-only
bash 03_pipeline/run_pipeline.sh --config <your.yaml> -n   # dry run: prepare_GEX_ACC must NOT appear
```

If `prepare_GEX_ACC_multiome` or `prepare_GEX_ACC_non_multiome` shows up in the
dry-run plan, stop — the run is about to destroy the paired object.

---

## The rule graph

14 nodes: `all` plus 13 compute rules. Two rules are conditional and resolve to
exactly one form at parse time, so 13 always run.

```
download_genome_annotations (221)
   │  genome_annotation, chromsizes
   ├──────────────────────────────┐
   │                              ▼
   │                     motif_enrichment_dem (93 | 147)
   │                              │  dem_result_fname, dem_html
   ▼                              │
get_search_space (237) ◄── paired MuData
   │  search_space                │
   │                              │
   │     motif_enrichment_cistarget (50)
   │              │  ctx_result_fname, ctx_html
   │              ▼               ▼
   │          prepare_menr (197) ◄┘   ◄── paired MuData
   │              │  tf_names, cistromes_direct, cistromes_extended
   │              │
   ├──► region_to_gene (283) ◄── paired MuData
   │        │  region_to_gene_adjacencies
   │        │
   └──► tf_to_gene (260) ◄── paired MuData, tf_names
            │  tf_to_gene_adjacencies
            ▼
    ┌───────┴────────┐
    ▼                ▼
eGRN_direct (306)   eGRN_extended (347)      both also read input_data.ctx_db_fname
    │                │
    ▼                ▼
AUCell_direct (388) AUCell_extended (404)    both also read the paired MuData
    │                │
    └────────┬───────┘
             ▼
      scplus_mudata (420)  ──►  output_data.scplus_mdata   ← rule all (5)
```

### Per-rule reference

| rule (line) | CLI invoked | key inputs | outputs |
|---|---|---|---|
| `prepare_GEX_ACC_multiome` (11) | `prepare_data prepare_GEX_ACC` | cisTopic obj, GEX anndata | `combined_GEX_ACC_mudata` |
| `prepare_GEX_ACC_non_multiome` (28) | `prepare_data prepare_GEX_ACC` | cisTopic obj, GEX anndata | `combined_GEX_ACC_mudata` |
| `motif_enrichment_cistarget` (50) | `grn_inference motif_enrichment_cistarget` | `region_set_folder`, **`ctx_db_fname`**, motif annotations | `ctx_result_fname`, ctx HTML |
| `motif_enrichment_dem` (93 \| 147) | `grn_inference motif_enrichment_dem` | `region_set_folder`, **`dem_db_fname`**, motif annotations | `dem_result_fname`, dem HTML |
| `prepare_menr` (197) | `prepare_data prepare_menr` | dem + ctx results, paired MuData | `tf_names`, `cistromes_direct`, `cistromes_extended` |
| `download_genome_annotations` (221) | `prepare_data download_genome_annotations` | — (network) | `genome_annotation`, `chromsizes` |
| `get_search_space` (237) | `prepare_data search_spance` [^typo] | paired MuData, annotation, chromsizes | `search_space` |
| `tf_to_gene` (260) | `grn_inference TF_to_gene` | paired MuData, `tf_names` | `tf_to_gene_adjacencies` |
| `region_to_gene` (283) | `grn_inference region_to_gene` | paired MuData, `search_space` | `region_to_gene_adjacencies` |
| `eGRN_direct` (306) | `grn_inference eGRN` | tf2g, r2g, `cistromes_direct`, **`ctx_db_fname`** | `eRegulons_direct` |
| `eGRN_extended` (347) | `grn_inference eGRN` | tf2g, r2g, `cistromes_extended`, **`ctx_db_fname`** | `eRegulons_extended` |
| `AUCell_direct` (388) | `grn_inference AUCell` | `eRegulons_direct`, paired MuData | `AUCell_direct` |
| `AUCell_extended` (404) | `grn_inference AUCell` | `eRegulons_extended`, paired MuData | `AUCell_extended` |
| `scplus_mudata` (420) | `grn_inference create_scplus_mudata` | both AUCell + both eRegulon tables, paired MuData | `scplus_mdata` |

[^typo]: `search_spance` is spelled that way in the Snakefile (line 250). It is a
genuine upstream typo in the registered subcommand name, not a transcription
error here — it must match the CLI, so do not "fix" it.

---

## Four properties worth checking before you trust a run

Each is a claim this repo relies on. Each is stated here with how to re-derive
it, so none has to be taken on faith. `docs/verify_claims.py` checks all four
automatically once `fetch_snakefile.sh` has run.

**1. The cisTarget database is read by three rules, not one.**
`input_data.ctx_db_fname` appears at lines 53 (`motif_enrichment_cistarget`),
311 (`eGRN_direct`) and 352 (`eGRN_extended`). This is why database size drives
memory in three separate places, not just motif enrichment —
`slurm/scenicplus.sbatch:48`.

```bash
grep -n 'ctx_db_fname' docs/reference/Snakefile.v1.0a2
```

**2. DEM is not optional.** The `dem_balance_number_of_promoters` conditional
(line 92/146) swaps *which* DEM rule is defined, but **both branches declare
`dem_result_fname` and `output_fname_dem_html` as unconditional `output:`
entries**, and `prepare_menr` (line 197) takes `dem_result_fname` as an
`input:`. So DEM always runs and always gates the rest of the DAG. This is why
`04_db/download_precomputed_db.sh` refuses `--skip-scores`: the run would die
at the DEM rule *after* cisTarget had already spent its hours.

**3. No rule declares a `conda:` directive** — zero, across all 437 lines.
Every rule shells out to a bare `scenicplus` on `$PATH`, so snakemake does
**not** manage the environment. The env must be active in the shell before
snakemake starts, which is why `slurm/scenicplus.sbatch` activates it rather
than relying on a profile (`environment.yml:285-289`).

```bash
grep -c 'conda:' docs/reference/Snakefile.v1.0a2    # must print 0
```

**4. The config surface is exactly 67 keys over 151 lookups.**
`run_pipeline.sh:220` enumerates all 67 and fails fast on a missing one, so a
stale config is caught in seconds rather than inside snakemake hours later.
The enumeration has been diffed against the Snakefile and matches exactly —
no missing keys, no phantom keys.

---

## Reading the run while it is going

`run_pipeline.sh` prints the resolved Snakefile path and the full snakemake
command line before invoking (line 721), so the log records which file actually
ran. After a failure, the surviving outputs are the fastest orientation:

```bash
python 03_pipeline/pipeline_status.py --config <your.yaml> --workdir <dir>
```

Snakemake resumes from completed outputs, so a multi-hour partial run is not
lost. Logs are under `<workdir>/.snakemake/log/`.
