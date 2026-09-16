# Review map — where to spend attention

This codebase is ~19,000 lines of code. **The science is about 300 of them.**

Nothing in the repo previously said which 300, so a reviewer reads all of it at
uniform attention and drowns. This file is the triage: three tiers, what each
one is for, and what to actually ask of each file.

| tier | what it is | files | review |
|---|---|---|---|
| **1 — Science** | choices that change the result | 8 | **read closely** |
| **2 — Correctness** | plumbing where a silent bug corrupts results | 6 | spot-check the named risk |
| **3 — Support** | read-only, diagnostic, or recovery | 18 | skip unless you use it |

**How to review this repo:** read tier 1 and answer the sign-off questions at
the bottom. Tier 2 only needs the one named risk checked per file. Tier 3 cannot
change a number in a results table — it inspects, sizes, reports or recovers.

A note on what you are reviewing: this code was AI-generated. Reading it line by
line to look for "AI mistakes" is the wrong frame and will not scale to 19k
lines. What matters is whether the **method** is right and whether the code
**does what it claims** — which is why tier 1 is scoped to decisions, and why
each entry names the question to ask rather than the lines to admire.

---

## Tier 1 — Science

Each of these makes a choice that moves the results. Read the cited function,
then answer the question.

### `02_pair/glue_metacells.py` — the pairing algorithm
**The core of the whole project.** Everything downstream inherits it.

- **Read:** `build_metacells_for_group`, lines 75–158. That is ~50 lines of real
  algorithm: farthest-point (maximin) anchor selection over the pooled latent
  coordinates, then k nearest neighbours pulled from *each* modality around each
  anchor. Then `aggregate` (158–166) — a plain column mean.
- **Ask:** is maximin the right anchor rule? It tiles the manifold rather than
  sampling its density, so metacells over-represent sparse regions *by design*.
  Is that what you want, or should anchors follow cell density?
- **Ask:** neighbours are taken by Euclidean distance in GLUE latent space with
  no distance cap. A cell type whose modalities are poorly aligned still gets
  k neighbours — just far ones. `median_crossmodal_gap` (line ~140) reports
  this but nothing enforces it. Should there be a hard cutoff?
- **⚠ Default vs recommendation:** `--cells-per-metacell` defaults to **25**,
  but README and PROTOCOL both say **use 50**. Running the bare command gives
  you the un-recommended value. Decide which is right and make them agree.

### `docs/benchmark_pairing.py` — the evidence for the central claim
- **Read:** all 225 lines. It is the only support for "GLUE-anchored pairing
  beats stock" (AUROC 1.000 vs 0.529).
- **Ask:** the simulation builds cell-type centroids plus 40 independent
  within-type programs, with peak *i* truly driving gene *i*. Is that generative
  model fair to the stock method, or does it construct exactly the structure
  maximin sampling is good at? This is the single most contestable claim in the
  repo.
- Good sign: it imports `build_metacells_for_group` from the shipped file
  (line 39), so it benchmarks the real code, not a copy.

### `02_pair/qc_paired.py` — the positive control on REAL data
The benchmark above is a simulation. This is the empirical evidence, and it
tests something the diagnostics structurally cannot.

- **Read:** the docstring (lines 1–40) and `main` (90). The logic: for each
  marker gene, correlate its expression against mean accessibility of peaks
  within `--window` of its TSS, across metacells — scored against a null built
  from **distant peaks matched on accessibility**.
- **Why it matters:** `median_crossmodal_gap` is a property of the *embedding*.
  It can look excellent while the object is useless — if labels were
  misassigned, or the right-distance cells were the wrong cells, the gap would
  not notice. This check would.
- **Ask:** the reported result is 30/32 markers against 4.5 expected by chance.
  Is the null fair? It matches on accessibility but the distant peaks are still
  drawn from the same metacells, so any global covariance structure (e.g. depth
  or cell-count effects per metacell) inflates both arms. Confirm that is
  handled, because this number is the strongest real-data claim in the repo.
- **Ask:** marker TSS coordinates come from `00_inspect/marker_tss_hg38.tsv`.
  Confirm the build matches the peaks (GRCh38 throughout).

### `01_cistopic/region_sets_from_metacells.py` — what motif enrichment can see
Region sets bound every eRegulon. A TF with no site in these sets cannot appear.

- **Read:** `_mwu_from_ranks` (161) and the selection logic in `main` (272).
- **Ask:** DARs are called per cell-type label, so the sets are label-driven by
  construction — shared or continuous programs are found poorly. Is a
  label-driven region set the right input for a method whose selling point is
  recovering *within*-label structure?
- **Ask:** `--min-independent` defaults to 5, and `--assume-oversample` adjusts
  the effective N. Metacells share cells (~2× reuse), so they are not
  independent observations. Is the correction right?
- Runs `--self-test`.

### `01_cistopic/choose_dar_threshold.py` — the effect-size cutoff
- **Read:** `refine` (53). It searches for a `min_log2fc` that lands the region
  count between `--min-usable` (2,000) and `--max-regions` (20,000).
- **Ask:** the threshold is chosen to hit a *count target*, not a biological
  effect size. That is defensible for a fixed compute budget, but it means the
  cutoff is set by how many regions you wanted, not by what counts as a real
  difference. Is that acceptable, and is it stated in the methods?

### `04_db/peak_overlap_audit.py` — which peaks survive into the database
- **Read:** lines 20–45 only. The rest is I/O and indexing.
- **The rule:** a peak is kept if `Overlap_query > frac` **OR**
  `Overlap_target > frac`, where query is normalised by the *database* region
  width and target by *your peak* width. So a 501 bp peak containing a 272 bp
  cCRE scores 1.0 and is kept.
- **Ask:** this mirrors pycistarget's own filter, so it is the right rule to
  audit with. Confirm it still matches the installed pycistarget version — if
  upstream changes the filter, this audit silently measures the wrong thing.
- Runs `--self-test`.

### `05_report/celltype_rho.py` — a statistic SCENIC+ does not compute
This invents a number. It needs the most scrutiny of anything in tier 1.

- **Read:** `_spearman_cols` (101) and `min_detectable_rho` (116).
- **Ask:** `min_detectable_rho(n, alpha=0.05, power=0.80)` bakes in a power
  calculation. Are those the right α and power, and is the formula right?
- **Ask:** `--oversample` defaults to **8**, and the independent-N floor uses
  `k // 8`. That is an assumption that 8 metacells carry 1 independent
  observation. Where does 8 come from, and does it match the actual cell reuse
  reported in `pairing_diagnostics.csv`? **This one number scales every
  significance claim in the table.**

### `05_report/make_usable_subset.py` — which cell types get published
- **Read:** `main` (74), especially lines 104–107.
- **⚠ Inconsistent threshold:** `--min-metacells` defaults to **400** here but
  **30** in `celltype_rho.py` — a >13× difference applied to the same table.
  Plausibly intentional (a detection floor vs a publication-grade subset), but
  **an analyst must confirm and document that**, because a reader of the two
  scripts cannot tell which is the operative inclusion rule.

---

## Tier 2 — Correctness-critical plumbing

Not scientific choices, but a silent bug here corrupts results. Check only the
named risk.

| file | the risk to check |
|---|---|
| `02_pair/attach_atac_labels.py` | A **positional** merge of labels onto cells. If row order ever differs, every cell is silently mislabelled and nothing downstream notices. Confirm the verification it claims actually fires. |
| `02_pair/aggregate_atac_sparse.py` | Aggregates sparse ATAC without densifying. Confirm it produces the same values as the dense path in `glue_metacells.aggregate`. |
| `03_pipeline/validate_h5mu.py` | The gate before a 12–24 h job. Confirm it would actually *reject* a malformed object, not just pass a well-formed one. |
| `03_pipeline/build_genome_annotation.py` | TSS positions set the region-to-gene search space. A wrong build or coordinate convention shifts every link. Confirm the genome build matches the peaks. |
| `05_report/export_peak_gene_links.py` | `--min-abs-rho` (0.05) and `--top-n-per-gene` (20) filter the **main deliverable**. Confirm the published table states them. |
| `04_db/fetch_db_regions.py` | Reads region IDs from a remote feather footer via range request. Confirm the IDs match the database actually downloaded. |

---

## Tier 3 — Support code (skip)

Cannot change a number in a results table. Read only if you are using the tool
or debugging it.

**Inspectors / structural readers (read-only):** `00_inspect/inspect_anndata.py`,
`inspect_h5ad_lite.py`, `compare_labels.py`, `03_pipeline/peek_h5mu.py`,
`04_db/peaks_to_bed.py`

**Sizing / status / diagnostics:** `03_pipeline/pipeline_status.py`,
`probe_region_to_gene_memory.py`, `size_cistarget_memory.py`, `_check_pins.py`

**Post-hoc characterisation (reports on dropped peaks, changes nothing):**
`04_db/characterize_dropped.py`, `reanalyze_dropped.py`

**Reporting / export (renders results computed elsewhere):**
`05_report/summarize_eregulons.py`, `extract_for_plots.py`,
`export_accessibility_tracks.py`, `export_browser_tracks.py`

**Meta:** `docs/verify_claims.py` (checks the docs against their sources),
`docs/benchmark_oversample.py` (a supporting sweep, not the main claim)

**Note:** `01_cistopic/run_cistopic.py` (1,067 lines) is mostly tier 3 — it
wraps pycisTopic's canonical path with memory management. The scientific
content is the topic-count choice and the region-set export, which are covered
by the two `01_cistopic` files in tier 1.

---

## Sign-off questions

The review is done when these are answered. None requires reading the
implementation — they are method questions the code merely instantiates.

1. Is maximin anchor selection the right way to place metacells, given it
   deliberately over-samples sparse regions of the manifold?
2. Should poorly-aligned cell types be **excluded automatically**, rather than
   reported in `pairing_diagnostics.csv` and left to the reader?
3. Is the benchmark's generative model a fair test, or does it favour the
   proposed method by construction?
4. Is `k=25` or `k=50` the operative default? (They currently disagree.)
5. Where does `--oversample 8` come from, and does it match observed cell reuse?
6. Which `--min-metacells` is the inclusion rule — 30 or 400?
7. Are label-driven DAR region sets appropriate for a method that claims to
   recover within-label structure?
8. Is the positive control's null (distant peaks, accessibility-matched, same
   metacells) free of shared structure that would inflate both arms?

Questions 4, 5 and 6 are inconsistencies found while building this map. They
are not necessarily errors, but each needs a decision recorded.

---

## What this map does not give you

It scopes the reading. It does not prove the code does what it claims — nothing
here is executable. There are **no tests** in this repository; two scripts carry
a `--self-test` flag (`peak_overlap_audit.py`, `region_sets_from_metacells.py`)
and that is the whole of it.

The next step, if this review needs to conclude anything stronger than "the
method is sound as described," is a test suite: the analyst states in plain
English what must be true, and those become runnable assertions. Tests written
by whoever wrote the code — human or model — encode the same misunderstandings,
so the assertions have to come from the reviewer.
