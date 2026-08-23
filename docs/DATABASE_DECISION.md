# cisTarget database: precomputed SCREEN, or custom build?

## The audit result (real peaks, 2026-08-23)

`peak_overlap_audit.py`, 393,832 consensus peaks against the 1,837,304-region
hg38 SCREEN cCRE catalogue, `fraction_overlap > 0.4` reciprocal:

| | |
|---|---|
| representable in DB | **319,802 (81.20%)** |
| silently dropped | **74,030 (18.80%)** |
| — zero overlap with any cCRE | 61,072 (**82% of the loss**) |
| — overlaps but below threshold | 12,958 (18% of the loss) |
| unique DB regions recruited | 319,726 |
| peak → DB region collapse | 1.00 |
| median best-overlap fraction | 1.000 |

The tool called this BORDERLINE, which is correct: 18.8% is too large to ignore
and too small to settle the question by itself.

## What the numbers already rule out

**A threshold or width artifact.** 82% of the loss is peaks with *zero* overlap
with any cCRE — not near-misses. Loosening `fraction_overlap` recovers at most
the other 18%. The width mismatch is not the problem either. SCREEN regions are
150–350 bp (median 272) against 500 bp peaks, and pycistarget passes a pair if
*either* overlap fraction clears 0.4, so the two routes demand different
absolute overlaps:

| route | overlap required |
|---|---|
| query side (overlap / DB width) | 0.4 × 150–350 = **60–140 bp**, median 109 |
| target side (overlap / peak width) | 0.4 × 500 = **200 bp**, flat |

Because every SCREEN region is narrower than 500 bp, the query-side requirement
is the smaller of the two for **100%** of regions — so the effective bar is a
median 109 bp of overlap, not 200. Narrow database regions make representability
*easier*, not harder. Confirmed empirically: the median best-overlap fraction
among retained peaks is 1.000, i.e. the typical retained peak swallows its
matched region whole.

(Only 20.5% of SCREEN regions are themselves ≤200 bp; that figure is about the
width distribution and is not what the threshold compares against.)

**A localised annotation gap.** Per-chromosome loss is nearly uniform: worst is
chr4 at 22.5%, only 3.7 points above the 18.80% mean. For comparison, a fixture
with 20% of peaks placed at random showed an 11.2-point spread. No chromosome is
disproportionately unrepresented.

**A collapse problem.** 1.00 peaks per recruited DB region means retained peaks
map essentially one-to-one; the database is not merging distinct peaks into
shared regions.

So the loss is real, unbiased across the genome, and consists of accessible
regions ENCODE never nominated as candidate cCREs.

## What the numbers cannot settle

Whether those 61,072 unannotated peaks matter. SCREEN is built from a fixed
panel of biosamples; elements specific to cell types outside that panel are
absent by construction, and bone-marrow progenitor subsets are exactly the
populations that can fall in that gap. "Not in SCREEN" is therefore not evidence
of "not real."

**An 18.8% loss is cheap if the dropped peaks are weak, ubiquitous, or
technical. It is expensive if they are the cell-type-restricted distal elements
an eGRN analysis exists to find.** Nothing in the audit distinguishes these.

## The test that decides it

`04_db/characterize_dropped.py` asks the data rather than the annotation, on
three axes:

1. **prevalence** — fraction of cells in which the peak is accessible
2. **strength** — mean accessibility
3. **specificity** — max-group mean / mean over groups

Decision rule, encoded in the script:

- dropped peaks **more** cell-type-specific (rank-biserial > +0.1) →
  **custom database justified**; the precomputed one is discarding the
  informative tail
- dropped peaks **less prevalent and no more specific** →
  **precomputed is fine**; the loss is weak or technical peaks
- neither → **genuinely borderline**; decide on cost and report the dropped
  fraction

Validated on two fixtures with known ground truth: when dropped peaks were
constructed as real restricted enhancers it returned "custom justified"
(specificity effect +0.298, p = 7e-39); when constructed as sparse noise it
returned "precomputed is fine" (prevalence effect −1.000, specificity −0.333).

Streams the sparse matrix in 2,000-cell row blocks. The only dense allocation is
the per-group sum, 24 × 393,832 × 8 B = 0.07 GB, so peak memory stays under
~1.5 GB at any cell count — login-node safe.

```bash
source setenv.sh
python 04_db/characterize_dropped.py \
    --atac      "$REF/atac.h5ad" \
    --per-peak  audit_screen.per_peak.csv \
    --obs-tsv   "$LABELS" \
    --group-key predicted_CellType_Broad \
    --out       dropped_characterization
```

## The characterisation ran — and its verdict was wrong

`characterize_dropped.py` on the real object (163,969 cells, labels from the
sidecar TSV, all rows matching):

| metric | retained | dropped | effect | p |
|---|---|---|---|---|
| prevalence | 0.0030 | 0.0014 | **−0.488** | ~0 |
| strength | 0.0053 | 0.0025 | **−0.497** | ~0 |
| specificity | 4.1626 | 4.8003 | **+0.195** | ~0 |

It returned "CUSTOM DATABASE JUSTIFIED" on the specificity line. **Do not act on
that.** The decision rule I encoded fired on a confounded metric.

`specificity = max-group mean / mean over groups` is not prevalence-invariant. A
peak accessible in 230 cells spread over 24 groups of very unequal size will
concentrate somewhere by chance. Simulating peaks with **zero** true cell-type
preference at this reference's real group sizes (69 to 31,618 cells):

| prevalence | apparent specificity under the null |
|---|---|
| 0.0030 (retained median) | 2.549 |
| 0.0014 (dropped median) | 3.268 |

The null alone predicts dropped peaks scoring **1.282×** higher. The observed
ratio is **1.153×** — *below* the null expectation. And the prevalence effect
(−0.488) is 2.5× larger than the specificity effect and points the other way:
dropped peaks are mostly just sparser.

## The corrected test

`04_db/reanalyze_dropped.py` re-analyses the per-peak CSV that
`characterize_dropped.py` already wrote — no second pass over the matrix — and
compares specificity **at matched prevalence**.

Quantile stratification alone is insufficient, which is worth recording because
it looked sufficient: on a fixture whose only difference was prevalence, ten
strata still left a pooled effect of +0.047 that cleared the permutation null —
a false positive. Even inside a decile the dropped peaks sit at the low end and
the retained at the high end. So the primary test is 1:1 nearest-neighbour
matching on prevalence (2% relative tolerance, without replacement), which
drives the matched prevalence ratio to 1.000.

Validated on two fixtures identical but for one thing:

| fixture | naive effect | matched effect | verdict |
|---|---|---|---|
| sparser only, no real specificity difference | +0.345 | **+0.009** (p=0.24) | precomputed adequate |
| sparser *and* genuinely more specific | +0.862 | **+0.673** (p≈0) | custom justified |

```bash
python 04_db/reanalyze_dropped.py \
    --per-peak dropped_characterization.per_peak.csv \
    --out      dropped_stratified
```

**Prediction, stated before the run so it can be wrong:** since the observed
specificity ratio (1.153) is below what prevalence alone predicts (1.282), I
expect the matched effect to be near zero or negative — "precomputed is
adequate". That is inference from summary statistics; the matched test on the
per-peak data decides.

## Cost of each branch, for when the verdict lands

**Precomputed SCREEN.** One 32.8 GiB download, no compute. Loses 18.8% of peaks
in motif enrichment, and that loss must be stated in the methods. Note the
database is memory-resident during `pycistarget` — the 32.8 GiB feather is read
whole.

**Custom build on the consensus peak set.** `create_cistarget_motif_databases.py`
scores every one of the 393,832 peaks, so nothing is dropped. Needs the `cbust`
binary, an hg38 FASTA, the v10 motif collection, and many-core hours; the
resulting database is also smaller than SCREEN's, because 393,832 regions is 21%
of 1,837,304 (4.67× fewer) — which reduces the memory footprint of every
downstream rule that holds it resident.

That last point cuts against the intuition that the custom route is the
expensive one: it costs CPU once and saves memory on every subsequent run.
