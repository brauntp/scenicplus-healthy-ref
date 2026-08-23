# Cell-type grouping for this reference

Derived from the verified label merge (`02_pair/attach_atac_labels.py`), not
assumed. Numbers are cells carrying a `predicted_CellType_Broad` label.

## The label situation, settled

- `atac.h5ad` carries **no cell-type column** — only `Clusters` (C1–C32),
  `Sample`, and 23 QC fields.
- The transferred labels live in `atac_metadata_with_transferred_labels.tsv`,
  which has **no cell-barcode column**. The only possible join is positional.
- That join is **verified, not assumed**: all 23 columns shared between the TSV
  and `atac.h5ad.obs` (`Clusters`, `TSSEnrichment`, `FRIP`, `nFrags`, `PassQC`,
  `Sample`, `RNA_id`, `domain`, `balancing_weight`, …) agree on **all 163,969
  rows, zero mismatches**. Row *i* of the TSV is cell *i*.
- `predicted_CellType` has 48 distinct values plus **259 blank** (written as NA
  and dropped).
- `predicted_CellType_Broad` exists only on the RNA side, so ATAC's broad labels
  are derived through RNA's fine→broad correspondence (54 mappings).

## Support per group

`independent_metacell_equiv` = limiting-modality cells / k. This is the honest
denominator: with `--oversample` the pipeline produces more metacells than this,
but they share cells and are not independent observations.

| support | groups | interpretation |
|---|---|---|
| strong (>20) | 17 | eRegulons are interpretable |
| moderate (6–20) | 3 | EoBasoMast Precursor, Plasma Cell, Stromal — treat as provisional |
| weak (3–5) | 2 | Megakaryocyte Precursor (5), Late GMP (3) |
| very weak (≤2) | 2 | cDC (2), Pro-Monocyte (1) — **do not report eRegulons** |

**ATAC is the limiting modality in 21 of 24 groups**, so ATAC cell counts —
not the much larger RNA counts — set the resolution ceiling everywhere.

Three groups are limited by RNA instead: Monocyte (2,976 RNA vs 3,021 ATAC),
Plasma Cell (503 vs 561), and Stromal (490 vs 6,215). Stromal is the notable
one: 6,215 ATAC cells are wasted against 490 RNA cells, a 12.7× imbalance.

## Unlabelled cells

The label TSV leaves **259 ATAC cells** without a `predicted_CellType`, and
**3,588 RNA cells** have no `predicted_CellType_Broad`. These are excluded from
pairing and reported by name at the top of every run.

This was a bug, caught by reading the first real dry run: `.astype(str)` turns
NaN into the literal string `"nan"`, so the unlabelled cells were being pooled
into a 25th group called `nan` and given 60 metacells — a mixture of whatever
those cells are, handed to SCENIC+ as if it were a lineage. Blank strings from
the TSV arrive the same way. Both are now treated as missing.

## Recommended parameters

```
--cells-per-metacell   50
--oversample           8
--min-metacells-per-group 60
--min-cells-per-group  50      # drops nothing at broad level; Pro-Monocyte (69) survives
```

**Measured on the real data** (`--dry-run`, 2026-08-23): 25,323 metacells,
40.4 GB dense output, 14.8 GB sparse inputs, so peak RSS ~55 GB. `--mem=80G`
in `slurm/pairing.sbatch` leaves headroom. Compare `impute_accessibility` on
this reference: 241 GB, 481 GB with the transpose.

(The dry run initially reported 25,383 — the extra 60 were the spurious `nan`
group described above.)

`--oversample 16` doubles both (50,510 metacells, ~80 GB) and buys ranking but
not independence — see `docs/oversampling_tradeoff.png`.

## What NOT to do

**Do not use the 53-level `initial_CellType` for grouping.** Its rare levels
(ASDC = 7 cells, GMP-Cycle = 28, GMP-Neut = 82) cannot support metacells at any
k, and 6 of its levels have no ATAC counterpart at all.

**Do not lower k to rescue the weak groups.** Dropping to k=10 would give
Pro-Monocyte 6 independent metacells instead of 1, but k=50 was measured to give
the largest effect size on true links (`docs/pairing_sensitivity.csv`), and
lowering it globally to accommodate 69 cells degrades the 17 groups that are
well supported. If those lineages matter, the fix is more ATAC cells.
