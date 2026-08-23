# QC result for the paired object (`ACC_GEX.h5mu`, 2026-08-23)

**Verdict: usable.** Locus linkage is 9.4x chance (p = 5e-23). The tool's first
report said WARN, which was a threshold bug, not a property of the data --
explained below because the same mistake is easy to repeat when reading these
numbers.

## Structure: exactly as predicted

25,323 metacells in both modalities, identical `obs_names`, 24 groups, no junk
group, no all-zero metacells. Matches the `--dry-run` projection.

## Locus linkage: the pairing works

Expression vs mean accessibility within ±150 kb of each marker's TSS, across
metacells, against accessibility-matched peaks on other chromosomes.

| quantity | value |
|---|---|
| median rho, true locus | 0.2188 |
| median rho, matched decoys | 0.1802 |
| markers beating their own null p95 | **31/66 = 47%** |
| enrichment over chance | **9.4x** (5% expected) |
| binomial p vs chance | **5.3e-23** |

### Why 47% is a pass, not a warning

The first run called this WARN against a hardcoded 50% threshold. That threshold
was calibrated on a simulation whose decoy median rho was ~0. Here it is 0.1802,
and that difference is the whole story:

metacells are grouped **by cell type**, and peaks and genes share lineage
programs. A peak that is open in lymphoid cells correlates with a gene expressed
in lymphoid cells whatever chromosome it sits on. So the accessibility-matched
distal null is not a null for *linkage* -- it is a null for *cell-type structure
exists*, which it does.

Verified by simulation (no locus-specific linkage at all, peaks and genes driven
by shared lineage programs): decoy median rho reached 0.79 while the per-gene
"beats its own null p95" rate stayed at 8.8%. So shared structure inflates the
absolute correlation of true and decoy loci **alike**, and leaves the per-gene
test calibrated. The enrichment over chance transfers between datasets; the raw
correlation does not.

`qc_paired.py` now reports the binomial test instead of comparing to a fixed
fraction. 47% at 9.4x chance passes; the row-shuffled control fails at 1.8x
(p = 0.11).

**Practical reading:** links are real, but the locus-specific increment over the
lineage baseline is modest (0.2188 vs 0.1802 on the median). Proximal links
deserve more confidence than distal ones, and a link's support should be judged
against `independent_metacell_equiv`, not `n_metacells`.

## Marker specificity: 76% correct, and the misses are informative

Of 74 marker-group pairs, 76% place the marker in a group it is listed for
(strict rank-1: 68%). Fifteen miss, and they are not random:

**11 of 15 land one maturation stage away within the same lineage.** MPO, ELANE,
PRTN3 and CEBPA all mark "Early GMP" in the table but peak in Late GMP; CD79A
marks Pre-B and peaks in B; GATA1 marks MEP and peaks in Megakaryocyte
Precursor. These are the marker table being more precise than granulopoietic or
B-cell maturation actually is -- not evidence of mislabelled cells. Reference
markers separate lineages well and adjacent stages of one lineage poorly.

**4 of 15 cross lineages, and every one involves a thinly supported group:**

| marker | listed | observed peak | thin group (independent metacells) |
|---|---|---|---|
| SATB1 | LMPP | Naive T | -- both strong; SATB1 is genuinely T-lineage |
| TFRC | Early Erythroid | EoBasoMast Precursor | EoBasoMast, 11 |
| AHSP | Late Erythroid | Stromal | Stromal, 9 |
| XBP1 | Plasma Cell | Late GMP | Plasma Cell 10, Late GMP 3 |

SATB1 is the one genuine marker-table error: it is a T-cell factor, and finding
it highest in Naive T is correct biology. The other three route through groups
with 3-11 independent metacells, where a group mean is unstable regardless of
labelling.

**Conclusion on labels:** no evidence of a systematic labelling problem. The
lineage-level assignment is sound; adjacent-stage boundaries within a lineage are
soft, which is expected from label transfer and matters when interpreting an
eRegulon attributed to Early vs Late GMP.

## Cross-modal gaps: an erythroid-specific finding

From `pairing_diagnostics.csv`, worst 8 of 24 groups:

| group | gap | RNA | ATAC | independent metacells |
|---|---|---|---|---|
| Early Erythroid | **0.3888** | 10,832 | 2,001 | 40 |
| Late Erythroid | **0.3812** | 9,759 | 4,565 | 91 |
| Late GMP | 0.1867 | 853 | 153 | 3 |
| Early Lymphoid | 0.1357 | 21,656 | 4,353 | 87 |
| Early GMP | 0.1217 | 29,846 | 4,376 | 87 |
| Pre-B | 0.1120 | 29,660 | 9,523 | 190 |
| Megakaryocyte Precursor | 0.1107 | 433 | 293 | 5 |
| Stromal | 0.1073 | 490 | 6,215 | 9 |

**This is not a sampling artifact.** Gap vs independent metacells: rho = +0.13
(none). Gap vs RNA/ATAC imbalance: rho = +0.55, but the two worst groups are
well supported (40 and 91 independent metacells) and Stromal — the most
imbalanced group in the reference at 0.08 RNA/ATAC — sits 8th.

The two erythroid groups are 1st and 2nd of 24 at ~3.5x the typical group.

**Scale.** Gaps are distances between L2-normalised 50-dim GLUE vectors, where
unrelated cells sit ~1.41 apart. So 0.39 is 27% of the random-pair distance:
RNA and ATAC erythroid cells are still in broadly the same region, but far more
separated than any other lineage.

**A plausible mechanism, not a verified one.** Late erythroid maturation
involves global chromatin condensation and transcriptional shutdown before
enucleation, so the RNA and ATAC manifolds may genuinely diverge along that
trajectory — and GLUE aligns modalities through a peak-to-gene prior, which has
less to anchor on when accessibility and transcription decouple. This is
consistent with the observation but not established by it; ruling in requires
looking at the GLUE embedding itself, not the diagnostics.

**Convergent evidence is absent, and that is worth stating.** Erythroid/Mk
markers do miss more often in section 2 (4/14 = 29% vs 11/60 = 18% elsewhere),
but Fisher exact gives OR 1.78, p = 0.30 — not significant at n=14. The gap
finding rests on the diagnostics alone.

**Note on the prediction.** Before seeing this, the expectation was that Stromal,
Pro-Monocyte and cDC would be worst, reasoning from modality imbalance and low
cell counts. That was wrong: Stromal is 8th and Pro-Monocyte is not in the top 8.
Cross-modal gap and cell abundance turned out to be nearly uncorrelated, so the
two must be read as independent axes of caution — a group can be abundant and
badly aligned (Early Erythroid), or thin and well aligned (Pro-Monocyte).

## What to distrust downstream

1. **Thin support** — Late GMP (3 independent metacells), Pro-Monocyte (1),
   Stromal (9), Plasma Cell (10), EoBasoMast (11). Any eRegulon resting on these
   is supported by a handful of independent observations no matter how many
   metacells were emitted.
2. **Poor cross-modal alignment — a SEPARATE axis.** Early Erythroid (gap 0.389)
   and Late Erythroid (0.381) are abundant and well sampled but the worst-aligned
   groups in the reference. Their metacells pair RNA and ATAC cells that sit
   further apart in the GLUE space than any other lineage, so erythroid
   region-to-gene links carry more pairing error than their metacell counts
   suggest. Late GMP is bad on both axes and should be treated as the least
   reliable group overall.
3. **Distal links** more than proximal ones -- the locus-specific increment over
   the lineage baseline is small.
4. **p-values from region-to-gene** as calibrated significance. Oversampling at
   8x improves ranking and inflates naive significance
   (`docs/oversample_null.csv`). Use them to rank.
