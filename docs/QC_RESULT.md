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

## What to distrust downstream

1. **Late GMP (3 independent metacells), Pro-Monocyte (1), Stromal (9),
   Plasma Cell (10), EoBasoMast (11).** Any eRegulon resting on these is
   supported by a handful of independent observations no matter how many
   metacells were emitted.
2. **Distal links** more than proximal ones -- the locus-specific increment over
   the lineage baseline is small.
3. **p-values from region-to-gene** as calibrated significance. Oversampling at
   8x improves ranking and inflates naive significance
   (`docs/oversample_null.csv`). Use them to rank.
