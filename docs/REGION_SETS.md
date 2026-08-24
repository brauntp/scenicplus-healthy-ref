# Region sets: the cap was choosing the regions, not the thresholds

## What the first real run produced

`region_sets_from_metacells.py` on `ACC_GEX.h5mu`, at the defaults
`--min-log2fc 0.25 --fdr 0.05 --max-regions 20000`:

- 21 BEDs written (24 groups, minus Pro-Monocyte, cDC and Late GMP on
  independent-observation count — the diagnostics path was used, confirmed by
  `independence_source` in the summary)
- **20 of the 21 groups capped at exactly 20,000 regions.** Only Stromal came in
  under, at 9,796.
- **409,796 BED lines total — 1.04× the entire 393,832-peak set.**

## Why that is a defect and not a result

When almost every group caps, the cap is doing the selecting. Each set is an
arbitrary top-N slice by effect size, ranked among far more regions that also
passed, and 20,000 is a number with no biological meaning. Two groups whose true
biology differs in breadth come out identical in size.

**The FDR filter cannot bind at this scale.** With 25,323 metacells:

| shift in P(superiority) | z | p |
|---|---|---|
| +0.02 | 2.1 | 1.6e-02 |
| +0.05 | 5.4 | 4.0e-08 |
| +0.10 | 10.7 | 3.5e-27 |
| +0.20 | 21.5 | 1.4e-102 |

A 2% shift is already significant, and oversampling at 8× makes the p-values
anticonservative on top of that (`docs/oversample_null.csv`: nominal 0.05 inflates
to 0.338). So `--min-log2fc` is the only filter doing work — and 0.25 was a guess
made before any effect-size distribution existed.

That a group is "differentially accessible" at 5% of all peaks is also not a
useful claim. It is most of the accessible genome.

## The fix: choose the threshold from the measured distribution

Two additions, and neither needs a second pass over the matrix:

1. `region_sets_from_metacells.py --dump-stats <file.npz>` writes per-group
   effect sizes and q-values (two 21 × 393,832 float32 arrays, 63 MB uncompressed, less on disk since savez_compressed is used). It also now warns
   explicitly when most groups cap, rather than printing `(capped)` and moving
   on.
2. `01_cistopic/choose_dar_threshold.py` reads that dump and sweeps thresholds
   offline — regions per group, how many still cap, how much of the peak set the
   union covers — then rewrites the BEDs at the chosen value.

```bash
# one more generator run, this time keeping the statistics
sbatch slurm/region_sets.sbatch          # with DUMP_STATS set

# then, offline and instant:
python 01_cistopic/choose_dar_threshold.py --stats 01_atac/dar_stats.npz
python 01_cistopic/choose_dar_threshold.py --stats 01_atac/dar_stats.npz \
    --min-log2fc <chosen> --out-dir 01_atac/region_sets --write
```

The rewrite deletes the previous threshold's BEDs first: SCENIC+ globs the
directory, so a leftover file is a silent extra enrichment run at the old
setting.

## Validated on a fixture built to reproduce the defect

200 strongly-shifted regions per group plus 4,000 weakly-shifted ones, so
`log2fc 0.25` admits far more than the cap:

| min_log2fc | median/group | capped | total BED | union % of peaks |
|---|---|---|---|---|
| 0.25 | 1,961 | **6 of 6** | 9,000 | 78.4% |
| 0.50 | 341 | 0 | 2,041 | 24.5% |
| 1.00 | 200 | 0 | 1,200 | 14.1% |
| 3.00 | 174 | 0 | 1,050 | 13.1% |

The plateau at 200 is exactly the planted strong signal, and the rewrite at
`log2fc 1.0` recovered **1,200 of 1,200 planted regions with 0 false
positives** — the weak-shift regions correctly excluded. The tool finds the
lowest threshold that neither caps nor empties a group, and says which
direction each error goes.

## The sweep ran, and no shared threshold works

Real sweep on `01_atac/dar_stats.npz` (21 groups, 393,832 regions, FDR 0.05):

| min_log2fc | min/group | median | max | capped |
|---|---|---|---|---|
| 0.25 | 9,796 | 75,548 | 264,353 | 20 |
| 1.00 | 8,156 | 45,881 | 145,252 | 17 |
| 1.50 | 2,811 | 24,540 | 64,028 | 14 |
| 2.00 | **473** | 13,247 | 34,341 | 5 |
| 3.00 | **7** | 1,818 | 11,700 | 0 |

**The cap "clears" at 3.00 only because the thin groups collapse.** At that
threshold the smallest group has 7 regions. The spread between smallest and
largest never closes — 27× at 0.25, 1,671× at 3.00 — because the groups' effect
size distributions genuinely differ that much.

### The missing constraint: how small is too small

`ctx_auc_threshold` is 0.005, so cisTarget integrates its recovery curve over
the top 0.5% of the 1,837,304-region ranking — the top 9,186 ranks. A region set
of size N puts N × 0.005 regions in that window under the null:

| set size | expected in AUC window | usable? |
|---|---|---|
| 7 | 0.04 | no — statistic undefined in practice |
| 473 | 2.4 | no — NES 3.0 is noise at this count |
| 2,811 | 14.1 | yes |
| 20,000 | 100.0 | yes |

So ~2,000 regions is where the recovery curve stabilises. The sweep's original
`zero` column only counted sets at exactly 0, which meant a 7-region set was
reported as healthy. It now counts sets below `--min-usable` (default 2,000) and
labels the column **starved**.

Reading the real table against that floor: 2.00 fails both ways (473 is
unusable, 5 groups still cap) and 3.00 is worse. There is no threshold where
nothing caps and every group stays usable.

### The alternative: equal-N per group

`--top-n N` takes the top N FDR-passing regions per group by effect size,
ignoring `--min-log2fc` entirely. Every group gets the same number, so no group
starves and none caps.

```bash
python 01_cistopic/choose_dar_threshold.py --stats 01_atac/dar_stats.npz \
    --top-n 5000 --out-dir 01_atac/region_sets --write
```

**Stromal's 9,796 is the binding constraint** — it is the smallest FDR-passing
count at 0.25, so any `--top-n` up to that is feasible for all 21 groups.

**What equal-N costs, stated plainly.** Set membership becomes rank-based rather
than effect-based: a region admitted in one group would be rejected in another,
because the implied log2fc cutoff differs per group. The tool prints that span
so the cost is visible. In exchange, every set is the same size, which makes NES
values comparable across groups — under a shared threshold they are not, since a
20,000-region set and a 473-region set have very different null distributions.

## What this does not settle

The right threshold on the real data is an empirical question the sweep answers,
not one this document can. Two things to weigh when reading it:

- **cisTarget runtime scales with total BED lines**, each region scored against
  1,837,304 database regions. 409,796 lines is the expensive end.
- **The thinner groups starve first.** Stromal (9 independent metacells) was the
  only group under the cap at 0.25; it will be the first to empty as the
  threshold rises. A group that empties contributes no eRegulon at all, which is
  a stronger statement than a small region set.
