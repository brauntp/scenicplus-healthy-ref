#!/usr/bin/env python3
"""
Per-cell-type region-to-gene correlation. The score SCENIC+ does not compute.

WHY THIS EXISTS
---------------
`scenicplus.enhancer_to_gene.calculate_regions_to_genes_relationships` scores
every (region, gene) pair ONCE over all metacells:
`_score_regions_to_single_gene(acc, exp, ...)` takes plain arrays and has no
grouping parameter, and `cli/commands.py:infer_region_to_gene` calls it a single
time on `mdata["scATAC"].to_df()` / `mdata["scRNA"].to_df()`. Grep either module
for `group_by` / `cell_type` / `per_group`: there is nothing. So the shipped
`rho` in `region_to_gene_adj.tsv` answers "do this peak and this gene covary
across the whole reference", which is dominated by cell-type identity and cannot
say WHICH cell type a link belongs to.

This script recomputes ONLY the Spearman correlation, per cell-type group, over
the pairs the pipeline already nominated.

WHAT IT REUSES vs RECOMPUTES
----------------------------
REUSED (cell-type-independent, do not recompute):
  * the candidate pair set -- `search_space.tsv`, or equivalently the pairs
    present in `region_to_gene_adj.tsv`. It is +/-150 kb TAD-capped genomic
    distance, so it does not depend on cell type.
  * `importance` -- the GBM feature score. Deliberately NOT refit per group:
    that is the expensive half (it was the pipeline's whole region_to_gene
    rule), and it is a feature-selection score rather than a directional claim,
    so it is more defensible as a cell-type-agnostic prior. The output carries
    the global value under `importance_global` so nothing looks recomputed that
    was not.

RECOMPUTED per group:
  * `rho` -- Spearman, that group's metacells only.

WHAT THIS DOES NOT DO
---------------------
It does not make a link "cell-type specific" on its own. A high within-group rho
in one group is not evidence of specificity unless it is also LOW elsewhere; the
output therefore reports every group's rho side by side plus a specificity
contrast, and leaves thresholding to the caller. Groups too small to detect a
correlation are reported with their detection floor rather than a misleading
near-zero rho -- see `--min-metacells` and the `min_detectable_rho` column.

MEMORY
------
The dense accessibility matrix is 393,832 x 25,323 float32 = 37 GB. It is read
in blocks of `--block` regions from a BACKED object, exactly as
01_cistopic/region_sets_from_metacells.py does, so peak RSS is set by the block
rather than the object. Ranks are computed once per block with the same
`_ranks_and_ties` primitive (one stable argsort, ranks kept in the input dtype);
using scipy.rankdata here instead peaked at 9.59 GB on a 20,000-region block and
was OOM-killed, which is why that helper exists.
"""
import argparse
import pathlib
import sys

try:
    import numpy as np
except ImportError:                                              # pragma: no cover
    sys.exit("ERROR: this python has no 'numpy'.\n"
             f"       interpreter: {sys.executable}\n"
             "       Activate the pairing env (conda activate scplus-pairing).")

import pandas as pd


# ---------------------------------------------------------------- rank helper
def _ranks(X):
    """Average ranks per COLUMN, same dtype as X, from one stable argsort.

    Lifted from 01_cistopic/region_sets_from_metacells.py:_ranks_and_ties minus
    the tie-term accumulator (Spearman needs the ranks; the tie correction there
    was for the Mann-Whitney variance). Keeping ranks in the input dtype is
    load-bearing: scipy.stats.rankdata returns float64 and argsorts int64
    internally, which doubled peak RSS and OOM-killed a 16G job.
    """
    order = np.argsort(X, axis=0, kind="stable")
    Xs = np.take_along_axis(X, order, axis=0)
    neq = np.empty(Xs.shape, dtype=bool)
    neq[0] = True
    np.not_equal(Xs[1:], Xs[:-1], out=neq[1:])
    del Xs
    grp = np.cumsum(neq, axis=0, dtype=np.int32)
    del neq
    G = grp.max(0)
    ranks_sorted = np.empty_like(X)
    for j in range(X.shape[1]):
        cnt = np.bincount(grp[:, j], minlength=int(G[j]) + 1)[1:]
        ends = np.cumsum(cnt)
        starts = ends - cnt
        ranks_sorted[:, j] = np.repeat(
            ((starts + ends + 1) / 2.0).astype(X.dtype), cnt)
    del grp
    R = np.empty_like(X)
    np.put_along_axis(R, order, ranks_sorted, axis=0)
    return R


def _spearman_cols(Ra, rg):
    """Spearman between each column of rank-matrix Ra and rank-vector rg.

    Both inputs are ALREADY ranks, so this is Pearson on ranks -- the definition
    of Spearman. float64 accumulators: ranks are float32 to halve memory, and
    summing thousands of values near 2.5e4 loses precision otherwise.
    """
    A = Ra.astype(np.float64, copy=False) - Ra.mean(0, dtype=np.float64)
    g = rg.astype(np.float64, copy=False) - rg.mean(dtype=np.float64)
    num = A.T @ g
    den = np.sqrt((A ** 2).sum(0) * (g ** 2).sum())
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(den > 0, num / den, np.nan)


def min_detectable_rho(n, alpha=0.05, power=0.80):
    """Smallest |rho| detectable at n observations, via Fisher z."""
    from scipy.stats import norm
    if n <= 3:
        return np.nan
    z = norm.ppf(1 - alpha / 2) + norm.ppf(power)
    return float(np.tanh(z / np.sqrt(n - 3)))


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--h5mu", type=pathlib.Path, required=True,
                    help="paired object (ACC_GEX.h5mu). Read BACKED.")
    ap.add_argument("--adj", type=pathlib.Path, required=True,
                    help="region_to_gene_adj.tsv -- supplies the candidate pairs "
                         "AND the global rho/importance carried through")
    ap.add_argument("--group-key", default="predicted_CellType_Broad")
    ap.add_argument("--out", type=pathlib.Path, required=True,
                    help="output table (.tsv, or .csv.gz if the name says so)")
    ap.add_argument("--block", type=int, default=4_000,
                    help="regions per streaming block (default 4000; peak RSS is "
                         "~5-6.5x the block's float32 size, measured 2.5 GB at 4000)")
    ap.add_argument("--min-metacells", type=int, default=30,
                    help="skip groups with fewer metacells than this (default 30). "
                         "Below ~30 the detection floor exceeds |rho|=0.5 and a "
                         "near-zero rho would be indistinguishable from no power.")
    ap.add_argument("--oversample", type=int, default=8,
                    help="metacell oversampling factor used during pairing "
                         "(default 8). Only used to report an independence-"
                         "corrected detection floor -- it does not change any rho.")
    ap.add_argument("--dry-run", action="store_true",
                    help="plan only: report groups, blocks, memory. Writes nothing.")
    args = ap.parse_args()

    for p in (args.h5mu, args.adj):
        if not p.exists():
            sys.exit(f"ERROR: not found: {p}")

    import mudata

    print("=" * 74)
    print("per-cell-type region-to-gene correlation")
    print("=" * 74)

    adj = pd.read_table(args.adj)
    need = {"target", "region", "importance", "rho"}
    if not need <= set(adj.columns):
        sys.exit(f"ERROR: {args.adj} lacks {sorted(need - set(adj.columns))}\n"
                 f"       columns present: {list(adj.columns)}")
    print(f"  candidate pairs : {len(adj):,} over {adj.target.nunique():,} genes "
          f"and {adj.region.nunique():,} regions")

    md = mudata.read(str(args.h5mu), backed=True)
    acc, gex = md["scATAC"], md["scRNA"]
    if args.group_key not in gex.obs.columns:
        sys.exit(f"ERROR: --group-key {args.group_key!r} not in scRNA.obs\n"
                 f"       available: {list(gex.obs.columns)}")
    # Align by NAME, never position: the modalities are written independently and
    # a positional assumption would scramble every correlation while still
    # producing plausible numbers.
    if not (acc.obs_names == gex.obs_names).all():
        sys.exit("ERROR: scATAC and scRNA obs_names differ in order or content; "
                 "this script assumes the paired object's shared index.")
    groups = gex.obs[args.group_key].astype(str)
    n_mc = len(groups)
    sizes = groups.value_counts()
    keep = sizes[sizes >= args.min_metacells]
    skip = sizes[sizes < args.min_metacells]
    print(f"  metacells       : {n_mc:,}   groups: {len(sizes)}")
    print(f"  groups kept     : {len(keep)}  (>= {args.min_metacells} metacells)")
    for g, k in skip.items():
        print(f"    skip {g!r}: {k} metacells "
              f"(floor |rho| >= {min_detectable_rho(k):.2f})")

    # Restrict to pairs whose region and gene are both present in the object.
    reg_idx = pd.Index(acc.var_names)
    gen_idx = pd.Index(gex.var_names)
    ok = adj.region.isin(reg_idx) & adj.target.isin(gen_idx)
    if not ok.all():
        print(f"  dropping {(~ok).sum():,} pairs whose region or gene is absent "
              f"from the object")
    adj = adj[ok].copy()
    if adj.empty:
        sys.exit("ERROR: no candidate pair survives the region/gene intersection. "
                 "Wrong object, or region ids in a different format?")

    regions = pd.Index(adj.region.unique())
    n_reg = len(regions)
    n_blocks = int(np.ceil(n_reg / args.block))
    blk_gb = args.block * n_mc * 4 / 1024 ** 3
    print(f"  regions touched : {n_reg:,} in {n_blocks} blocks of {args.block:,}")
    print(f"  block resident  : {blk_gb:.2f} GB "
          f"(full matrix {n_reg * n_mc * 4 / 1024 ** 3:.0f} GB); "
          f"peak ~{blk_gb * 6.5:.1f} GB")

    if args.dry_run:
        print()
        print("  --dry-run: nothing written.")
        print(f"  recommended: --mem {max(8, int(blk_gb * 6.5 * 1.3) + 4)}G")
        print("=" * 74)
        return

    # Row positions per group, and the gene expression ranks (small: genes are
    # 34k not 394k, so the whole gene block per group fits comfortably).
    pos = {g: np.flatnonzero((groups == g).to_numpy()) for g in keep.index}
    reg_pos = pd.Series(np.arange(len(reg_idx)), index=reg_idx)
    gen_pos = pd.Series(np.arange(len(gen_idx)), index=gen_idx)

    genes_needed = pd.Index(adj.target.unique())
    print(f"  reading expression for {len(genes_needed):,} genes ...", flush=True)
    Xg = np.asarray(gex[:, genes_needed].X, dtype=np.float32)
    gcol = pd.Series(np.arange(len(genes_needed)), index=genes_needed)

    # Per-group gene ranks, computed once (not per block).
    grank = {g: _ranks(Xg[pos[g]]) for g in keep.index}
    del Xg

    # Pair bookkeeping: for each region block, which (pair-row, gene-col) to hit.
    adj = adj.reset_index(drop=True)
    adj["_rblk"] = (regions.get_indexer(adj.region) // args.block)
    out = {g: np.full(len(adj), np.nan, dtype=np.float32) for g in keep.index}

    for b in range(n_blocks):
        sel = adj.index[adj._rblk == b]
        if not len(sel):
            continue
        blk_regions = regions[b * args.block:(b + 1) * args.block]
        Xa = np.asarray(acc[:, blk_regions].X, dtype=np.float32)
        bcol = pd.Series(np.arange(len(blk_regions)), index=blk_regions)
        rcol = bcol.reindex(adj.loc[sel, "region"]).to_numpy()
        gc = gcol.reindex(adj.loc[sel, "target"]).to_numpy()
        for g in keep.index:
            Ra = _ranks(Xa[pos[g]])
            Rg = grank[g]
            # Correlate each needed (region, gene) pair. Group by gene so the
            # column slice is contiguous per gene rather than per pair.
            for gi in np.unique(gc):
                m = gc == gi
                r = _spearman_cols(Ra[:, rcol[m]], Rg[:, gi])
                out[g][sel[m]] = r
            del Ra
        del Xa
        if b % 5 == 0:
            print(f"    block {b + 1}/{n_blocks} "
                  f"({(b + 1) * args.block:,}/{n_reg:,} regions)", flush=True)

    res = adj[["region", "target", "Distance"]].copy() if "Distance" in adj.columns \
        else adj[["region", "target"]].copy()
    res["importance_global"] = adj["importance"].to_numpy()
    res["rho_global"] = adj["rho"].to_numpy()
    for g in keep.index:
        res[f"rho__{g}"] = out[g]

    rho_cols = [c for c in res.columns if c.startswith("rho__")]
    M = res[rho_cols].to_numpy(dtype=np.float64)

    # Specificity must NOT be computed on raw |rho|. Raw correlations are not
    # comparable across groups of different size: measured on a pure-noise
    # fixture, the 95th percentile of |rho| was 0.054 at n=1200 but 0.202 at
    # n=90 -- nearly 4x. A best-minus-second contrast on raw |rho| therefore
    # structurally favours the SMALLEST group, which would have manufactured
    # "specificity" for exactly the groups with the least evidence.
    #
    # Fisher's z divided by its standard error is a standard normal under H0
    # regardless of n. On the same fixture all three group sizes gave median
    # |z| 0.65-0.72 and 95th percentile 1.84-1.91 (expected 0.67 / 1.96), i.e.
    # one shared scale. The contrast is computed on THAT.
    nvec = np.array([int(keep[c.removeprefix("rho__")]) for c in rho_cols],
                    dtype=np.float64)
    Z = np.arctanh(np.clip(M, -0.999999, 0.999999)) * np.sqrt(
        np.maximum(nvec - 3.0, 1.0))
    for c, z in zip(rho_cols, Z.T):
        res[c.replace("rho__", "z__")] = z

    absZ = np.abs(Z)
    filled = np.where(np.isnan(absZ), -np.inf, absZ)
    srt = np.sort(filled, axis=1)
    best, second = srt[:, -1], srt[:, -2]
    bi = np.argmax(filled, axis=1)
    res["best_group"] = [rho_cols[i].removeprefix("rho__") for i in bi]
    res["best_abs_z"] = np.where(np.isfinite(best), best, np.nan)
    res["second_abs_z"] = np.where(np.isfinite(second), second, np.nan)
    # Standardised contrast: how much more evidence the best group carries than
    # the runner-up, on a scale where 1.96 is the usual two-sided 5% point.
    res["specificity_z"] = res.best_abs_z - res.second_abs_z
    # The raw rho of the winning group, for interpretation (effect size), kept
    # separate from the ranking statistic so neither is mistaken for the other.
    res["best_rho"] = M[np.arange(len(M)), bi]
    res["n_groups_scored"] = np.isfinite(M).sum(1)

    if str(args.out).endswith((".csv.gz", ".csv")):
        res.to_csv(args.out, index=False)
    else:
        res.to_csv(args.out, sep="\t", index=False)

    floors = pd.DataFrame({
        "group": keep.index,
        "metacells": keep.to_numpy(),
        "min_detectable_rho": [min_detectable_rho(k) for k in keep],
        "min_detectable_rho_indep": [
            min_detectable_rho(max(k // args.oversample, 4)) for k in keep]})
    fo = pathlib.Path(str(args.out).rsplit(".", 1)[0]
                      .replace(".csv", "") + ".detection_floors.csv")
    floors.to_csv(fo, index=False)

    print()
    print(f"  wrote {args.out}  ({len(res):,} pairs x {len(rho_cols)} groups)")
    print(f"  wrote {fo}")
    print()
    print("  READ THE FLOORS FILE. A near-zero rho in a small group means "
          "'undetectable',")
    print("  not 'absent': at 60 metacells the floor is |rho| ~ 0.36 raw and "
          "~0.89 after")
    print("  the oversampling correction.")
    print()
    print("  RANK ON `specificity_z`, NOT on raw rho differences. Raw |rho| is "
          "not")
    print("  comparable across groups -- under pure noise its 95th percentile "
          "was 0.054")
    print("  at n=1200 but 0.202 at n=90, so a raw contrast favours the "
          "smallest group.")
    print("  `specificity_z` is best-minus-second on Fisher z/SE, which shares "
          "one null")
    print("  scale across sizes. `best_rho` is the effect size for "
          "interpretation.")
    print("  A uniformly-coupled (housekeeping) link scores near zero by "
          "design.")
    print("=" * 74)


if __name__ == "__main__":
    main()
