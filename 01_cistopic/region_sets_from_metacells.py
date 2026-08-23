#!/usr/bin/env python3
"""
Region sets for motif enrichment, computed from the PAIRED object -- no
pycisTopic, no MALLET, no LDA.

WHY THIS EXISTS
---------------
`input_data.region_set_folder` is a required config entry: SCENIC+ globs
`region_set_folder/<family>/*.bed` and runs one cisTarget/DEM enrichment per
BED. Nothing in this pipeline had produced it, and the canonical source is
`01_cistopic/run_cistopic.py` -- which at this reference's scale means a
~40 GB Matrix Market export of 393,832 x 163,969, then MALLET LDA (hours to
days, Java heap controllable only through MALLET_MEMORY, serialized corpus on
scratch).

The config's own comments name `DARs_cell_type` as a valid region-set family
alongside `topics_otsu`. Per-cell-type differentially accessible regions are
computable directly from `ACC_GEX.h5mu` -- 25,323 metacells already grouped by
`predicted_CellType_Broad` -- in minutes, streaming, without the ATAC matrix
ever being fully resident.

WHAT YOU GIVE UP
----------------
This is not equivalent to topic modelling, and the difference matters:

  * Topics are unsupervised co-accessibility programs. They can span cell
    types, split one label into sub-programs, or capture a gradient that no
    label names. DARs are label-driven by construction, so any regulatory
    program not aligned to `--group-key` has no region set of its own.
  * eRegulons can only be discovered for programs that some region set
    represents. A DAR-only run therefore finds cell-type-associated eRegulons
    well and shared/continuous ones poorly.
  * Region sets are also the runtime driver (one enrichment per BED), so a DAR
    run with 24 groups is much cheaper than a topic run with 100+ topics --
    part of why it finds less.

For a healthy hematopoietic reference where the question is "which TFs drive
which cell type", DARs are the aligned choice. If the question becomes
"what shared programs exist across the hierarchy", add topics later: the two
families coexist as separate subdirectories and the DAG runs both.

METHOD
------
For each group in `--group-key`, a one-vs-rest comparison on the metacell
accessibility matrix:

  * effect size  = mean(in-group) - mean(rest), on the matrix as stored
  * significance = Mann-Whitney U, per region, BH-corrected within the group
  * a region enters the group's BED if it clears BOTH `--min-log2fc` and
    `--fdr`, capped at `--max-regions` by effect size

Mann-Whitney rather than a t-test because metacell accessibility is bounded
and right-skewed, not normal. Ranks are computed ONCE per block and every
group's U is derived from its rank sum -- scipy's mannwhitneyu re-ranks on each
call, which at 24 groups meant ranking the same data 24 times (4.8 h projected,
against a 4 h walltime). The closed form includes the tie and continuity
corrections and agrees with scipy to 5e-08; accessibility is tie-heavy, and
without the tie correction the disagreement reaches 7e-03.

Streaming: the matrix is read in blocks of `--block` regions from a backed
object, so peak RSS is set by the block, not by the object. At 393,832 peaks
x 25,323 metacells the whole matrix would be 37 GB; a 20,000-region block is
1.9 GB.

Usage
-----
    python 01_cistopic/region_sets_from_metacells.py \\
        --h5mu ACC_GEX.h5mu \\
        --group-key predicted_CellType_Broad \\
        --out-dir 01_atac/region_sets

    # then point the config at the PARENT directory:
    #   region_set_folder: <repo>/01_atac/region_sets
    # which will contain  DARs_cell_type/<group>.bed
"""
import argparse
import json
import re
import sys
from pathlib import Path

try:
    import numpy as np
except ImportError:                                              # pragma: no cover
    sys.exit("ERROR: needs numpy")

_COORD = re.compile(r"^(chr[0-9A-Za-z_]+)[:\-](\d+)[\-_](\d+)$")


def parse_regions(names):
    """var_names -> (chrom, start, end) arrays. Accepts chr:start-end and
    chr-start-end. Anything unparseable is reported, not silently dropped."""
    chrom, start, end, bad = [], [], [], []
    for i, n in enumerate(names):
        m = _COORD.match(str(n))
        if not m:
            bad.append((i, n))
            continue
        chrom.append(m.group(1))
        start.append(int(m.group(2)))
        end.append(int(m.group(3)))
    return (np.array(chrom), np.array(start, dtype=np.int64),
            np.array(end, dtype=np.int64), bad)


def bh(p):
    """Benjamini-Hochberg, returning q in the input order."""
    p = np.asarray(p, dtype=float)
    n = p.size
    order = np.argsort(p)
    q = np.empty(n, dtype=float)
    ranked = p[order] * n / np.arange(1, n + 1)
    # enforce monotonicity from the largest p downwards
    q[order] = np.minimum.accumulate(ranked[::-1])[::-1]
    return np.clip(q, 0, 1)


def _tie_terms(X):
    """sum(t**3 - t) over tie groups, per column.

    Needed for the Mann-Whitney variance: without it, tied values inflate p.
    Metacell accessibility is tie-heavy (many exact zeros), so this is not a
    negligible correction -- omitting it disagreed with scipy by up to 7e-03.
    """
    Xs = np.sort(X, axis=0)
    out = np.zeros(X.shape[1], dtype=np.float64)
    n = X.shape[0]
    for j in range(X.shape[1]):
        col = Xs[:, j]
        idx = np.flatnonzero(np.diff(col)) + 1
        sizes = np.diff(np.concatenate(([0], idx, [n])))
        out[j] = float(((sizes.astype(np.float64) ** 3) - sizes).sum())
    return out


def _mwu_from_ranks(Rk, mask, tie, norm):
    """One-sided (greater) Mann-Whitney p from precomputed per-column ranks.

    Equivalent to scipy.stats.mannwhitneyu(a, b, axis=0,
    alternative="greater") with the normal approximation, tie correction and
    continuity correction -- verified to 5e-08 both with and without ties.
    Derived from ranks so the expensive sort happens once per block instead of
    once per group.
    """
    n1 = int(mask.sum())
    n2 = Rk.shape[0] - n1
    N = n1 + n2
    U = Rk[mask].sum(0) - n1 * (n1 + 1) / 2.0
    mu = n1 * n2 / 2.0
    var = n1 * n2 / 12.0 * ((N + 1) - tie / (N * (N - 1)))
    sd = np.sqrt(np.maximum(var, 1e-12))
    return norm.sf((U - mu - 0.5) / sd)


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--h5mu", type=Path, required=True)
    ap.add_argument("--modality", default="scATAC")
    ap.add_argument("--group-key", required=True)
    ap.add_argument("--out-dir", type=Path, required=True,
                    help="PARENT region-set directory; a 'DARs_cell_type' "
                         "subdirectory is created inside it, because SCENIC+ "
                         "ignores loose .bed files at the top level.")
    ap.add_argument("--family", default="DARs_cell_type",
                    help="subdirectory name (default DARs_cell_type)")
    ap.add_argument("--min-log2fc", type=float, default=0.25,
                    help="minimum log2 fold change, in-group vs rest")
    ap.add_argument("--fdr", type=float, default=0.05)
    ap.add_argument("--max-regions", type=int, default=20_000,
                    help="cap per group, taken by effect size (default 20000). "
                         "Runtime of cisTarget scales with this.")
    ap.add_argument("--min-independent", type=int, default=5,
                    help="skip groups with fewer INDEPENDENT observations than "
                         "this (default 5). Raw metacell count overstates "
                         "support by the oversample factor -- at oversample=8, "
                         "Pro-Monocyte's 22 metacells are ~3 independent "
                         "observations. The divisor is read from "
                         "uns['glue_pairing']['oversample']; pass "
                         "--assume-oversample if the object predates that.")
    ap.add_argument("--diagnostics", type=Path, default=None,
                    help="pairing_diagnostics.csv from the pairing job. Its "
                         "independent_metacell_equiv column is the AUTHORITATIVE "
                         "per-group independent-observation count -- computed "
                         "from the limiting modality's cell count, which this "
                         "script cannot see. Strongly preferred over the "
                         "oversample-division fallback.")
    ap.add_argument("--assume-oversample", type=float, default=None,
                    help="fallback when no --diagnostics is given and the object "
                         "records no oversample factor (1.0 = treat metacells "
                         "as independent)")
    ap.add_argument("--block", type=int, default=20_000,
                    help="regions per streaming block (default 20000)")
    ap.add_argument("--dry-run", action="store_true",
                    help="report groups, counts and footprint; write nothing")
    args = ap.parse_args()

    if not args.h5mu.exists():
        sys.exit(f"ERROR: {args.h5mu} not found")

    try:
        import mudata
        from scipy.stats import rankdata, norm
    except ImportError as e:
        sys.exit(f"ERROR: needs mudata and scipy ({e})\n"
                 "       conda activate scplus-pairing")

    md = mudata.read(str(args.h5mu), backed=True)
    if args.modality not in md.mod:
        sys.exit(f"ERROR: modality '{args.modality}' not in "
                 f"{list(md.mod)}")
    A = md[args.modality]
    if args.group_key not in A.obs.columns:
        sys.exit(f"ERROR: '{args.group_key}' not in {args.modality}.obs "
                 f"(have: {list(A.obs.columns)})")

    labels = A.obs[args.group_key].astype(str).to_numpy()
    n_mc, n_reg = A.shape
    groups, counts = np.unique(labels, return_counts=True)

    # Filter on INDEPENDENT observations, not raw metacell count. With
    # oversample > 1 the metacells within a group share cells, so the raw count
    # overstates how much independent evidence a one-vs-rest test actually has
    # -- by exactly the oversample factor. Filtering on the raw count let every
    # group through: the smallest here is 22 metacells, but that is ~3
    # independent observations at oversample=8.
    # Preferred source: the diagnostics CSV the pairing job wrote. Its
    # independent_metacell_equiv is derived from the LIMITING modality's cell
    # count, which is not recoverable from the paired object -- dividing the
    # emitted metacell count by the oversample factor is close but not the same
    # number (Pro-Monocyte: floor(22/8)=2 vs the diagnostics' 1).
    indep_src = None
    indep_map = {}
    if args.diagnostics is not None:
        if not args.diagnostics.exists():
            sys.exit(f"ERROR: --diagnostics {args.diagnostics} not found")
        import csv as _csv
        with open(args.diagnostics) as fh:
            for row in _csv.DictReader(fh):
                key = row.get("group") or row.get("Group")
                val = row.get("independent_metacell_equiv")
                if key is None or val in (None, ""):
                    continue
                indep_map[str(key)] = int(float(val))
        if not indep_map:
            sys.exit(f"ERROR: {args.diagnostics} has no "
                     "'independent_metacell_equiv' column -- is it the "
                     "diagnostics CSV from the pairing job?")
        missing = [g for g in groups if g not in indep_map]
        if missing:
            sys.exit("ERROR: --diagnostics is missing these groups: "
                     f"{missing}\n       It does not describe this object. "
                     "Refusing to guess.")
        indep = np.array([indep_map[g] for g in groups], dtype=int)
        indep_src = f"{args.diagnostics.name} (independent_metacell_equiv)"
        ovs = None
    else:
        ovs = args.assume_oversample
        if ovs is None:
            ovs = (md.uns.get("glue_pairing", {}) or {}).get("oversample")
        if ovs is None:
            sys.exit("ERROR: cannot determine independent-observation counts.\n"
                     "       Pass --diagnostics pairing_diagnostics.csv (best), "
                     "or --assume-oversample.\n       Raw metacell counts "
                     "overstate support by the oversample factor, so\n"
                     "       filtering on them lets thinly-supported groups "
                     "through.")
        ovs = float(ovs)
        if ovs < 1.0:
            sys.exit(f"ERROR: oversample={ovs} < 1 makes no sense")
        indep = np.maximum(1, np.floor(counts / ovs)).astype(int)
        indep_src = (f"metacells / oversample={ovs:g} "
                     f"(APPROXIMATE -- prefer --diagnostics)")
    keep = indep >= args.min_independent
    skipped = [(g, int(c), int(i))
               for g, c, i, k in zip(groups, counts, indep, keep) if not k]

    chrom, start, end, bad = parse_regions(A.var_names)
    print("=" * 74)
    print("region sets from metacell DARs")
    print("=" * 74)
    print(f"  object          : {args.h5mu}  [{args.modality}]")
    print(f"  metacells       : {n_mc:,}")
    print(f"  regions         : {n_reg:,}  parsed {len(chrom):,}")
    if bad:
        print(f"  UNPARSEABLE     : {len(bad)} var_names, e.g. {bad[:3]}")
        sys.exit("ERROR: every var_name must be chr:start-end -- refusing to "
                 "write region sets from a partially parsed peak set.")
    print(f"  independence    : {indep_src}")
    print(f"  groups          : {len(groups)}  "
          f"({int(keep.sum())} kept, {len(skipped)} below "
          f"--min-independent={args.min_independent})")
    for g, c, i in skipped:
        print(f"    skip {g!r}: {c} metacells = ~{i} independent")
    blk_gb = args.block * n_mc * 4 / 1024**3
    print(f"  streaming block : {args.block:,} regions -> {blk_gb:.2f} GB "
          f"resident (full matrix would be "
          f"{n_reg * n_mc * 4 / 1024**3:.1f} GB)")

    if args.dry_run:
        print("\n--dry-run: nothing written.")
        print("=" * 74)
        return

    kept_groups = groups[keep]
    masks = {g: (labels == g) for g in kept_groups}
    # accumulate per group
    eff = {g: np.empty(n_reg, dtype=np.float32) for g in kept_groups}
    pv = {g: np.ones(n_reg, dtype=np.float64) for g in kept_groups}

    # RANK ONCE PER BLOCK, not once per group. scipy's mannwhitneyu re-ranks on
    # every call, so calling it 24 times per block ranks the same data 24 times.
    # Measured: 1.84 ms per region per group that way -- 393,832 x 24 = 4.8 h,
    # against a 4 h walltime, so the job would have been killed with nothing
    # written. Ranking once and deriving each group's U from its rank sum is
    # 20x faster (~14 min) and agrees with scipy to 5e-08, ties included.
    for i0 in range(0, n_reg, args.block):
        i1 = min(i0 + args.block, n_reg)
        blk = A[:, i0:i1].to_memory().X
        blk = np.asarray(blk.todense()) if hasattr(blk, "todense") else np.asarray(blk)
        Rk = rankdata(blk, axis=0)
        tie = _tie_terms(blk)
        for g, m in masks.items():
            a_mean = blk[m].mean(0)
            b_mean = blk[~m].mean(0)
            # log2 fold change on a pseudocount, since accessibility can be 0
            eff[g][i0:i1] = (np.log2(a_mean + 1e-6) - np.log2(b_mean + 1e-6))
            pv[g][i0:i1] = _mwu_from_ranks(Rk, m, tie, norm)
        if i0 % (args.block * 5) == 0:
            print(f"    ... {i1:,}/{n_reg:,} regions")

    outdir = args.out_dir / args.family
    outdir.mkdir(parents=True, exist_ok=True)
    summary = []
    print()
    for g in kept_groups:
        q = bh(pv[g])
        sel = (eff[g] >= args.min_log2fc) & (q <= args.fdr)
        idx = np.flatnonzero(sel)
        if idx.size > args.max_regions:
            idx = idx[np.argsort(-eff[g][idx])[:args.max_regions]]
            capped = True
        else:
            capped = False
        idx = idx[np.lexsort((start[idx], chrom[idx]))]
        safe = re.sub(r"[^0-9A-Za-z_.-]+", "_", g)
        path = outdir / f"{safe}.bed"
        with open(path, "w") as fh:
            for j in idx:
                fh.write(f"{chrom[j]}\t{start[j]}\t{end[j]}\t"
                         f"{safe}_{j}\t{eff[g][j]:.4f}\t.\n")
        summary.append(dict(group=g, file=path.name, n_regions=int(idx.size),
                            capped=bool(capped),
                            n_metacells=int(masks[g].sum()),
                            n_independent=int(
                                dict(zip(groups, indep))[g])))
        print(f"  {g:<26} {idx.size:>6,} regions"
              + ("  (capped)" if capped else ""))

    (args.out_dir / f"{args.family}.summary.json").write_text(
        json.dumps(dict(source=str(args.h5mu), group_key=args.group_key,
                        min_log2fc=args.min_log2fc, fdr=args.fdr,
                        max_regions=args.max_regions,
                        independence_source=indep_src,
                        min_independent=args.min_independent,
                        skipped_groups=[dict(group=g, n_metacells=c,
                                             n_independent=i)
                                        for g, c, i in skipped],
                        sets=summary), indent=2))
    empty = [s["group"] for s in summary if s["n_regions"] == 0]
    print()
    print(f"  wrote {len(summary)} BEDs to {outdir}")
    if empty:
        print(f"  WARNING: {len(empty)} group(s) produced ZERO regions: "
              f"{empty}")
        print("           A zero-region BED contributes no eRegulon. Loosen "
              "--min-log2fc\n           or --fdr, or accept that those groups "
              "are undetectable here.")
    print(f"  point the config at the PARENT: region_set_folder: "
          f"{args.out_dir}")
    print("=" * 74)


if __name__ == "__main__":
    main()
