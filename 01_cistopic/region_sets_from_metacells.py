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


def _ranks_and_ties(X):
    """Average ranks (SAME dtype as X) and sum(t**3 - t) per column, from ONE
    argsort.

    Why not scipy.stats.rankdata + a separate np.sort:

      * rankdata returns float64 -- 2x a float32 input -- and argsorts int64
        internally. On a 20,000-region block of 25,323 metacells (1.89 GB) that
        peaked at 9.59 GB measured, and the job was OOM-killed at 16G (exit 137)
        against my "flat ~2 GB" projection.
      * computing tie terms from a second np.sort materialised another full copy
        of the block.

    Doing both from one stable argsort, keeping ranks in the input dtype, peaks
    at 2.45 GB for a 4,000-region block and gives ranks IDENTICAL to
    rankdata (max |diff| 0.0, ties included).
    """
    n = X.shape[0]
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
    tie = np.zeros(X.shape[1], dtype=np.float64)
    for j in range(X.shape[1]):
        cnt = np.bincount(grp[:, j], minlength=int(G[j]) + 1)[1:]
        ends = np.cumsum(cnt)
        starts = ends - cnt
        ranks_sorted[:, j] = np.repeat(
            ((starts + ends + 1) / 2.0).astype(X.dtype), cnt)
        c = cnt.astype(np.float64)
        tie[j] = float(((c ** 3) - c).sum())
    del grp
    R = np.empty_like(X)
    np.put_along_axis(R, order, ranks_sorted, axis=0)
    return R, tie


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
    # float64 accumulator: ranks are stored in the block's dtype (float32) to
    # halve memory, and summing ~3k float32 ranks near 2.5e4 loses precision.
    U = Rk[mask].sum(0, dtype=np.float64) - n1 * (n1 + 1) / 2.0
    mu = n1 * n2 / 2.0
    var = n1 * n2 / 12.0 * ((N + 1) - tie / (N * (N - 1)))
    sd = np.sqrt(np.maximum(var, 1e-12))
    return norm.sf((U - mu - 0.5) / sd)


def _build_parser():
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
                         "Pro-Monocyte's 22 metacells are floor(22/8)=2 "
                         "independent observations by the fallback rule, and 1 "
                         "by the diagnostics. Prefer --diagnostics; the "
                         "fallback divisor comes from "
                         "uns['glue_pairing']['oversample'].")
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
    ap.add_argument("--block", type=int, default=4_000,
                    help="regions per streaming block (default 4000). Peak RSS "
                         "is 5-6.5x the block's float32 size: the block, its "
                         "ranks, a stable argsort (int64) and the gathered "
                         "sorted copy. 4,000 x 25,323 measured 2.45 GB; 20,000 "
                         "measured 9.59 GB and was OOM-killed at 16G.")
    ap.add_argument("--dump-stats", type=Path, default=None,
                    help="write per-group effect sizes and q-values to this "
                         ".npz alongside the BEDs. 21 groups x 393,832 regions "
                         "of float32, twice (effects and q-values), is 63 MB "
                         "uncompressed and less on disk. It makes every threshold "
                         "question answerable OFFLINE instead of costing "
                         "another 14-minute pass. Strongly recommended on the "
                         "first run for a new object.")
    ap.add_argument("--dry-run", action="store_true",
                    help="report groups, counts and footprint; write nothing")
    ap.add_argument("--self-test", action="store_true",
                    help="verify the help text's worked example against the "
                         "code's own arithmetic, then exit")
    return ap

def _self_test():
    """Check the --min-independent help text against what the code actually does.

    This exists because that help string said "~3 independent" while the same
    script printed "~2 independent" for the same input: an example in prose
    drifted from the arithmetic beside it, and reading the text could not catch
    it. So compute the example rather than trusting it.
    """
    import re
    ap = _build_parser()
    txt = " ".join(a.help or "" for a in ap._actions)
    txt = re.sub(r"\s+", " ", txt)
    m = re.search(r"(\d+) metacells are floor\((\d+)/(\d+)\)=(\d+)", txt)
    if not m:
        print("SELF-TEST FAIL: no worked example found in --min-independent help")
        return 1
    mc, num, div, stated = (int(g) for g in m.groups())
    if num != mc:
        print(f"SELF-TEST FAIL: example says {mc} metacells but divides {num}")
        return 1
    computed = int(max(1, np.floor(mc / div)))
    if computed != stated:
        print(f"SELF-TEST FAIL: help says floor({mc}/{div})={stated}, "
              f"the fallback rule gives {computed}")
        return 1
    print(f"SELF-TEST PASS: help example floor({mc}/{div})={stated} matches the "
          f"fallback rule")
    return 0


def main():
    ap = _build_parser()
    if "--self-test" in sys.argv:
        sys.exit(_self_test())
    args = ap.parse_args()

    if not args.h5mu.exists():
        sys.exit(f"ERROR: {args.h5mu} not found")

    try:
        import mudata
        from scipy.stats import norm
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
    print(f"  streaming block : {args.block:,} regions = {blk_gb:.2f} GB, "
          f"peak ~{blk_gb * 6.5:.1f} GB "
          f"(measured factor; full matrix is "
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
        Rk, tie = _ranks_and_ties(blk)
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

    n_capped = sum(1 for x in summary if x["capped"])
    if n_capped:
        print()
        print(f"  WARNING: {n_capped} of {len(summary)} groups hit "
              f"--max-regions={args.max_regions:,}.")
        if n_capped > len(summary) / 2:
            print("           When most groups cap, the CAP is selecting the "
                  "regions, not your")
            print("           thresholds -- each set is an arbitrary top-N "
                  "slice by effect size,")
            print("           and the number 20,000 has no biological meaning. "
                  "At 25,323")
            print("           metacells an FDR filter cannot bind (a 2% shift "
                  "in probability-")
            print("           of-superiority gives p~1e-2, a 5% shift p~4e-8), "
                  "so --min-log2fc")
            print("           is the only real filter. Raise it until the caps "
                  "clear:")
            print("           re-run with --dump-stats to choose it from the "
                  "measured")
            print("           distribution rather than by guessing.")

    if args.dump_stats is not None:
        args.dump_stats.parent.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(
            args.dump_stats,
            groups=np.array(list(kept_groups)),
            eff=np.vstack([eff[g] for g in kept_groups]).astype(np.float32),
            q=np.vstack([bh(pv[g]) for g in kept_groups]).astype(np.float32),
            chrom=chrom, start=start, end=end)
        print(f"  stats   : {args.dump_stats} "
              f"({args.dump_stats.stat().st_size / 1024**2:.1f} MB) -- "
              f"threshold sweeps need no re-run")

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
