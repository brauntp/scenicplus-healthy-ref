#!/usr/bin/env python3
"""
Choose --min-log2fc from the measured effect-size distribution, then rewrite the
BEDs. No second pass over the accessibility matrix.

WHY THIS EXISTS
---------------
The first real run capped 20 of 21 groups at --max-regions 20,000 (409,796 BED
lines total, 1.04x the entire 393,832-peak set). When almost every group caps,
the cap is choosing the regions and the thresholds are not: each set is an
arbitrary top-N slice by effect size, and 20,000 has no biological meaning.

The FDR filter cannot bind at this scale. With 25,323 metacells a 2% shift in
probability-of-superiority gives p ~ 1e-2 and a 5% shift p ~ 4e-8; oversampling
at 8x makes the p-values anticonservative on top of that (docs/
oversample_null.csv). So --min-log2fc is the only filter doing work, and 0.25
was a guess made before any distribution was in hand.

This reads the .npz that region_sets_from_metacells.py --dump-stats writes and
answers, for a grid of thresholds: how many regions each group keeps, how many
still cap, and how much of the peak set the union covers. Then it writes the
BEDs at your chosen threshold directly from the dump.

Usage
-----
    # what would each threshold give?
    python 01_cistopic/choose_dar_threshold.py --stats 01_atac/dar_stats.npz

    # rewrite the BEDs at a chosen threshold
    python 01_cistopic/choose_dar_threshold.py --stats 01_atac/dar_stats.npz \\
        --min-log2fc 1.0 --out-dir 01_atac/region_sets --write

WHAT TO AIM FOR
---------------
No hard rule, but two anchors. cisTarget runtime scales with total BED lines, so
409,796 across 21 sets is the expensive end. And a region set that is 5% of all
peaks is not "differentially accessible in this cell type" in any useful sense --
it is most of the accessible genome. A few thousand regions per group, with no
group capped, is the regime the method was designed for.
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


def refine(sig, eff, lo, hi, max_regions, step=0.01):
    """Lowest threshold in [lo, hi] that caps no group and empties none.

    The coarse grid overshoots: on a 0.5 step it reported 1.50 where 1.05 also
    worked, discarding regions for nothing. Both the grid path and the
    auto-extend path refine through here so the reported threshold is the real
    minimum, not the first grid point that happened to clear.
    """
    import numpy as _np
    for t in _np.arange(lo, hi + 1e-9, step):
        n = (sig & (eff >= t)).sum(1)
        if n.min() > 0 and n.max() <= max_regions:
            return float(t), n
    return None, None


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--stats", type=Path, required=True,
                    help=".npz from region_sets_from_metacells.py --dump-stats")
    ap.add_argument("--fdr", type=float, default=0.05)
    ap.add_argument("--max-regions", type=int, default=20_000)
    ap.add_argument("--min-usable", type=int, default=2_000,
                    help="a region set smaller than this is reported as "
                         "STARVED, not merely small (default 2000). cisTarget's "
                         "AUC window is ctx_auc_threshold (0.005) of the "
                         "1,837,304-region database, so a set of N regions puts "
                         "N*0.005 regions in that window under the null: 7 "
                         "regions gives 0.04 expected and 473 gives 2.4, at "
                         "which point a NES of 3.0 is noise. ~2,000 is where "
                         "the recovery curve stabilises.")
    ap.add_argument("--top-n", type=int, default=None,
                    help="ignore --min-log2fc and take the top N regions per "
                         "group by effect size (FDR still applied). Use when no "
                         "single threshold serves every group -- which is the "
                         "usual outcome when effect-size distributions differ "
                         "by more than ~20x across groups.")
    ap.add_argument("--grid", type=float, nargs="+",
                    default=[0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0])
    ap.add_argument("--min-log2fc", type=float, default=None,
                    help="chosen threshold; required with --write")
    ap.add_argument("--out-dir", type=Path, default=None,
                    help="PARENT region-set directory (a --family subdirectory "
                         "is created inside it)")
    ap.add_argument("--family", default="DARs_cell_type")
    ap.add_argument("--write", action="store_true",
                    help="rewrite the BEDs at --min-log2fc")
    args = ap.parse_args()

    if not args.stats.exists():
        sys.exit(f"ERROR: {args.stats} not found\n"
                 "       Re-run region_sets_from_metacells.py with "
                 "--dump-stats to produce it.")
    z = np.load(args.stats, allow_pickle=False)
    groups = [str(g) for g in z["groups"]]
    eff, q = z["eff"], z["q"]
    chrom, start, end = z["chrom"], z["start"], z["end"]
    n_grp, n_reg = eff.shape

    print("=" * 78)
    print("DAR threshold sweep")
    print("=" * 78)
    print(f"  stats     : {args.stats}")
    print(f"  groups    : {n_grp}   regions: {n_reg:,}")
    print(f"  FDR       : {args.fdr}   cap: {args.max_regions:,}")
    print()

    sig = q <= args.fdr
    print(f"{'min_log2fc':>11}{'median/grp':>12}{'min':>9}{'max':>9}"
          f"{'capped':>8}{'starved':>9}{'total BED':>12}{'union % peaks':>15}")
    print("-" * 78)
    table = []
    for t in args.grid:
        sel = sig & (eff >= t)
        n = sel.sum(1)
        kept = np.minimum(n, args.max_regions)
        union = sel.any(0).sum()
        row = dict(min_log2fc=float(t), median=int(np.median(n)),
                   min=int(n.min()), max=int(n.max()),
                   n_capped=int((n > args.max_regions).sum()),
                   n_zero=int((n < args.min_usable).sum()),
                   total=int(kept.sum()),
                   union_frac=float(union / n_reg))
        table.append(row)
        print(f"{t:>11.2f}{row['median']:>12,}{row['min']:>9,}{row['max']:>9,}"
              f"{row['n_capped']:>8}{row['n_zero']:>9}{row['total']:>12,}"
              f"{row['union_frac']:>14.1%}")

    print()
    # The first threshold that caps nothing and starves nothing.
    clean = [r for r in table if r["n_capped"] == 0 and r["n_zero"] == 0]
    if clean:
        rec = clean[0]
        # The grid point that cleared is an upper bound, not the answer: refine
        # down to the true minimum so no regions are discarded for nothing.
        below = [x for x in args.grid if x < rec["min_log2fc"]]
        lo = max(below) if below else min(args.grid)
        t_ref, n_ref = refine(sig, eff, lo, rec["min_log2fc"],
                              args.max_regions)
        if t_ref is not None and t_ref < rec["min_log2fc"] - 1e-9:
            union = (sig & (eff >= t_ref)).any(0).sum()
            head = 1 - n_ref.max() / args.max_regions
            print(f"  FLOOR   {t_ref:.2f}  -- lowest that caps nothing: "
                  f"median {int(np.median(n_ref)):,}/group, "
                  f"{int(n_ref.sum()):,} lines, union {union / n_reg:.1%}")
            if head < 0.10:
                print(f"          but the largest group is {head:.0%} from the "
                      f"cap, so this is the")
                print(f"          boundary, not a safe operating point.")
            print(f"  Pick from the table above by what the region sets should "
                  f"MEAN, not by")
            print(f"  the floor: {rec['min_log2fc']:.2f} gives "
                  f"{rec['median']:,}/group with "
                  f"{1 - max(r['max'] for r in table if r['min_log2fc'] == rec['min_log2fc']) / args.max_regions:.0%}"
                  f" headroom.")
            rec = dict(min_log2fc=t_ref)
        else:
            print(f"  lowest threshold with NO capped and NO empty group: "
                  f"{rec['min_log2fc']:.2f}")
            print(f"    median {rec['median']:,} regions/group, "
                  f"{rec['total']:,} BED lines total, "
                  f"union {rec['union_frac']:.1%} of peaks")
        print("  Lower keeps more marginal regions but lets the cap choose; "
              "higher is")
        print("  stricter but starves the thinner groups first.")
    else:
        # Auto-extend rather than telling the user to re-run with a wider grid:
        # the dump is already in memory and each extra point costs a boolean
        # reduction, so there is no reason to make this a second invocation.
        print("  No grid point both clears the cap and keeps every group "
              "non-empty -- extending automatically:")
        t = max(args.grid)
        rec = None
        for _ in range(24):
            t += 0.5
            sel = sig & (eff >= t)
            n = sel.sum(1)
            n_cap = int((n > args.max_regions).sum())
            n_zero = int((n == 0).sum())
            union = sel.any(0).sum()
            print(f"{t:>11.2f}{int(np.median(n)):>12,}{int(n.min()):>9,}"
                  f"{int(n.max()):>9,}{n_cap:>8}{n_zero:>6}"
                  f"{int(np.minimum(n, args.max_regions).sum()):>12,}"
                  f"{union / n_reg:>14.1%}")
            if n_cap == 0 and n_zero == 0:
                t_ref, n_ref = refine(sig, eff, t - 0.5, t, args.max_regions)
                if t_ref is None:
                    t_ref, n_ref = float(t), n
                un = (sig & (eff >= t_ref)).any(0).sum()
                rec = dict(min_log2fc=float(t_ref),
                           median=int(np.median(n_ref)),
                           total=int(np.minimum(n_ref, args.max_regions).sum()),
                           union_frac=float(un / n_reg))
                break
            if n_zero > 0:
                # Groups start emptying before the cap clears: the two
                # conditions cannot both be met, which is a finding about the
                # thin groups rather than a parameter to keep tuning.
                empt = [groups[i] for i in np.flatnonzero(n == 0)]
                print()
                if n_cap == 0:
                    # Caps cleared at the same coarse step the thin group
                    # emptied. Do NOT tell the user to retry at a finer step
                    # without checking one exists -- refine here and report
                    # whichever way it comes out.
                    # Search from the LOWEST grid point upward, not just the
                    # last coarse step: the starving group can empty below the
                    # window where caps clear, and a satisfying threshold, if
                    # one exists, may sit under it.
                    lo, hi = min(args.grid), t
                    hit, nf = refine(sig, eff, lo, hi, args.max_regions)
                    if hit is not None:
                        print(f"  refined between {lo:.2f} and {hi:.2f} at 0.01:")
                        print(f"  lowest threshold with NO capped and NO empty "
                              f"group: {hit:.2f}")
                        print(f"    median {int(np.median(nf)):,} regions/group, "
                              f"{int(nf.sum()):,} BED lines total")
                        rec = dict(min_log2fc=hit)
                    else:
                        print(f"  NO threshold exists that satisfies both. "
                              f"Refined {lo:.2f}-{hi:.2f} at 0.01 resolution:")
                        print(f"        {', '.join(empt[:5])} empt{'ies' if len(empt) == 1 else 'y'} "
                              f"before the broad groups stop capping.")
                        print("        The two conditions are incompatible on "
                              "this data -- that is a")
                        print("        finding about the thin groups, not a "
                              "parameter to keep tuning.")
                        print("        Options: accept the cap for the broad "
                              "groups (their sets are")
                        print("        top-N slices), drop the starving groups "
                          "via --min-independent, or")
                        print("        use per-group thresholds and record that "
                              "sizes are not comparable.")
                    break
                else:
                    print(f"  STOP: {n_zero} group(s) empty at {t:.2f} while "
                          f"{n_cap} group(s) still cap.")
                    print(f"        emptied: {', '.join(empt[:5])}")
                    print("        No single threshold serves every group: the "
                          "thin groups run out")
                    print("        before the broad ones stop capping. Either "
                          "accept the cap for the")
                    print("        broad groups, or use a per-group threshold "
                          "and record that the")
                    print("        sets are not comparable in size.")
                break
        if rec:
            print()
            print(f"  lowest threshold with NO capped and NO empty group: "
                  f"{rec['min_log2fc']:.2f}")
            print(f"    median {rec['median']:,} regions/group, "
                  f"{rec['total']:,} BED lines total, "
                  f"union {rec['union_frac']:.1%} of peaks")
        elif rec is None and t >= max(args.grid) + 12:
            print()
            print("  Extended to +12 log2fc without clearing the cap. Something "
                  "is wrong with")
            print("  the effect-size distribution -- inspect the dump directly "
                  "before proceeding.")

    if args.top_n is not None:
        print()
        print(f"  --top-n {args.top_n:,}: equal-sized sets, no shared threshold.")
        n_avail = (sig.sum(1))
        short = [(groups[i], int(n_avail[i]))
                 for i in np.flatnonzero(n_avail < args.top_n)]
        print(f"    groups with fewer than {args.top_n:,} FDR-passing regions: "
              f"{len(short)}")
        for g, n in short[:8]:
            print(f"      {g:<26}{n:>8,}")
        eff_at = np.array([np.sort(eff[i][sig[i]])[::-1][
            min(args.top_n, int(n_avail[i])) - 1] if n_avail[i] else np.nan
            for i in range(n_grp)])
        print(f"    implied log2fc cutoff spans "
              f"{np.nanmin(eff_at):.2f} to {np.nanmax(eff_at):.2f} across groups")
        print("    -- that spread is the cost: a region admitted in one group "
              "would be")
        print("       rejected in another, so set membership is rank-based, not "
              "effect-based.")

    if not args.write:
        print()
        print("  (sweep only -- pass --min-log2fc T --out-dir DIR --write to "
              "rewrite the BEDs)")
        print("=" * 78)
        return

    if args.out_dir is None:
        sys.exit("ERROR: --write needs --out-dir")
    if args.min_log2fc is None and args.top_n is None:
        sys.exit("ERROR: --write needs either --min-log2fc T or --top-n N")

    outdir = args.out_dir / args.family
    outdir.mkdir(parents=True, exist_ok=True)
    # Remove BEDs from the previous threshold: SCENIC+ globs the directory, so a
    # leftover file becomes a silent extra enrichment run at the old setting.
    stale = sorted(outdir.glob("*.bed"))
    for f in stale:
        f.unlink()
    if stale:
        print(f"\n  removed {len(stale)} BED(s) from the previous threshold")

    print()
    summary = []
    for i, g in enumerate(groups):
        if args.top_n is not None:
            idx = np.flatnonzero(sig[i])
            capped = False
            if idx.size > args.top_n:
                idx = idx[np.argsort(-eff[i][idx])[:args.top_n]]
        else:
            sel = sig[i] & (eff[i] >= args.min_log2fc)
            idx = np.flatnonzero(sel)
            capped = idx.size > args.max_regions
            if capped:
                idx = idx[np.argsort(-eff[i][idx])[:args.max_regions]]
        idx = idx[np.lexsort((start[idx], chrom[idx]))]
        safe = re.sub(r"[^0-9A-Za-z_.-]+", "_", g)
        path = outdir / f"{safe}.bed"
        with open(path, "w") as fh:
            for j in idx:
                fh.write(f"{chrom[j]}\t{start[j]}\t{end[j]}\t"
                         f"{safe}_{j}\t{eff[i][j]:.4f}\t.\n")
        summary.append(dict(group=g, file=path.name, n_regions=int(idx.size),
                            capped=bool(capped)))
        print(f"  {g:<26} {idx.size:>7,} regions"
              + ("  (capped)" if capped else ""))

    (args.out_dir / f"{args.family}.summary.json").write_text(
        json.dumps(dict(source=str(args.stats), min_log2fc=args.min_log2fc,
                        top_n=args.top_n, min_usable=args.min_usable,
                        fdr=args.fdr, max_regions=args.max_regions,
                        rewritten_by="01_cistopic/choose_dar_threshold.py",
                        sets=summary), indent=2))
    tot = sum(s["n_regions"] for s in summary)
    print()
    print(f"  wrote {len(summary)} BEDs, {tot:,} lines total, to {outdir}")
    n_cap = sum(1 for s in summary if s["capped"])
    if n_cap:
        print(f"  WARNING: {n_cap} group(s) still capped at "
              f"{args.max_regions:,} -- raise --min-log2fc further")
    print("=" * 78)


if __name__ == "__main__":
    main()
