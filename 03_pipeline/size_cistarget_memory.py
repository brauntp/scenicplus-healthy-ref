#!/usr/bin/env python3
"""
Size --mem and --cpus-per-task for motif_enrichment_cistarget, from the actual
region sets and the actual database.

WHY THIS EXISTS
---------------
`motif_enrichment_cistarget` was OOM-killed with no exit status while
`motif_enrichment_dem` -- which loads a database the same way -- completed. My
earlier sizing assumed ONE database load per rule. The code says otherwise
(scenicplus/cli/commands.py::run_motif_enrichment_cistarget):

    cistarget_results = joblib.Parallel(n_jobs=n_cpu, temp_folder=temp_dir)(
        joblib.delayed(_run_cistarget_single_region_set)(
            ..., cistarget_db_fname=cistarget_db_fname, ...)
        for key in region_set_dict)

and `_run_cistarget_single_region_set` constructs its own `cisTargetDatabase`
from that filename. joblib's default backend is loky -- separate PROCESSES -- so
there is no shared mapping: **n_cpu region sets are in flight at once, each
holding its own slice of the database.**

That also explains why DEM survived: its per-worker slice is a foreground set
plus a capped background, not a whole cell type's DAR set.

The peak is therefore driven by the LARGEST region sets, not by the total, and
the lever is `--cpus-per-task` as much as `--mem`. This script reads the real
.bed files and reports both.

Read-only, and cheap BY DEFAULT: it counts lines in the .bed files and does
arithmetic. Nothing here touches the 33 GB database unless you pass --db, and
even then it is one bounded column read (~131 KB) rather than a full map.

That caveat exists because the first version of this script got it wrong. It
used `pa.ipc.open_file(...).get_batch(i).num_rows` to count motifs, with
`read_table(columns=[])` as a "safe" fallback. Measured on a 3.04 GB fixture in
a fresh process, those peak at 3341 MB and 3508 MB respectively -- they map the
whole file -- against 251 MB for a single-column read. On the real 33 GB
rankings database that is tens of GB resident on a shared login node, which is
exactly the slowness that prompted this fix.

Usage:
    python 03_pipeline/size_cistarget_memory.py \
        --region-set-folder 01_atac/region_sets \
        --db resources/cistarget_db/hg38_screen_v10_clust.regions_vs_motifs.rankings.feather
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

BYTES_PER_RANK = 4          # int32 rankings
TRANSIENT_FACTOR = 2.0      # pyarrow Table + pandas DataFrame both briefly live

# THE DOMINANT TERM, and the one this script originally missed entirely.
#
# cisTarget.run_ctx() calls ctxcore's recovery(), which pre-allocates
#
#     rccs = np.empty(shape=(n_features, rank_threshold))
#
# with n_features = ALL motifs in the database and rank_threshold =
# int(ctx_rank_threshold * total_regions_in_database). np.empty's default dtype
# is float64. It is then wrapped in df_rccs, a DataFrame with a MultiIndex over
# every curve point, so two copies are live.
#
# Note what this does NOT depend on: the size of the region set. That is why
# rewriting the sets to an equal 5,000 regions did not stop the OOM.
RECOVERY_BYTES = 8          # np.empty -> float64
RECOVERY_COPIES = 2.0       # rccs + df_rccs


def count_bed(p: Path) -> int:
    n = 0
    with open(p, "rb") as fh:
        for line in fh:
            if line.strip() and not line.startswith((b"#", b"track", b"browser")):
                n += 1
    return n


def n_motifs_in_db(db: Path) -> int | None:
    """Row count of the rankings feather = number of motifs.

    MEASURED, not assumed. Three ways to get this number, on a 3.04 GB / 400k-
    column fixture in a fresh process:

        approach                    time    peak RSS
        get_batch(i).num_rows       1.43 s    3341 MB   <- maps the whole file
        read_table(columns=[])      0.82 s    3508 MB   <- worse
        read_table(columns=[one])   0.40 s     251 MB   <- bounded
        schema only                 0.20 s     148 MB   (but gives COLUMNS,
                                                         not rows)

    The first two scale with FILE SIZE: on the real 33 GB rankings database that
    is tens of GB resident. This script's earlier version used `get_batch` with
    `columns=[]` as its "safe" fallback -- both of them the same hazard, which is
    why the run the user reported took a long time on the login node.

    Reading ONE column is bounded: the schema plus 32,765 int32 values, i.e.
    ~131 KB of data regardless of how wide the file is.
    """
    try:
        import pyarrow as pa
        import pyarrow.feather as pf
    except ImportError:
        return None
    try:
        # Schema first (cheap), then exactly one column.
        with pa.memory_map(str(db), "rb") as src:
            names = pa.ipc.open_file(src).schema.names
        if not names:
            return None
        # The motif/track id column is a string column; any DATA column gives
        # the same row count and is smaller. Prefer a non-id column when the
        # schema exposes one.
        pick = next((n for n in names if n not in ("motifs", "tracks")), names[0])
        return pf.read_table(db, columns=[pick]).num_rows
    except Exception:
        return None


def human(gb: float) -> str:
    return f"{gb:,.1f} GB"


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--region-set-folder", type=Path, required=True)
    ap.add_argument("--db", type=Path, default=None,
                    help="OPTIONAL rankings feather, to confirm the motif count "
                         "from the file itself. NOT needed: the count is a fixed "
                         "property of the motif collection (v10nr_clust = 32765) "
                         "and --n-motifs covers it. Reading it costs one bounded "
                         "column read, but the file may live on a filesystem "
                         "where even that is slow -- skip it on a login node.")
    ap.add_argument("--n-motifs", type=int, default=32765,
                    help="Fallback motif count if --db is not given or not "
                         "readable (v10nr_clust: 32765).")
    ap.add_argument("--cpus", type=int, nargs="+", default=[1, 2, 4, 8, 16],
                    help="Worker counts to tabulate.")
    ap.add_argument("--total-db-regions", type=int, default=1_837_304,
                    help="Regions in the cisTarget database (hg38 SCREEN v10 "
                         "= 1,837,304). Sets the recovery-curve length, which "
                         "is the dominant memory term.")
    ap.add_argument("--rank-threshold", type=float, default=0.05,
                    help="params_motif_enrichment.ctx_rank_threshold from the "
                         "config. Curve length = this x --total-db-regions.")
    ap.add_argument("--auc-threshold", type=float, default=0.005,
                    help="params_motif_enrichment.ctx_auc_threshold. Sets the "
                         "FLOOR on --rank-threshold: ctxcore asserts "
                         "rank_cutoff <= curve length.")
    ap.add_argument("--mem-limit-gb", type=float, default=128,
                    help="The allocation to judge against, i.e. what --mem says "
                         "in slurm/scenicplus.sbatch. Default 128.")
    args = ap.parse_args()

    if not args.region_set_folder.is_dir():
        sys.exit(f"ERROR: not a directory: {args.region_set_folder}")

    beds = sorted(args.region_set_folder.rglob("*.bed"))
    if not beds:
        sys.exit(f"ERROR: no .bed files under {args.region_set_folder}\n"
                 f"       cisTarget globs <folder>/<family>/*.bed")

    n_motifs = args.n_motifs
    src = "fallback"
    if args.db and args.db.is_file():
        got = n_motifs_in_db(args.db)
        if got:
            n_motifs, src = got, "read from the database"

    sizes = [(p.parent.name + "/" + p.stem, count_bed(p)) for p in beds]
    sizes.sort(key=lambda t: -t[1])

    print("=" * 74)
    print("motif_enrichment_cistarget memory sizing")
    print("=" * 74)
    print(f"  region sets : {len(sizes)}")
    print(f"  motifs      : {n_motifs:,}  ({src})")
    print()

    print("-- region sets, largest first --------------------------------------")
    for name, n in sizes[:8]:
        print(f"  {name:<44} {n:>9,} regions"
              f"   {human(n*n_motifs*BYTES_PER_RANK/1024**3):>10}")
    if len(sizes) > 8:
        print(f"  ... {len(sizes)-8} more, down to {sizes[-1][1]:,} regions")
    print()

    curve_pts = int(args.rank_threshold * args.total_db_regions)
    rec_gb = (n_motifs * curve_pts * RECOVERY_BYTES
              * RECOVERY_COPIES / 1024**3)

    print("-- the recovery-curve term (dominates; independent of set size) -----")
    print(f"  ctx_rank_threshold {args.rank_threshold} x "
          f"{args.total_db_regions:,} db regions")
    print(f"    = {curve_pts:,} curve points per motif")
    print(f"  rccs: {n_motifs:,} motifs x {curve_pts:,} points, float64 "
          f"= {human(rec_gb/RECOVERY_COPIES)}")
    print(f"  + df_rccs (DataFrame copy with a MultiIndex over every point)")
    print(f"  = {human(rec_gb)} PER WORKER, whatever the region set contains")
    print()

    print("-- concurrent peak: joblib runs n_cpu sets at once ------------------")
    print("  Each worker is a separate PROCESS: its own database slice AND its")
    print("  own recovery matrix. Nothing is shared.")
    print()
    print(f"  {'--cpus':>7}  {'db slice':>10}  {'recovery':>10}  "
          f"{'total':>10}  {'suggested --mem':>16}")
    rows = []
    for c in args.cpus:
        top = [n for _, n in sizes[:c]]
        db_gb = sum(top) * n_motifs * BYTES_PER_RANK / 1024**3 * TRANSIENT_FACTOR
        rc_gb = rec_gb * c
        peak = db_gb + rc_gb
        mem = int(peak * 1.3) + 8
        rows.append((c, db_gb, peak, mem))
        print(f"  {c:>7}  {human(db_gb):>10}  {human(rc_gb):>10}  "
              f"{human(peak):>10}  {str(mem)+'G':>16}")
    print()
    print("  (--mem adds 30% + 8 GB for the interpreter, motif annotations and")
    print("   the result objects.)")
    print()

    # The floor is NOT simply auc_threshold. scenicplus passes
    # int(rank_threshold * total) as the curve length while ctxcore computes
    # rank_cutoff = round(auc_threshold * total) and asserts
    # rank_cutoff <= curve_length. At rank_threshold == auc_threshold those
    # differ by the truncation: int(0.005*1837304) = 9186 vs round(...) = 9187,
    # so the assert FAILS. Derive the smallest value that actually passes.
    _cut = round(args.auc_threshold * args.total_db_regions)
    floor = args.auc_threshold
    while int(floor * args.total_db_regions) < _cut:
        floor = round(floor + 0.0005, 4)
    print("-- lowering ctx_rank_threshold: the real lever ----------------------")
    print(f"  FLOOR: ctx_rank_threshold >= {floor} (NOT simply "
          f"ctx_auc_threshold = {args.auc_threshold}).")
    print(f"  scenicplus passes int(rank_threshold x total) as the curve length;")
    print(f"  ctxcore computes rank_cutoff = round(auc_threshold x total) = "
          f"{_cut:,} and")
    print(f"  asserts rank_cutoff <= curve length. At rank_threshold = "
          f"{args.auc_threshold} those are")
    print(f"  {int(args.auc_threshold*args.total_db_regions):,} vs {_cut:,} -- "
          f"the assert FAILS on the truncation.")
    print()
    print(f"  {'rank_thr':>9}  {'points':>10}  {'per worker':>11}  " +
          "  ".join(f"{'x'+str(c):>7}" for c in args.cpus))
    grid = sorted({0.05, 0.02, 0.01, floor}, reverse=True)
    for rt in grid:
        if rt < floor:
            continue
        # Belt-and-braces: every value printed must actually pass ctxcore's
        # assert, or the table is recommending a crash.
        assert int(rt * args.total_db_regions) >= _cut, rt
        pts = int(rt * args.total_db_regions)
        one = n_motifs * pts * RECOVERY_BYTES * RECOVERY_COPIES / 1024**3
        # 4 decimals, not 3: the derived floor can be 0.0055, and printing it as
        # "0.005" would name a value that fails ctxcore's assert.
        line = f"  {rt:>9.4f}  {pts:>10,}  {one:>10.1f}G  "
        line += "  ".join(f"{one*c:>6.0f}G" for c in args.cpus)
        print(line)
    print()
    print("  What this does NOT change: NES, AUC, and which motifs are called")
    print("  enriched. The AUC integrates rccs[:, :rank_cutoff] where")
    print(f"  rank_cutoff = {round(floor*args.total_db_regions)-1:,}, fixed by "
          f"ctx_auc_threshold.")
    print("  Any rank_threshold at or above the floor still covers it.")
    print()
    print("  What it DOES change, and this is more than cosmetic:")
    print("    rank_at_max = argmax(rcc - avg2stdrcc) over the curve, and")
    print("    leading_edge() returns every region ranked at or before it as a")
    print("    MOTIF HIT. Those hits become the cistromes, which are the input")
    print("    to prepare_menr and hence to eGRN construction. Truncating the")
    print("    curve caps rank_at_max, so an enriched motif whose recovery")
    print("    keeps climbing past the cut loses the hits beyond it -- smaller")
    print("    cistromes, not just a different diagnostic.")
    print()
    print("    The motif CALLS (NES/AUC) are unaffected, so this does not")
    print("    invent or lose regulons. It can shrink the region membership of")
    print("    the ones it finds.")
    print()

    print("-- if the sets were equal-sized (--top-n) --------------------------")
    print("  With one shared size the peak is a PRODUCT rather than a sum over")
    print("  the largest few, so it is predictable and does not move when a")
    print("  cell type gains regions. 01_cistopic/choose_dar_threshold.py")
    print("  --top-n N --write rewrites the sets that way.")
    print()
    print(f"  {'top-n':>8}  {'per worker':>11}  " +
          "  ".join(f"{'x'+str(c):>7}" for c in args.cpus))
    for N in (2_000, 5_000, 10_000, 20_000):
        one = N * n_motifs * BYTES_PER_RANK / 1024**3 * TRANSIENT_FACTOR
        line = f"  {N:>8,}  {one:>10.1f}G  "
        line += "  ".join(f"{one*c:>6.0f}G" for c in args.cpus)
        print(line)
    print()
    print("  (per worker includes the transient 2x. Note the trade the equal-N")
    print("   mode makes, which is scientific rather than computational:")
    print("   membership becomes rank-based, so a region admitted for one cell")
    print("   type would be rejected for another -- see docs/REGION_SETS.md.)")
    print()

    print("-- what this means -------------------------------------------------")
    print("  The recovery term dominates and is n_cpu-fold, so BOTH levers")
    print("  work, and they multiply:")
    print()
    budget = args.mem_limit_gb * 0.9
    best = None
    for rt in (0.05, 0.02, 0.01, floor):
        if rt < floor:
            continue
        pts = int(rt * args.total_db_regions)
        one = n_motifs * pts * RECOVERY_BYTES * RECOVERY_COPIES / 1024**3
        ok = [c for c in args.cpus
              if one * c + (sum(n for _, n in sizes[:c]) * n_motifs
                            * BYTES_PER_RANK / 1024**3 * TRANSIENT_FACTOR)
              <= budget]
        top = max(ok) if ok else None
        print(f"    ctx_rank_threshold {rt:<6} -> largest --cpus-per-task that "
              f"fits {args.mem_limit_gb:g} GB: "
              f"{top if top else 'NONE (even 1 worker exceeds it)'}")
        if top and (best is None or top > best[1]):
            best = (rt, top)
    print()
    if best:
        rt, c = best
        print(f"  Best throughput inside {args.mem_limit_gb:g} GB: "
              f"ctx_rank_threshold {rt}, --cpus-per-task {c}.")
    print(f"  There are {len(sizes)} region sets, so more than {len(sizes)} "
          f"workers buys nothing.")
    print()
    print("  WHICH LEVER FIRST: --cpus-per-task, not the threshold.")
    print("  Fewer workers costs wall time in this rule and changes NOTHING")
    print("  about the result. Lowering ctx_rank_threshold is cheaper in time")
    print("  but can shrink cistrome membership (see above), and cistromes feed")
    print("  eGRN construction. Reach for the threshold only if the wall time")
    print("  at an acceptable worker count does not fit the walltime.")
    print()
    n_sets = len(sizes)
    print(f"  Wave count at each worker level ({n_sets} region sets, run in")
    print(f"  batches of --cpus-per-task):")
    for c in args.cpus:
        peak = rec_gb * c + (sum(n for _, n in sizes[:c]) * n_motifs
                             * BYTES_PER_RANK / 1024**3 * TRANSIENT_FACTOR)
        fits = "fits" if peak <= budget else "OOM"
        print(f"    --cpus-per-task {c:>2}: {-(-n_sets // c):>2} waves, "
              f"peak {human(peak):>9}  [{fits}]")
    print("=" * 74)


if __name__ == "__main__":
    main()
