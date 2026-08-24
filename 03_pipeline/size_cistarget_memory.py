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

    print("-- concurrent peak: joblib runs n_cpu sets at once ------------------")
    print("  Each worker is a separate PROCESS holding its own database slice,")
    print("  so the peak is the sum over the n_cpu LARGEST sets running")
    print("  together -- not the largest single set, and not the total.")
    print()
    print(f"  {'--cpus':>7}  {'resident':>12}  {'with transient 2x':>18}"
          f"  {'suggested --mem':>16}")
    rows = []
    for c in args.cpus:
        top = [n for _, n in sizes[:c]]
        gb = sum(top) * n_motifs * BYTES_PER_RANK / 1024**3
        peak = gb * TRANSIENT_FACTOR
        mem = int(peak * 1.3) + 8
        rows.append((c, gb, peak, mem))
        print(f"  {c:>7}  {human(gb):>12}  {human(peak):>18}  {str(mem)+'G':>16}")
    print()
    print("  (transient 2x: subset_to_pandas builds a pyarrow Table and then a")
    print("   pandas DataFrame, both live for a moment. --mem adds 30% + 8 GB")
    print("   for the interpreter, motif annotations and the result objects.)")
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
    # BUG FIXED HERE: this was `next(...)`, which on an ascending --cpus list
    # returns the SMALLEST passing value while the text below calls it "the
    # largest that fits". On equal-N sets where every worker count fits, it
    # printed "--cpus-per-task 1" -- advising a 16x slowdown for no reason.
    fits = [c for c, _, peak, _ in rows if peak <= args.mem_limit_gb * 0.9]
    largest_ok = max(fits) if fits else None
    print("  Halving --cpus-per-task roughly halves the peak, at a proportional")
    print("  cost in wall time for THIS rule only -- the other rules keep")
    print(f"  whatever --cores you pass. There are {len(sizes)} region sets, so")
    print(f"  more than {len(sizes)} workers buys nothing.")
    if largest_ok:
        peak_at = next(pk for c, _, pk, _ in rows if c == largest_ok)
        print(f"  Under a {args.mem_limit_gb} GB allocation, the LARGEST worker "
              f"count that fits")
        print(f"  is --cpus-per-task {largest_ok} "
              f"(peak {human(peak_at)}, i.e. "
              f"{peak_at/args.mem_limit_gb:.0%} of the allocation).")
        if largest_ok == max(args.cpus):
            print(f"  That is the largest value tabulated -- there is headroom "
                  f"for more,")
            print(f"  though more than {len(sizes)} workers buys nothing "
                  f"(one per region set).")
    else:
        print("  No tabulated worker count fits in 128 GB; either raise --mem or")
        print("  cut the region sets down (01_cistopic/choose_dar_threshold.py")
        print("  --top-n gives equal-sized sets, which also makes this peak")
        print("  predictable rather than driven by the largest cell type).")
    print("=" * 74)


if __name__ == "__main__":
    main()
