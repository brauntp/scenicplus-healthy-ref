#!/usr/bin/env python3
"""
Per-cell-type accessibility tracks (bedGraph -> bigWig) from the paired object.

WHAT THIS IS, AND WHAT IT IS NOT
--------------------------------
This writes a PEAK-LEVEL step function, not fragment coverage. Read that
sentence twice before interpreting a track.

A bigWig is base-pair signal, and there are two ways to produce one:

  1. From FRAGMENTS -- true per-bp coverage, shows footprints, shoulders,
     nucleosome structure. Requires per-sample `fragments.tsv.gz`.
  2. From a PEAK MATRIX -- constant within each peak, zero between peaks.

**This pipeline has no fragments.** The ATAC side of the reference is
`atac.h5ad`, a fixed-width peak x cell matrix on an inherited consensus peak
set; fragments were never an input. So option 1 is unavailable here, and this
script does option 2.

The consequence is real and must travel with the files: a track from this script
is FLAT ACROSS EACH PEAK and ZERO between peaks. It answers "how accessible is
this peak in this cell type" quantitatively and correctly. It cannot answer
anything sub-peak -- do not read a footprint, a summit position, or a peak
boundary off it, because the boundaries are the input peak set's and the shape
inside them is an artefact of the representation.

If sub-peak resolution is needed, that is a different job: get the fragments and
build coverage per cell-type barcode set.

NORMALISATION
-------------
Values are the mean of the paired object's `scATAC` layer over the metacells of
each group. That layer is already mean-per-metacell accessibility (max observed
1.14 on this reference, i.e. a rate rather than a count), so a group's track is
directly comparable to another group's without further scaling. No CPM or
library-size correction is applied, and none is appropriate: the quantity is
already an average of averages, not a count.

`--min-metacells` skips groups too small for a stable mean. The default (400)
matches the threshold established for the per-cell-type correlation analysis,
where route agreement fell from ~36% above it to ~1% below.

MEMORY
------
Streams the peak matrix in blocks. The full dense matrix is 25,323 x 393,832
float32 = 37.2 GB, so a naive `.X[:]` would OOM. A 20,000-peak block is 1.9 GB
resident; the default 8,000 is 0.76 GB, with the per-group sums held as
float64 vectors (393,832 x 13 x 8 B = 41 MB, negligible).

BIGWIG CONVERSION
-----------------
bedGraph is written unconditionally (it is a text format and needs nothing).
bigWig conversion needs ONE of:
  * `pyBigWig` (pip install pyBigWig) -- used automatically if importable;
  * UCSC `bedGraphToBigWig` on PATH -- used as a fallback.
If neither is present the bedGraphs are still written and the exact conversion
command is printed. A bedGraph loads in IGV and the UCSC browser directly; the
bigWig is only needed for remote hosting or very large tracks.

Usage:
    python 05_report/export_accessibility_tracks.py \\
        --h5mu ACC_GEX.h5mu \\
        --chromsizes anno/chromsizes.tsv \\
        --out-dir browser_tracks \\
        --group-key predicted_CellType_Broad --block 8000
"""
from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

try:
    import numpy as np
    import pandas as pd
except ImportError:
    sys.exit(
        "ERROR: this python lacks numpy/pandas.\n"
        f"       interpreter: {sys.executable}\n"
        "       conda activate scplus-pairing")


def parse_regions(names) -> pd.DataFrame:
    s = pd.Series(list(names), dtype=object)
    ex = s.str.extract(r"^(?P<chrom>[^:]+):(?P<start>\d+)-(?P<end>\d+)$")
    bad = ex.chrom.isna()
    if bad.any():
        raise SystemExit(
            f"ERROR: {int(bad.sum())} peak ids are not 'chr:start-end', "
            f"e.g. {s[bad].head(3).tolist()}")
    return pd.DataFrame({"chrom": ex.chrom.to_numpy(),
                         "start": ex.start.astype(np.int64).to_numpy(),
                         "end": ex.end.astype(np.int64).to_numpy()})


def to_bigwig(bg: Path, chromsizes: Path, bw: Path) -> str:
    """Returns the method used, or '' if neither backend is available."""
    try:
        import pyBigWig  # noqa: F401
    except ImportError:
        pass
    else:
        import pyBigWig
        cs = pd.read_csv(chromsizes, sep="\t")
        order = [(str(r.Chromosome), int(r.End)) for r in cs.itertuples()]
        d = pd.read_csv(bg, sep="\t", header=None,
                        names=["chrom", "start", "end", "value"])
        keep = set(c for c, _ in order)
        d = d[d.chrom.isin(keep)]
        # pyBigWig requires entries sorted by the header's chromosome order.
        rank = {c: i for i, (c, _) in enumerate(order)}
        d = d.assign(_r=d.chrom.map(rank)).sort_values(["_r", "start"])
        h = pyBigWig.open(str(bw), "w")
        h.addHeader(order)
        # pyBigWig's C extension wants native python int/float and rejects
        # numpy scalars with an unhelpful "valid set of entries" RuntimeError,
        # so coerce explicitly rather than relying on .tolist() dtypes.
        h.addEntries([str(x) for x in d.chrom],
                     [int(x) for x in d.start],
                     ends=[int(x) for x in d.end],
                     values=[float(x) for x in d.value])
        h.close()
        return "pyBigWig"

    exe = shutil.which("bedGraphToBigWig")
    if exe:
        cs2 = bg.with_suffix(".chromsizes.txt")
        cs = pd.read_csv(chromsizes, sep="\t")
        cs[["Chromosome", "End"]].to_csv(cs2, sep="\t", header=False, index=False)
        r = subprocess.run([exe, str(bg), str(cs2), str(bw)],
                           capture_output=True, text=True)
        cs2.unlink(missing_ok=True)
        if r.returncode == 0:
            return "bedGraphToBigWig"
        print(f"    bedGraphToBigWig failed: {r.stderr.strip()[:200]}")
    return ""


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--h5mu", type=Path, required=True)
    ap.add_argument("--chromsizes", type=Path, required=True)
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--group-key", default="predicted_CellType_Broad")
    ap.add_argument("--block", type=int, default=8000,
                    help="peaks per streaming block")
    ap.add_argument("--min-metacells", type=int, default=400,
                    help="skip groups below this; 400 matches the threshold "
                         "established for the per-cell-type correlations")
    ap.add_argument("--groups", nargs="*", default=None)
    ap.add_argument("--no-bigwig", action="store_true",
                    help="write bedGraph only")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    for p in (args.h5mu, args.chromsizes):
        if not p.exists():
            sys.exit(f"ERROR: not found: {p}")
    try:
        import mudata
    except ImportError:
        sys.exit("ERROR: this python lacks mudata.\n"
                 f"       interpreter: {sys.executable}\n"
                 "       conda activate scplus-pairing")

    print("=" * 74)
    print("per-cell-type accessibility tracks")
    print("=" * 74)

    md = mudata.read(str(args.h5mu), backed="r")
    A = md.mod["scATAC"]
    if args.group_key not in A.obs.columns:
        # labels may live on the RNA modality
        alt = md.mod.get("scRNA")
        if alt is not None and args.group_key in alt.obs.columns:
            grp = alt.obs[args.group_key].astype(str).to_numpy()
        else:
            sys.exit(f"ERROR: --group-key '{args.group_key}' not in scATAC.obs "
                     f"or scRNA.obs\n       scATAC has: {list(A.obs.columns)}")
    else:
        grp = A.obs[args.group_key].astype(str).to_numpy()

    n_mc, n_peak = A.shape
    co = parse_regions(A.var_names)
    sizes = pd.Series(grp).value_counts()
    groups = args.groups or sorted(sizes.index)
    keep = [g for g in groups if sizes.get(g, 0) >= args.min_metacells]
    skip = [g for g in groups if g not in keep]

    print(f"  object     : {args.h5mu.name}  [scATAC {n_mc:,} x {n_peak:,}]")
    print(f"  group key  : {args.group_key}")
    print(f"  groups     : {len(sizes)} present, {len(keep)} kept "
          f"(>= {args.min_metacells} metacells)")
    for g in skip:
        print(f"    skip '{g}': {int(sizes.get(g,0))} metacells")
    blk_gb = args.block * n_mc * 4 / 1024**3
    print(f"  block      : {args.block:,} peaks = {blk_gb:.2f} GB resident "
          f"(full matrix {n_mc*n_peak*4/1024**3:.1f} GB)")
    print(f"  blocks     : {int(np.ceil(n_peak/args.block))}")

    if not keep:
        sys.exit(
            f"\nERROR: no group has >= {args.min_metacells} metacells, so there is "
            "nothing to write.\n"
            f"       largest group: '{sizes.index[0]}' at {int(sizes.iloc[0]):,}.\n"
            "       Lower --min-metacells, but note the default (400) is the "
            "threshold\n"
            "       below which per-cell-type correlations stopped agreeing with an\n"
            "       independent route (~1% vs ~36%). A track from a smaller group "
            "is a\n"
            "       noisier mean, not a wrong one -- but treat it as indicative.")

    if args.dry_run:
        print()
        print("--dry-run: nothing written.")
        print("=" * 74)
        return

    # -- stream: accumulate per-group mean over peaks -------------------------
    masks = {g: (grp == g) for g in keep}
    counts = {g: int(m.sum()) for g, m in masks.items()}
    sums = {g: np.zeros(n_peak, dtype=np.float64) for g in keep}
    for i0 in range(0, n_peak, args.block):
        i1 = min(i0 + args.block, n_peak)
        X = A.X[:, i0:i1]
        X = np.asarray(X.todense() if hasattr(X, "todense") else X, dtype=np.float32)
        for g, m in masks.items():
            sums[g][i0:i1] = X[m].sum(0)
        if (i0 // args.block) % 5 == 0:
            print(f"    block {i0//args.block + 1}/{int(np.ceil(n_peak/args.block))} "
                  f"({min(i1, n_peak):,}/{n_peak:,} peaks)", flush=True)
        del X

    args.out_dir.mkdir(parents=True, exist_ok=True)
    made = []
    for g in keep:
        v = sums[g] / counts[g]
        safe = g.replace(" ", "_").replace("/", "-")
        d = (pd.DataFrame({"chrom": co.chrom, "start": co.start,
                           "end": co.end, "value": np.round(v, 5)})
             .query("value > 0")
             .sort_values(["chrom", "start"]))
        bg = args.out_dir / f"{safe}.accessibility.bedGraph"
        with open(bg, "w") as fh:
            fh.write(f'track type=bedGraph name="{g} accessibility" '
                     f'description="mean scATAC over {counts[g]:,} metacells; '
                     f'PEAK-LEVEL, flat within peaks" '
                     'visibility=full autoScale=on alwaysZero=on\n')
            d.to_csv(fh, sep="\t", header=False, index=False)
        note = ""
        if not args.no_bigwig:
            bw = args.out_dir / f"{safe}.accessibility.bw"
            how = to_bigwig(bg, args.chromsizes, bw)
            note = f" -> {bw.name} ({how})" if how else " (bedGraph only)"
        made.append((safe, len(d), counts[g]))
        print(f"  {safe:<26} {len(d):>8,} intervals, {counts[g]:>5,} metacells{note}")

    if not args.no_bigwig and not any(
            (args.out_dir / f"{s}.accessibility.bw").exists() for s, _, _ in made):
        print()
        print("  No bigWig backend found. bedGraphs are written and load directly")
        print("  in IGV and the UCSC browser. To convert:")
        print("      pip install pyBigWig      # then re-run")
        print("  or, with the UCSC tool:")
        print("      cut -f1,3 <chromsizes.tsv> | tail -n +2 > cs.txt")
        print("      bedGraphToBigWig <in>.bedGraph cs.txt <out>.bw")

    print()
    print("READ THIS BEFORE INTERPRETING A TRACK: these are PEAK-LEVEL step")
    print("functions, not fragment coverage. Flat within each peak, zero between.")
    print("This reference has no fragments -- the ATAC side is a peak x cell")
    print("matrix on an inherited consensus peak set. Quantitative comparison")
    print("between cell types at a peak is valid; anything sub-peak is not.")
    print()
    print(f"  wrote {args.out_dir}")
    print("=" * 74)


if __name__ == "__main__":
    main()
