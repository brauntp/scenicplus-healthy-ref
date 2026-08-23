#!/usr/bin/env python3
"""
Measure -- do not project -- the memory the region-to-gene step will need.

WHY THIS EXISTS
---------------
Two projections in this project were wrong in opposite directions:

  pairing job : predicted ~55 GB, actual 63.5 GB   (15% UNDER)
  QC job      : predicted ~80.7 GB, actual 43.2 GB (1.87x OVER)

The QC overshoot came from assuming a scan materialises a full second copy; it
streams in chunks instead. The pairing undershoot came from inferring sparse
input size from a compressed file size. Both were arithmetic on an assumed
access pattern, and the assumption was the error.

So this does not estimate. It reads the real `.h5mu`, calls the same
`.to_df()` that scenicplus calls, and reports peak RSS -- first on a column
subset to prove the shape of the growth, then optionally on the full object.

WHAT SCENIC+ ACTUALLY DOES
--------------------------
`scenicplus.enhancer_to_gene.calculate_regions_to_genes_relationships` takes
DataFrames, and the CLI path reaches it after reading the mudata and calling
`.to_df()` on each modality.

MEASURED, not assumed: on a dense float32 `.h5mu`, `.to_df()` returns a frame
that SHARES memory with the AnnData's X (`np.shares_memory(...)` is True) and
adds +0.00 GB resident. There is no second copy. The earlier "to_df() doubles
the footprint" projection -- which produced the ~96G figure -- was wrong in the
same way the QC estimate was: assuming an access pattern rather than measuring
one.

Two reasons that does NOT settle the --mem line, and why this script exists
rather than a corrected arithmetic:

  1. The no-copy behaviour depends on the pandas/anndata versions in the
     SCENIC+ environment, which is not the environment the observation above was
     made in. A version that consolidates blocks would copy.
  2. `region_to_gene` then forks joblib workers that memory-map per-worker
     slices from `temp_dir`. The parent's peak is a FLOOR, not a ceiling.

Run it with --full in the SCENIC+ env, in a generously sized job, and set
--mem from what it reports.

Usage
-----
    # cheap: subsample columns, extrapolate the slope (seconds, <5 GB)
    python 03_pipeline/probe_region_to_gene_memory.py --h5mu ACC_GEX.h5mu

    # authoritative: the real thing, in a job sized generously
    python 03_pipeline/probe_region_to_gene_memory.py --h5mu ACC_GEX.h5mu --full

Report a `--full` run's peak RSS against the `#SBATCH --mem` in
slurm/scenicplus.sbatch before trusting that number.
"""
import argparse
import gc
import json
import resource
import sys
import time
from pathlib import Path

try:
    import numpy as np
    import h5py
except ImportError as e:                                         # pragma: no cover
    sys.exit(f"ERROR: needs numpy and h5py ({e})")


def peak_rss_gb():
    """Peak RSS in GB. ru_maxrss is KiB on Linux, bytes on macOS."""
    kb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return kb / 1024**2 if sys.platform != "darwin" else kb / 1024**3


def current_rss_gb():
    """Resident set now, in GB. Linux exposes this via /proc; elsewhere fall
    back to the high-water mark, which is an upper bound rather than 'now'."""
    try:
        with open("/proc/self/statm") as f:
            pages = int(f.read().split()[1])
        return pages * 4096 / 1024**3
    except Exception:
        return peak_rss_gb()          # not 'current', but never NaN


def shapes(path, mods=("scRNA", "scATAC")):
    out = {}
    with h5py.File(path, "r") as f:
        for m in mods:
            if m not in f.get("mod", {}):
                continue
            X = f["mod"][m]["X"]
            if isinstance(X, h5py.Dataset):
                out[m] = dict(shape=tuple(X.shape), dtype=str(X.dtype),
                              dense=True)
            else:
                enc = X.attrs.get("encoding-type", b"")
                enc = enc.decode() if isinstance(enc, bytes) else str(enc)
                sh = X.attrs.get("shape")
                out[m] = dict(shape=tuple(int(v) for v in sh) if sh is not None
                              else None, dtype=str(X["data"].dtype),
                              dense=False, encoding=enc,
                              nnz=int(X["data"].shape[0]))
    return out


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--h5mu", type=Path, required=True)
    ap.add_argument("--full", action="store_true",
                    help="read the whole object and call .to_df() for real")
    ap.add_argument("--cols", type=int, default=20_000,
                    help="ATAC columns to load in subset mode (default 20000)")
    ap.add_argument("--out", type=Path, default=None)
    args = ap.parse_args()

    if not args.h5mu.exists():
        sys.exit(f"ERROR: {args.h5mu} not found")

    sh = shapes(args.h5mu)
    print("=" * 74)
    print("region-to-gene memory probe")
    print("=" * 74)
    for m, d in sh.items():
        kind = "dense" if d["dense"] else f"sparse({d.get('encoding')})"
        print(f"  {m:<8} {d['shape']} {d['dtype']} {kind}"
              + ("" if d["dense"] else f" nnz={d['nnz']:,}"))
    rna, atac = sh.get("scRNA"), sh.get("scATAC")
    if not (rna and atac):
        sys.exit("ERROR: need both scRNA and scATAC modalities")

    n_mc = atac["shape"][0]
    dense_gb = {m: np.prod(d["shape"]) * 4 / 1024**3 for m, d in sh.items()}
    print(f"\n  on-disk dense equivalent: "
          + ", ".join(f"{m} {g:.1f} GB" for m, g in dense_gb.items())
          + f"  (total {sum(dense_gb.values()):.1f} GB)")

    try:
        import mudata
        import pandas as pd
    except ImportError as e:
        sys.exit(f"ERROR: needs mudata and pandas for the probe ({e})\n"
                 "       conda activate scplus-pairing (or the scenicplus env)")

    res = dict(h5mu=str(args.h5mu), shapes=sh,
               dense_equivalent_gb={k: float(v) for k, v in dense_gb.items()})

    if not args.full:
        n = min(args.cols, atac["shape"][1])
        print(f"\n[subset] loading scRNA whole + {n:,} of "
              f"{atac['shape'][1]:,} ATAC columns, then .to_df() on both")
        t0 = time.time()
        md = mudata.read(str(args.h5mu), backed=True)
        sub = md["scATAC"][:, :n].to_memory()
        rna_mem = md["scRNA"].to_memory()
        after_read = current_rss_gb()
        d_atac = sub.to_df()
        d_rna = rna_mem.to_df()
        peak = peak_rss_gb()
        el = time.time() - t0
        frac = n / atac["shape"][1]
        loaded = dense_gb["scRNA"] + dense_gb["scATAC"] * frac
        print(f"  loaded dense equivalent : {loaded:.2f} GB")
        print(f"  RSS after read          : {after_read:.2f} GB")
        print(f"  peak RSS incl. to_df()  : {peak:.2f} GB")
        print(f"  overhead factor         : {peak / loaded:.2f}x")
        print(f"  elapsed                 : {el:.1f} s")
        full_pred = sum(dense_gb.values()) * (peak / loaded)
        print(f"\n  EXTRAPOLATED full-object peak: {full_pred:.1f} GB "
              f"(dense {sum(dense_gb.values()):.2f} GB x "
              f"{peak / loaded:.2f})")
        print(f"  suggested --mem              : {full_pred * 1.25:.0f}G "
              f"(25% headroom)")
        print("\n  NOTE: this is an extrapolation from a column subset. The "
              "overhead\n  factor is what it measures; run --full in a job to "
              "confirm it holds\n  at scale, since that is exactly the step "
              "where my projections failed.")
        res.update(mode="subset", n_cols=int(n),
                   loaded_gb=float(loaded), peak_rss_gb=float(peak),
                   overhead_factor=float(peak / loaded),
                   extrapolated_full_gb=float(full_pred))
        del d_atac, d_rna, sub, rna_mem, md
        gc.collect()
    else:
        print(f"\n[full] reading the whole object unbacked and calling "
              f".to_df() on both modalities -- this is the real thing")
        t0 = time.time()
        md = mudata.read(str(args.h5mu))
        after_read = current_rss_gb()
        print(f"  RSS after mudata.read() : {after_read:.2f} GB")
        d_rna = md["scRNA"].to_df()
        after_rna = current_rss_gb()
        print(f"  RSS after scRNA.to_df() : {after_rna:.2f} GB "
              f"(+{after_rna - after_read:.2f})")
        d_atac = md["scATAC"].to_df()
        after_atac = current_rss_gb()
        peak = peak_rss_gb()
        el = time.time() - t0
        print(f"  RSS after scATAC.to_df(): {after_atac:.2f} GB "
              f"(+{after_atac - after_rna:.2f})")
        print(f"  PEAK RSS                : {peak:.2f} GB")
        print(f"  dense equivalent        : {sum(dense_gb.values()):.2f} GB")
        print(f"  overhead factor         : "
              f"{peak / sum(dense_gb.values()):.2f}x")
        print(f"  elapsed                 : {el:.1f} s")
        print(f"\n  suggested --mem for the SCENIC+ job: "
              f"{max(1.0, peak * 1.25):.0f}G  (measured peak + 25%)")
        print("  region_to_gene then forks joblib workers that memory-map "
              "slices\n  from temp_dir, so the parent's peak is the floor, not "
              "the ceiling.")
        res.update(mode="full", peak_rss_gb=float(peak),
                   rss_after_read_gb=float(after_read),
                   overhead_factor=float(peak / sum(dense_gb.values())),
                   suggested_mem_gb=float(peak * 1.25))
        del d_rna, d_atac, md
        gc.collect()

    print("=" * 74)
    if args.out:
        args.out.write_text(json.dumps(res, indent=2))
        print(f"  wrote {args.out}")


if __name__ == "__main__":
    main()
