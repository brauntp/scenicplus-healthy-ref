#!/usr/bin/env python3
"""
Extract a small plotting bundle on the cluster, so nothing large is transferred.

WHY THIS EXISTS
---------------
The AUCell objects (183 MB + 87 MB) already carry eRegulon activity per metacell,
and the metacell NAMES encode the cell type (`B_mc0`, `HSC MPP_mc12`), so they are
self-sufficient for per-cell-type plots -- no grouping table needed.

The one thing they lack is TF expression, which lives in the paired object's
scRNA matrix inside a 41 GB file. Extracting just the TF columns turns a 41 GB
transfer into ~22 MB:

    full paired object        40.4 GB
    its scRNA part alone       3.30 GB
    just the ~229 TF columns  22.1 MB

MEMORY: reads the paired object with h5py in row blocks and keeps only the TF
columns, so peak is one block (default 4,000 metacells x 34,112 genes float32 =
0.5 GB), not the matrix. Safe on a login node, unlike anything that calls
mudata.read() on this file.

Usage
-----
    python 05_report/extract_for_plots.py            # defaults are correct here
    python 05_report/extract_for_plots.py --out-dir /somewhere/else
"""

from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

try:
    import h5py
    import numpy as np
    import pandas as pd
except ImportError as exc:
    sys.exit(f"ERROR: missing {exc.name}.\n"
             f"       interpreter: {sys.executable}\n"
             "       conda activate scplus-pairing")


def decode(arr) -> np.ndarray:
    """HDF5 string datasets come back as bytes; anndata may also use categories."""
    out = np.asarray(arr)
    if out.dtype.kind == "S":
        return np.array([b.decode() for b in out])
    if out.dtype.kind == "O":
        return np.array([x.decode() if isinstance(x, bytes) else str(x) for x in out])
    return out.astype(str)


def read_index(h5: h5py.File, path: str) -> np.ndarray:
    """anndata stores the index name in an attr, and the values under it."""
    grp = h5[path]
    key = grp.attrs.get("_index", "_index")
    if isinstance(key, bytes):
        key = key.decode()
    if key not in grp:
        raise KeyError(f"{path}: no index dataset (attr _index={key!r}); "
                       f"keys present: {list(grp)[:8]}")
    return decode(grp[key][:])


def tf_list(direct: Path, extended: Path) -> list[str]:
    """TFs with at least one eRegulon. Read from the tables, not hardcoded."""
    tfs: set[str] = set()
    for p in (direct, extended):
        if p.exists():
            tfs |= set(pd.read_table(p, usecols=["TF"])["TF"].unique())
        else:
            print(f"  note: {p} absent, skipping")
    return sorted(tfs)


def extract_tf_expression(h5mu: Path, tfs: list[str], block: int) -> pd.DataFrame:
    with h5py.File(h5mu, "r") as h5:
        if "mod/scRNA" not in h5:
            raise KeyError(f"{h5mu}: no mod/scRNA (modalities: "
                           f"{list(h5.get('mod', []))})")
        g = h5["mod/scRNA"]
        genes = read_index(g, "var")
        cells = read_index(g, "obs")
        X = g["X"]
        if isinstance(X, h5py.Group):
            raise TypeError(f"{h5mu}: scRNA/X is sparse; this extractor assumes "
                            "the dense layout aggregate_atac_sparse.py writes.")

        pos = {gn: i for i, gn in enumerate(genes)}
        found = [t for t in tfs if t in pos]
        missing = [t for t in tfs if t not in pos]
        cols = np.array([pos[t] for t in found])
        order = np.argsort(cols)          # h5py needs increasing fancy indices
        cols_sorted = cols[order]
        names_sorted = [found[i] for i in order]

        print(f"  TFs requested {len(tfs)}, present in scRNA {len(found)}")
        if missing:
            print(f"  absent from the expression matrix ({len(missing)}): "
                  f"{', '.join(missing[:8])}{' ...' if len(missing) > 8 else ''}")
            print("    (a TF can have an eRegulon from motif evidence while its")
            print("     own transcript was filtered out of the RNA matrix)")

        out = np.empty((X.shape[0], len(cols_sorted)), dtype=np.float32)
        for i0 in range(0, X.shape[0], block):
            i1 = min(i0 + block, X.shape[0])
            out[i0:i1] = X[i0:i1, :][:, cols_sorted]
        df = pd.DataFrame(out, index=pd.Index(cells, name="metacell"),
                          columns=names_sorted)
    # The group is in the metacell name -- see 02_pair/aggregate_atac_sparse.py:344
    df.insert(0, "group", [n.rsplit("_mc", 1)[0] for n in df.index])
    return df


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--paired", type=Path, default=Path("ACC_GEX.h5mu"))
    ap.add_argument("--auc-direct", type=Path,
                    default=Path("03_pipeline/AUCell_direct.h5mu"))
    ap.add_argument("--auc-extended", type=Path,
                    default=Path("03_pipeline/AUCell_extended.h5mu"))
    ap.add_argument("--ereg-direct", type=Path,
                    default=Path("03_pipeline/eRegulon_direct.tsv"))
    ap.add_argument("--ereg-extended", type=Path,
                    default=Path("03_pipeline/eRegulons_extended.tsv"))
    ap.add_argument("--out-dir", type=Path, default=Path("05_report/plot_bundle"))
    ap.add_argument("--block", type=int, default=4000,
                    help="metacell rows per read; peak memory is one block")
    args = ap.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    print("=" * 74)
    print("plotting bundle")
    print("=" * 74)

    print("-- TF list from the eRegulon tables --------------------------------")
    tfs = tf_list(args.ereg_direct, args.ereg_extended)
    print(f"  {len(tfs)} TFs with at least one eRegulon")
    print()

    print("-- TF expression from the paired object (blocked read) -------------")
    if not args.paired.exists():
        sys.exit(f"ERROR: paired object not found: {args.paired}")
    tf_df = extract_tf_expression(args.paired, tfs, args.block)
    tf_out = args.out_dir / "tf_expression.parquet"
    tf_df.to_parquet(tf_out)
    print(f"  wrote {tf_out}  ({tf_out.stat().st_size / 1024**2:.1f} MB, "
          f"{tf_df.shape[0]:,} x {tf_df.shape[1] - 1})")
    print()

    print("-- copying the small files as-is -----------------------------------")
    for src in (args.auc_direct, args.auc_extended,
                args.ereg_direct, args.ereg_extended):
        if not src.exists():
            print(f"  MISSING {src}")
            continue
        dst = args.out_dir / src.name
        shutil.copy2(src, dst)
        print(f"  {dst.name:<26} {dst.stat().st_size / 1024**2:>8.1f} MB")

    total = sum(p.stat().st_size for p in args.out_dir.iterdir() if p.is_file())
    print()
    print(f"  bundle: {args.out_dir}  ({total / 1024**2:.0f} MB total)")
    print(f"  vs transferring the paired object: "
          f"{args.paired.stat().st_size / 1024**3:.0f} GB")
    print("=" * 74)


if __name__ == "__main__":
    main()
