#!/usr/bin/env python3
"""
Structural check on a .h5mu WITHOUT reading its matrices -- h5py only.

WHY THIS EXISTS
---------------
`validate_h5mu.py` calls `mudata.read()`, which materialises every matrix. On
the 41 GB paired object that is a multi-minute, tens-of-GB read: fine inside the
batch job that produced it, wrong on a login node. This reads only HDF5
metadata -- dataset shapes, dtypes, attributes, and a few scattered chunks --
so it answers "did the write complete and is the object well formed" in under a
second at any file size.

It cannot replace validate_h5mu.py: a full NaN scan needs every value. Use this
first, and run the full validator inside a batch job (slurm/qc_paired.sbatch).

Usage
-----
    python 03_pipeline/peek_h5mu.py ACC_GEX.h5mu
    python 03_pipeline/peek_h5mu.py ACC_GEX.h5mu --expect-metacells 25323 \\
        --expect-groups 24 --group-key predicted_CellType_Broad

Exit 0 if the structure is sound, 2 if something is definitely wrong.
"""
import argparse
import sys
from pathlib import Path

try:
    import h5py
except ImportError:
    sys.exit("ERROR: needs h5py.\n"
             f"       interpreter: {sys.executable}\n"
             "       bash 00_inspect/find_inspect_env.sh finds one that has it")
import numpy as np


def _strings(dset, limit=None):
    n = dset.shape[0] if limit is None else min(limit, dset.shape[0])
    out = dset[:n]
    return [v.decode() if isinstance(v, bytes) else str(v) for v in out]


def _index_key(grp):
    k = grp.attrs.get("_index", "_index")
    return k.decode() if isinstance(k, bytes) else str(k)


def _mat_info(g):
    """Shape/dtype of an .X that may be dense or sparse, without reading it."""
    if isinstance(g, h5py.Dataset):
        return {"kind": "dense", "shape": tuple(g.shape), "dtype": str(g.dtype),
                "nbytes_dense": int(np.prod(g.shape)) * g.dtype.itemsize}
    enc = g.attrs.get("encoding-type", b"")
    enc = enc.decode() if isinstance(enc, bytes) else str(enc)
    shape = g.attrs.get("shape")
    d = g["data"]
    return {"kind": enc or "sparse", "shape": tuple(int(x) for x in shape) if shape is not None else None,
            "dtype": str(d.dtype), "nnz": int(d.shape[0]),
            "indptr_len": int(g["indptr"].shape[0]) if "indptr" in g else None}


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("h5mu", type=Path)
    ap.add_argument("--expect-metacells", type=int, default=None)
    ap.add_argument("--expect-groups", type=int, default=None)
    ap.add_argument("--group-key", default="predicted_CellType_Broad")
    args = ap.parse_args()

    if not args.h5mu.exists():
        sys.exit(f"ERROR: {args.h5mu} does not exist")
    size_gb = args.h5mu.stat().st_size / 1024**3
    fail = []

    print(f"{args.h5mu}  ({size_gb:.1f} GB)")
    try:
        f = h5py.File(args.h5mu, "r")
    except OSError as e:
        # The single most useful thing this catches: a file whose write was
        # killed mid-flush is not a valid HDF5 container at all.
        sys.exit(f"ERROR: not a readable HDF5 file -- the write likely did not "
                 f"complete.\n       {e}")
    with f:
        mods = sorted(f["mod"].keys()) if "mod" in f else []
        print(f"modalities: {mods}")
        if not {"scRNA", "scATAC"} <= set(mods):
            fail.append(f"expected scRNA and scATAC, found {mods}")

        obs_names = {}
        for m in mods:
            g = f["mod"][m]
            xi = _mat_info(g["X"])
            ok_key = _index_key(g["obs"])
            names = _strings(g["obs"][ok_key], limit=3)
            n_obs = int(g["obs"][ok_key].shape[0])
            vk = _index_key(g["var"])
            n_var = int(g["var"][vk].shape[0])
            obs_names[m] = (n_obs, names)
            dense_gb = (xi.get("nbytes_dense") or
                        (np.prod(xi["shape"]) * 4 if xi.get("shape") else 0)) / 1024**3
            print(f"\n[{m}]")
            print(f"  X: {xi['kind']} {xi.get('shape')} {xi['dtype']}"
                  + (f"  nnz={xi['nnz']:,}" if "nnz" in xi else "")
                  + f"  (~{dense_gb:.1f} GB dense)")
            print(f"  obs: {n_obs:,} x {len(g['obs'].keys())} cols   var: {n_var:,}")
            print(f"  obs_names[:3]: {names}")
            if xi.get("shape") and xi["shape"][0] != n_obs:
                fail.append(f"{m}: X rows {xi['shape'][0]} != obs {n_obs}")
            if xi.get("shape") and xi["shape"][1] != n_var:
                fail.append(f"{m}: X cols {xi['shape'][1]} != var {n_var}")
            # Spot-check real values at three positions rather than scanning.
            X = g["X"]
            try:
                if isinstance(X, h5py.Dataset):
                    rows = sorted({0, X.shape[0] // 2, X.shape[0] - 1})
                    samp = np.concatenate([np.asarray(X[r, :min(2000, X.shape[1])])
                                           for r in rows])
                else:
                    d = X["data"]
                    take = min(200_000, d.shape[0])
                    samp = np.asarray(d[:take])
                finite = np.isfinite(samp)
                print(f"  spot-check {samp.size:,} values: "
                      f"nonfinite={int((~finite).sum())}, "
                      f"min={float(samp[finite].min()) if finite.any() else 'n/a'}, "
                      f"max={float(samp[finite].max()) if finite.any() else 'n/a'}")
                if not finite.all():
                    fail.append(f"{m}: NaN/inf found in the sampled values")
            except Exception as e:
                print(f"  spot-check failed: {type(e).__name__}: {e}")
                fail.append(f"{m}: could not read sample values ({type(e).__name__})")

            if m == "scRNA" and args.group_key in g["obs"]:
                og = g["obs"][args.group_key]
                if isinstance(og, h5py.Group) and "categories" in og:
                    cats = _strings(og["categories"])
                    codes = np.asarray(og["codes"])
                    used = np.unique(codes[codes >= 0])
                    print(f"  {args.group_key}: {len(used)} groups in use "
                          f"({len(cats)} categories)")
                    junk = [cats[i] for i in used
                            if str(cats[i]).strip().lower() in
                            {"nan", "", "na", "none", "unassigned", "unknown"}]
                    if junk:
                        fail.append(f"junk group(s) present: {junk}")
                    if args.expect_groups and len(used) != args.expect_groups:
                        fail.append(f"{len(used)} groups, expected {args.expect_groups}")
                    counts = {cats[i]: int((codes == i).sum()) for i in used}
                    top = sorted(counts.items(), key=lambda kv: -kv[1])
                    print("  largest: " + ", ".join(f"{k}={v:,}" for k, v in top[:4]))
                    print("  smallest: " + ", ".join(f"{k}={v:,}" for k, v in top[-3:]))

        if {"scRNA", "scATAC"} <= set(obs_names):
            nr, sr = obs_names["scRNA"]
            na, sa = obs_names["scATAC"]
            print(f"\nobs counts match: {nr == na} ({nr:,} vs {na:,})")
            print(f"first obs_names match: {sr == sa}")
            if nr != na:
                fail.append(f"metacell counts differ: RNA {nr} vs ATAC {na}")
            if sr != sa:
                fail.append("obs_names differ between modalities -- not paired")
            if args.expect_metacells and nr != args.expect_metacells:
                fail.append(f"{nr:,} metacells, expected {args.expect_metacells:,}")

    print()
    if fail:
        print("STRUCTURE PROBLEMS:")
        for x in fail:
            print(f"  - {x}")
        sys.exit(2)
    print("STRUCTURE OK (metadata only -- a full NaN scan needs "
          "03_pipeline/validate_h5mu.py in a batch job)")


if __name__ == "__main__":
    main()
