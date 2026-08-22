#!/usr/bin/env python3
"""
Aggregate raw sparse ATAC into metacells WITHOUT ever densifying per cell.

WHY THIS EXISTS
---------------
The canonical SCENIC+/pycisTopic route calls `impute_accessibility`, which
returns a DENSE `regions x cells` matrix. For this reference that is

    163,969 cells x 393,832 peaks x 4 bytes = 241 GB
    ... and ~481 GB with the transpose live during pairing.

No sensible partition has that. But the dense matrix is only ever consumed
after aggregation into metacells, and aggregation is a linear operation:

    metacell_profile = MEAN over the cells assigned to it

so it can be applied to the RAW SPARSE counts directly. The result is

    3,279 metacells (k=50) x 393,832 peaks x 4 bytes = 4.8 GB

a ~50x reduction, computed with a single sparse matrix product per cell type.
The averaging over 25-50 cells is itself the denoising that imputation was
providing; what is lost is the topic model's cross-cell smoothing, which is a
real difference and is stated in the pairing caveat rather than hidden.

Topic modelling (stage 4b) is still worth running for its BINARIZED TOPIC
REGION SETS, which pycistarget consumes. That path does not need the dense
imputed matrix -- only `impute_accessibility` does.

Usage
-----
    python aggregate_atac_sparse.py \
        --atac      atac.h5ad \
        --rna       rna.h5ad \
        --latent-key X_glue \
        --group-key  predicted_CellType_Broad \
        --cells-per-metacell 50 \
        --out       ACC_GEX.h5mu \
        --diagnostics pairing_diagnostics

This is a drop-in replacement for the `glue_metacells.py` ATAC input path when
imputation will not fit: it performs the SAME GLUE-anchored pairing, but reads
raw sparse counts instead of a dense imputed matrix.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from scipy import sparse


def _log(m):
    print(f"[aggregate_atac] {m}", flush=True)


def l2_normalize(Z):
    n = np.linalg.norm(Z, axis=1, keepdims=True)
    n[n == 0] = 1.0
    return Z / n


def maximin_anchors(Z, n_anchors, rng):
    """Farthest-point sampling: spread anchors over the manifold."""
    n = Z.shape[0]
    if n_anchors >= n:
        return np.arange(n)
    idx = [int(rng.integers(n))]
    d = np.linalg.norm(Z - Z[idx[0]], axis=1)
    for _ in range(1, n_anchors):
        nxt = int(np.argmax(d))
        idx.append(nxt)
        d = np.minimum(d, np.linalg.norm(Z - Z[nxt], axis=1))
    return np.asarray(idx)


def knn_rows(Z_query, Z_ref, k):
    """Indices of the k nearest Z_ref rows for each Z_query row, in blocks."""
    out = np.empty((Z_query.shape[0], k), dtype=np.int64)
    step = max(1, 2_000_000 // max(Z_ref.shape[0], 1))
    for s in range(0, Z_query.shape[0], step):
        e = min(s + step, Z_query.shape[0])
        d = Z_query[s:e] @ Z_ref.T                     # cosine, both L2-normed
        out[s:e] = np.argpartition(-d, kth=min(k, d.shape[1] - 1),
                                   axis=1)[:, :k]
    return out


def aggregate_sparse(X, member_idx):
    """
    Mean-aggregate rows of a sparse matrix into metacells.

    member_idx: (n_metacells, k) int array of row indices.
    Builds one (n_metacells x n_cells) averaging operator and applies it as a
    single sparse product -- never densifies the cell-level matrix.
    """
    n_mc, k = member_idx.shape
    rows = np.repeat(np.arange(n_mc), k)
    cols = member_idx.ravel()
    vals = np.full(rows.shape[0], 1.0 / k, dtype=np.float32)
    A = sparse.csr_matrix((vals, (rows, cols)), shape=(n_mc, X.shape[0]),
                          dtype=np.float32)
    return np.asarray((A @ X).todense(), dtype=np.float32)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--atac", type=Path, required=True, help="RAW sparse peak matrix")
    ap.add_argument("--rna", type=Path, required=True)
    ap.add_argument("--latent-key", default="X_glue")
    ap.add_argument("--group-key", required=True)
    ap.add_argument("--rna-layer", default=None)
    ap.add_argument("--use-raw", action="store_true")
    ap.add_argument("--cells-per-metacell", type=int, default=50)
    ap.add_argument("--cells-per-metacell-atac", type=int, default=None)
    ap.add_argument("--min-cells-per-group", type=int, default=50)
    ap.add_argument("--seed", type=int, default=666)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--diagnostics", type=Path, default=None)
    args = ap.parse_args()

    import anndata
    import mudata
    import pandas as pd

    rng = np.random.default_rng(args.seed)
    k_r = args.cells_per_metacell
    k_a = args.cells_per_metacell_atac or k_r

    _log(f"reading {args.rna} (backed)")
    rna = anndata.read_h5ad(args.rna)
    _log(f"reading {args.atac} (backed)")
    atac = anndata.read_h5ad(args.atac)

    for nm, ad_, in (("RNA", rna), ("ATAC", atac)):
        if args.latent_key not in ad_.obsm:
            sys.exit(f"ERROR: '{args.latent_key}' not in {nm}.obsm "
                     f"(have: {list(ad_.obsm.keys())})")
        if args.group_key not in ad_.obs.columns:
            sys.exit(f"ERROR: '{args.group_key}' not in {nm}.obs")

    Zr = l2_normalize(np.asarray(rna.obsm[args.latent_key], dtype=np.float32))
    Za = l2_normalize(np.asarray(atac.obsm[args.latent_key], dtype=np.float32))
    if Zr.shape[1] != Za.shape[1]:
        sys.exit(f"ERROR: latent dims differ ({Zr.shape[1]} vs {Za.shape[1]})")

    X_rna = (rna.raw.X if args.use_raw else
             rna.layers[args.rna_layer] if args.rna_layer else rna.X)
    rna_var = rna.raw.var.copy() if args.use_raw else rna.var.copy()
    X_atac = atac.X
    if not sparse.issparse(X_atac):
        _log("WARNING: ATAC X is dense; this script exists to avoid that")
        X_atac = sparse.csr_matrix(X_atac)
    X_rna = sparse.csr_matrix(X_rna) if not sparse.issparse(X_rna) else X_rna.tocsr()
    X_atac = X_atac.tocsr()

    gr = rna.obs[args.group_key].astype(str).to_numpy()
    ga = atac.obs[args.group_key].astype(str).to_numpy()
    groups = sorted(set(gr) & set(ga))
    only_r = sorted(set(gr) - set(ga))
    only_a = sorted(set(ga) - set(gr))
    if only_r:
        _log(f"RNA-only groups (dropped): {only_r}")
    if only_a:
        _log(f"ATAC-only groups (dropped): {only_a}")

    blocks_r, blocks_a, names, diags, dropped = [], [], [], {}, []
    for g in groups:
        ir = np.flatnonzero(gr == g)
        ia = np.flatnonzero(ga == g)
        if len(ir) < args.min_cells_per_group or len(ia) < args.min_cells_per_group:
            dropped.append({"group": g, "n_rna": int(len(ir)), "n_atac": int(len(ia))})
            continue
        n_mc = max(1, int(min(len(ir), len(ia)) / max(k_r, k_a) * 2))
        Zr_g, Za_g = Zr[ir], Za[ia]
        anchors = maximin_anchors(np.vstack([Zr_g, Za_g]), n_mc, rng)
        A = np.vstack([Zr_g, Za_g])[anchors]
        mem_r = knn_rows(A, Zr_g, min(k_r, len(ir)))
        mem_a = knn_rows(A, Za_g, min(k_a, len(ia)))

        blocks_r.append(aggregate_sparse(X_rna[ir], mem_r))
        blocks_a.append(aggregate_sparse(X_atac[ia], mem_a))
        names += [f"{g}_mc{i}" for i in range(A.shape[0])]

        d_r = 1.0 - (A * Zr_g[mem_r[:, 0]]).sum(1)
        d_a = 1.0 - (A * Za_g[mem_a[:, 0]]).sum(1)
        diags[g] = {"n_rna_cells": int(len(ir)), "n_atac_cells": int(len(ia)),
                    "n_metacells": int(A.shape[0]),
                    "mean_nn_dist_rna": float(np.mean(d_r)),
                    "mean_nn_dist_atac": float(np.mean(d_a)),
                    "median_crossmodal_gap": float(np.median(np.abs(d_r - d_a))),
                    "rna_cell_reuse_mean": float(mem_r.size / max(len(ir), 1)),
                    "atac_cell_reuse_mean": float(mem_a.size / max(len(ia), 1))}
        _log(f"group {g!r}: RNA {len(ir)}, ATAC {len(ia)} -> {A.shape[0]} metacells")

    if not blocks_r:
        sys.exit("ERROR: every group was dropped; lower --min-cells-per-group")

    M_rna = np.vstack(blocks_r)
    M_atac = np.vstack(blocks_a)
    _log(f"metacell matrices: RNA {M_rna.shape} "
         f"({M_rna.nbytes/1024**3:.2f} GB), ATAC {M_atac.shape} "
         f"({M_atac.nbytes/1024**3:.2f} GB)")

    obs = pd.DataFrame({"metacell": names,
                        args.group_key: [n.rsplit("_mc", 1)[0] for n in names]}
                       ).set_index("metacell")
    md = mudata.MuData({
        "scRNA": anndata.AnnData(X=M_rna, obs=obs.copy(), var=rna_var),
        "scATAC": anndata.AnnData(X=M_atac, obs=obs.copy(), var=atac.var.copy())})
    md.uns["glue_pairing"] = {
        "method": "GLUE-latent maximin-anchored kNN metacells, RAW SPARSE aggregation",
        "imputation": "NONE -- aggregation over k cells replaces topic imputation",
        "latent_key": args.latent_key, "group_key": args.group_key,
        "cells_per_metacell_rna": k_r, "cells_per_metacell_atac": k_a,
        "seed": args.seed,
        "caveat": ("Metacells pair cells by GLUE-latent proximity, not by "
                   "co-measurement. Accessibility is a mean of raw counts over k "
                   "cells, NOT topic-imputed: no cross-cell smoothing beyond the "
                   "averaging itself. Region-to-gene links are manifold-level "
                   "covariation, not single-nucleus co-occurrence.")}
    args.out.parent.mkdir(parents=True, exist_ok=True)
    md.write_h5mu(str(args.out))
    _log(f"wrote {args.out}: {M_rna.shape[0]} metacells, "
         f"{M_rna.shape[1]} genes x {M_atac.shape[1]} regions")

    if args.diagnostics:
        stem = str(args.diagnostics)
        payload = {"params": {k: (str(v) if isinstance(v, Path) else v)
                              for k, v in vars(args).items()},
                   "n_metacells_total": int(M_rna.shape[0]),
                   "groups_used": list(diags.keys()),
                   "groups_dropped": dropped,
                   "groups_rna_only": only_r, "groups_atac_only": only_a,
                   "per_group": diags}
        with open(f"{stem}.json", "w") as fh:
            json.dump(payload, fh, indent=2, default=str)
        pd.DataFrame([{"group": g, **d} for g, d in diags.items()]).to_csv(
            f"{stem}.csv", index=False)
        _log(f"wrote {stem}.json / {stem}.csv")
        for g, d in sorted(diags.items(),
                           key=lambda kv: -kv[1]["median_crossmodal_gap"])[:3]:
            _log(f"largest cross-modal gap: {g!r} "
                 f"median={d['median_crossmodal_gap']:.4f}")


if __name__ == "__main__":
    main()
