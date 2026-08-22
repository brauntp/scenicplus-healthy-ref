#!/usr/bin/env python3
"""
Build a paired pseudo-multiome MuData for SCENIC+ from an UNPAIRED, GLUE-integrated
RNA + ATAC reference -- pairing cells through the GLUE latent space.

WHY THIS EXISTS
---------------
SCENIC+'s own non-multiome path (`process_non_multiome_data` ->
`generate_pseudocells_for_numpy`) builds each metacell by drawing cells uniformly
at random *within a categorical label*, independently for RNA and ATAC. The RNA
cells and ATAC cells inside a metacell therefore share nothing but the label.
Region-to-gene regression downstream then sees only between-label variance plus
sampling noise: all within-cell-type covariation -- exactly the graded structure a
GLUE co-embedding recovers -- is destroyed, and links collapse to cell-type
resolution.

This script instead anchors each metacell at a point in the GLUE latent space and
takes neighbours from both modalities around that anchor, so the RNA profile and
the ATAC profile in a metacell describe cells that are genuinely close in the
shared embedding.

WHAT IT PRESERVES
-----------------
Output is a plain 2-modality MuData with matching obs_names -- the same shape and
schema `prepare_GEX_ACC` emits. Drop it at the pipeline's
`combined_GEX_ACC_mudata` path and every downstream SCENIC+ rule is unchanged.

WHAT IT CANNOT DO
-----------------
Pairing is an INFERENCE, not a measurement. A region-to-gene link from
paired-by-embedding data is evidence that accessibility and expression covary
across the manifold, NOT that they co-occur in single nuclei. Anything resting on
within-cell-type resolution inherits the GLUE integration's error. Report the
diagnostics this script writes (neighbour distances, cross-modal gap, cell reuse)
alongside any eGRN derived from its output.

Usage
-----
    python glue_metacells.py \
        --rna RNA.h5ad --atac ATAC_imputed.h5ad \
        --latent-key X_glue --group-key celltype \
        --cells-per-metacell 25 \
        --out ACC_GEX.h5mu --diagnostics pairing_diagnostics
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd


def _log(msg: str) -> None:
    print(f"[glue_metacells] {msg}", flush=True)


def _densify(X, rows=None):
    sub = X[rows] if rows is not None else X
    if hasattr(sub, "toarray"):
        sub = sub.toarray()
    return np.asarray(sub, dtype=np.float32)


def l2_normalize(M: np.ndarray) -> np.ndarray:
    """GLUE latents are compared by cosine geometry; normalize before Euclidean kNN."""
    n = np.linalg.norm(M, axis=1, keepdims=True)
    n[n == 0] = 1.0
    return M / n


# ------------------------------------------------------------------ pairing
def build_metacells_for_group(
    Z_rna: np.ndarray,
    Z_atac: np.ndarray,
    n_metacells: int,
    k_rna: int,
    k_atac: int,
    rng: np.random.Generator,
    anchor_mode: str = "both",
):
    """
    Anchor metacells in the shared latent space, then gather neighbours from BOTH
    modalities around each anchor.

    Anchors come from farthest-point (maximin) sampling over the pooled latent
    coordinates, so metacells tile the manifold this cell type occupies instead of
    oversampling its dense core (which uniform random draws do).

    Returns
    -------
    rna_idx, atac_idx : list[np.ndarray]   local indices per metacell
    diag              : dict               per-group diagnostics
    """
    from sklearn.neighbors import NearestNeighbors

    n_r, n_a = len(Z_rna), len(Z_atac)
    k_rna = int(min(k_rna, n_r))
    k_atac = int(min(k_atac, n_a))

    if anchor_mode == "rna":
        pool = Z_rna
    elif anchor_mode == "atac":
        pool = Z_atac
    else:
        pool = np.vstack([Z_rna, Z_atac])

    n_anchor = int(min(n_metacells, len(pool)))
    first = int(rng.integers(len(pool)))
    chosen = [first]
    d2 = np.sum((pool - pool[first]) ** 2, axis=1)
    for _ in range(1, n_anchor):
        nxt = int(np.argmax(d2))
        if d2[nxt] <= 0:                      # no distinct points left
            remaining = np.setdiff1d(np.arange(len(pool)), np.asarray(chosen))
            if not len(remaining):
                break
            nxt = int(rng.choice(remaining))
        chosen.append(nxt)
        d2 = np.minimum(d2, np.sum((pool - pool[nxt]) ** 2, axis=1))
    anchors = pool[np.asarray(chosen)]

    nn_r = NearestNeighbors(n_neighbors=k_rna).fit(Z_rna)
    nn_a = NearestNeighbors(n_neighbors=k_atac).fit(Z_atac)
    dr, ir = nn_r.kneighbors(anchors)
    da, ia = nn_a.kneighbors(anchors)

    rna_idx = [ir[i] for i in range(len(anchors))]
    atac_idx = [ia[i] for i in range(len(anchors))]

    diag = {
        "n_metacells": int(len(anchors)),
        "k_rna": k_rna,
        "k_atac": k_atac,
        "mean_nn_dist_rna": float(np.mean(dr)),
        "mean_nn_dist_atac": float(np.mean(da)),
        "max_nn_dist_rna": float(np.max(dr)),
        "max_nn_dist_atac": float(np.max(da)),
        # Cross-modal gap: |nearest ATAC - nearest RNA| distance per anchor. Large
        # values mean the modalities do not co-occupy that region of the latent
        # space, i.e. the pairing there is extrapolation, not interpolation.
        "median_crossmodal_gap": float(np.median(np.abs(da[:, 0] - dr[:, 0]))),
        "p95_crossmodal_gap": float(np.percentile(np.abs(da[:, 0] - dr[:, 0]), 95)),
        "n_rna_cells_available": int(n_r),
        "n_atac_cells_available": int(n_a),
        "rna_cell_reuse_mean": float(len(anchors) * k_rna / max(n_r, 1)),
        "atac_cell_reuse_mean": float(len(anchors) * k_atac / max(n_a, 1)),
        "frac_rna_cells_used": float(len(np.unique(np.concatenate(rna_idx))) / max(n_r, 1)),
        "frac_atac_cells_used": float(len(np.unique(np.concatenate(atac_idx))) / max(n_a, 1)),
        "per_metacell_nn_dist_rna": dr.mean(axis=1).tolist(),
        "per_metacell_nn_dist_atac": da.mean(axis=1).tolist(),
    }
    return rna_idx, atac_idx, diag


def aggregate(X, idx_list, how="mean"):
    out = np.zeros((len(idx_list), X.shape[1]), dtype=np.float32)
    for i, idx in enumerate(idx_list):
        sub = _densify(X, np.sort(np.asarray(idx)))
        out[i] = sub.mean(axis=0) if how == "mean" else sub.sum(axis=0)
    return out


# ------------------------------------------------------------------ main
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--rna", type=Path, required=True)
    ap.add_argument("--atac", type=Path, required=True,
                    help="topic-imputed accessibility, cells x regions")
    ap.add_argument("--latent-key", default="X_glue",
                    help="obsm key holding the shared GLUE embedding")
    ap.add_argument("--latent-key-atac", default=None, help="override if ATAC key differs")
    ap.add_argument("--group-key", required=True,
                    help="obs column in BOTH files; pairing is confined within it")
    ap.add_argument("--rna-layer", default=None, help="layer to aggregate (default X)")
    ap.add_argument("--use-raw", action="store_true", help="aggregate RNA from .raw.X")
    ap.add_argument("--atac-layer", default=None)
    ap.add_argument("--cells-per-metacell", type=int, default=25)
    ap.add_argument("--cells-per-metacell-atac", type=int, default=None,
                    help="ATAC is sparser; a larger k is often warranted")
    ap.add_argument("--target-metacells-per-group", default="auto",
                    help="'auto' (min(n_rna,n_atac)/k * 2, SCENIC+'s convention) or an int")
    ap.add_argument("--min-cells-per-group", type=int, default=50,
                    help="groups with fewer cells in EITHER modality are dropped")
    ap.add_argument("--anchor-mode", choices=["rna", "atac", "both"], default="both")
    ap.add_argument("--aggregate", choices=["mean", "sum"], default="mean")
    ap.add_argument("--seed", type=int, default=666)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--diagnostics", type=Path, default=None,
                    help="stem for <stem>.json / <stem>.csv")
    args = ap.parse_args()

    import anndata
    import mudata

    rng = np.random.default_rng(args.seed)
    lk_rna = args.latent_key
    lk_atac = args.latent_key_atac or args.latent_key
    k_r = args.cells_per_metacell
    k_a = args.cells_per_metacell_atac or args.cells_per_metacell

    _log(f"reading {args.rna}")
    rna = anndata.read_h5ad(args.rna)
    _log(f"reading {args.atac}")
    atac = anndata.read_h5ad(args.atac)

    for nm, adata, lk in (("RNA", rna, lk_rna), ("ATAC", atac, lk_atac)):
        if lk not in adata.obsm:
            sys.exit(f"ERROR: latent key '{lk}' not in {nm}.obsm "
                     f"(available: {list(adata.obsm.keys())})")
        if args.group_key not in adata.obs.columns:
            sys.exit(f"ERROR: group key '{args.group_key}' not in {nm}.obs "
                     f"(available: {list(adata.obs.columns)})")

    Zr_all = l2_normalize(np.asarray(rna.obsm[lk_rna], dtype=np.float32))
    Za_all = l2_normalize(np.asarray(atac.obsm[lk_atac], dtype=np.float32))
    if Zr_all.shape[1] != Za_all.shape[1]:
        sys.exit(f"ERROR: latent dims differ (RNA {Zr_all.shape[1]}, "
                 f"ATAC {Za_all.shape[1]}). These are not a shared space.")
    if np.isnan(Zr_all).any() or np.isnan(Za_all).any():
        sys.exit("ERROR: NaN in latent coordinates.")

    if args.use_raw:
        if rna.raw is None:
            sys.exit("ERROR: --use-raw given but RNA has no .raw")
        X_rna, rna_var = rna.raw.X, rna.raw.var.copy()
        _log("aggregating RNA from .raw.X")
    elif args.rna_layer:
        X_rna, rna_var = rna.layers[args.rna_layer], rna.var.copy()
        _log(f"aggregating RNA from layer '{args.rna_layer}'")
    else:
        X_rna, rna_var = rna.X, rna.var.copy()
        _log("aggregating RNA from X")
    X_atac = atac.layers[args.atac_layer] if args.atac_layer else atac.X

    g_rna = rna.obs[args.group_key].astype(str).to_numpy()
    g_atac = atac.obs[args.group_key].astype(str).to_numpy()
    shared = sorted(set(g_rna) & set(g_atac))
    only_r = sorted(set(g_rna) - set(g_atac))
    only_a = sorted(set(g_atac) - set(g_rna))
    if only_r:
        _log(f"WARNING dropping RNA-only groups: {', '.join(only_r)}")
    if only_a:
        _log(f"WARNING dropping ATAC-only groups: {', '.join(only_a)}")
    if not shared:
        sys.exit(f"ERROR: no shared levels of '{args.group_key}' between modalities")

    all_rna, all_atac, names, group_of = [], [], [], []
    diags, dropped = {}, {}
    for g in shared:
        ri = np.flatnonzero(g_rna == g)
        ai = np.flatnonzero(g_atac == g)
        if len(ri) < args.min_cells_per_group or len(ai) < args.min_cells_per_group:
            dropped[g] = {"n_rna": int(len(ri)), "n_atac": int(len(ai)),
                          "reason": f"under --min-cells-per-group "
                                    f"({args.min_cells_per_group})"}
            _log(f"WARNING dropping group '{g}': RNA {len(ri)}, ATAC {len(ai)} cells")
            continue
        if args.target_metacells_per_group == "auto":
            n_mc = max(1, round(min(len(ri), len(ai)) / max(k_r, 1)) * 2)
        else:
            n_mc = int(args.target_metacells_per_group)

        _log(f"group '{g}': RNA {len(ri)}, ATAC {len(ai)} -> {n_mc} metacells")
        lr, la, d = build_metacells_for_group(
            Zr_all[ri], Za_all[ai], n_metacells=n_mc,
            k_rna=k_r, k_atac=k_a, rng=rng, anchor_mode=args.anchor_mode)
        lr = [ri[x] for x in lr]
        la = [ai[x] for x in la]
        all_rna.append(aggregate(X_rna, lr, args.aggregate))
        all_atac.append(aggregate(X_atac, la, args.aggregate))
        names += [f"{g}_{i}" for i in range(len(lr))]
        group_of += [g] * len(lr)
        diags[g] = d

    if not all_rna:
        sys.exit("ERROR: every group was dropped; nothing to write.")

    meta_rna = np.vstack(all_rna)
    meta_atac = np.vstack(all_atac)
    obs = pd.DataFrame({args.group_key: pd.Categorical(group_of)}, index=names)

    mdata = mudata.MuData({
        "scRNA": anndata.AnnData(X=meta_rna, obs=obs.copy(), var=rna_var),
        "scATAC": anndata.AnnData(X=meta_atac, obs=obs.copy(), var=atac.var.copy()),
    })
    mdata.uns["glue_pairing"] = {
        "method": "GLUE-latent maximin-anchored kNN metacells",
        "latent_key_rna": lk_rna, "latent_key_atac": lk_atac,
        "group_key": args.group_key,
        "cells_per_metacell_rna": k_r, "cells_per_metacell_atac": k_a,
        "anchor_mode": args.anchor_mode, "aggregate": args.aggregate,
        "seed": args.seed,
        "caveat": ("Metacells pair cells by proximity in the GLUE latent space, not "
                   "by co-measurement. Region-to-gene links are manifold-level "
                   "covariation, not single-nucleus co-occurrence."),
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    mdata.write_h5mu(str(args.out))
    _log(f"wrote {args.out}: {meta_rna.shape[0]} metacells, "
         f"{meta_rna.shape[1]} genes x {meta_atac.shape[1]} regions")

    if args.diagnostics:
        stem = args.diagnostics
        stem.parent.mkdir(parents=True, exist_ok=True)
        params = {k: (str(v) if isinstance(v, Path) else v) for k, v in vars(args).items()}
        payload = {"params": params,
                   "n_metacells_total": int(meta_rna.shape[0]),
                   "groups_used": list(diags.keys()),
                   "groups_dropped": dropped,
                   "groups_rna_only": only_r, "groups_atac_only": only_a,
                   "per_group": diags}
        with open(f"{stem}.json", "w") as fh:
            json.dump(payload, fh, indent=2, default=str)
        rows = [{"group": g, **{k: v for k, v in d.items()
                                if not k.startswith("per_metacell")}}
                for g, d in diags.items()]
        pd.DataFrame(rows).to_csv(f"{stem}.csv", index=False)
        _log(f"wrote diagnostics to {stem}.json / {stem}.csv")
        for g, d in sorted(diags.items(),
                           key=lambda kv: -kv[1]["median_crossmodal_gap"])[:3]:
            _log(f"largest cross-modal gap: '{g}' "
                 f"median={d['median_crossmodal_gap']:.4f} "
                 f"(RNA nn {d['mean_nn_dist_rna']:.4f}, "
                 f"ATAC nn {d['mean_nn_dist_atac']:.4f})")


if __name__ == "__main__":
    main()
