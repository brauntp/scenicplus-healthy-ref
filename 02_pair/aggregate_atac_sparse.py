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
    ap.add_argument("--obs-tsv", type=Path, default=None,
                    help="TSV of per-cell metadata to merge into the ATAC obs, "
                         "keyed by cell id in the first column. This is how the "
                         "transferred labels reach the ATAC object when they live "
                         "in a sidecar file rather than in the .h5ad -- e.g. the "
                         "output of 02_pair/attach_atac_labels.py.")
    ap.add_argument("--obs-tsv-rna", type=Path, default=None,
                    help="same, for the RNA object")
    ap.add_argument("--rna-layer", default=None)
    ap.add_argument("--use-raw", action="store_true")
    ap.add_argument("--cells-per-metacell", type=int, default=50)
    ap.add_argument("--cells-per-metacell-atac", type=int, default=None)
    ap.add_argument("--min-cells-per-group", type=int, default=50)
    ap.add_argument("--oversample", type=float, default=2.0,
                    help="Metacells per group = (cells/k) * OVERSAMPLE. Values >1 "
                         "mean metacells share cells. This improves link RANKING "
                         "(more, better-conditioned observations for the GBM) but "
                         "does NOT add independent information -- see --help-stats. "
                         "4-8 is a reasonable range when compute is cheap; the "
                         "returns are largest for small populations. (default: 2)")
    ap.add_argument("--min-metacells-per-group", type=int, default=0,
                    help="Floor on metacells per group, raising the effective "
                         "oversample for SMALL populations only. This is where "
                         "oversampling earns its keep: a 300-cell population at "
                         "k=50 yields 6 non-overlapping metacells, too few for a "
                         "stable regression. 0 disables the floor.")
    ap.add_argument("--max-oversample", type=float, default=16.0,
                    help="Ceiling on the effective oversample the floor may reach. "
                         "Beyond ~16 the added metacells are near-duplicates.")
    ap.add_argument("--seed", type=int, default=666)
    ap.add_argument("--max-dense-gb", type=float, default=0.0,
                    help="Refuse if the metacell matrices would exceed this. Set "
                         "it to roughly 60%% of the job's --mem so the failure is "
                         "an immediate message rather than an OOM kill an hour "
                         "in. 0 disables the check.")
    ap.add_argument("--dry-run", action="store_true",
                    help="Plan the metacells, print the exact output shape and "
                         "memory, then stop without aggregating or writing. Run "
                         "this on a login node before submitting.")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--diagnostics", type=Path, default=None)
    args = ap.parse_args()

    # -- validate paths BEFORE importing or reading anything -----------------
    # An unset shell variable makes "$REF/rna.h5ad" expand to "/rna.h5ad", which
    # otherwise surfaces as an h5py traceback 20 frames deep. Catch it here and
    # name the likely cause.
    bad = []
    for flag, p in (("--rna", args.rna), ("--atac", args.atac),
                    ("--obs-tsv", args.obs_tsv), ("--obs-tsv-rna", args.obs_tsv_rna)):
        if p is None:
            continue
        if not p.exists():
            hint = ""
            # A path whose parent is the filesystem root is almost always an
            # empty variable: "$REF/atac.h5ad" with REF unset.
            if p.parent == Path("/"):
                hint = ("  <-- parent is '/', so an environment variable in this "
                        "path was EMPTY. Did you `export REF=...` in this shell? "
                        "Login shells do not inherit it from a previous session.")
            bad.append(f"  {flag} {p}{hint}")
    if bad:
        sys.exit("ERROR: input file(s) not found:\n" + "\n".join(bad) +
                 "\n\nNothing was read. Fix the paths and re-run; --dry-run is "
                 "free.")

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

    # -- sidecar metadata, merged BEFORE the group-key check ------------------
    for nm, ad_, tsv in (("RNA", rna, args.obs_tsv_rna),
                         ("ATAC", atac, args.obs_tsv)):
        if tsv is None:
            continue
        if not tsv.exists():
            sys.exit(f"ERROR: --obs-tsv{'-rna' if nm == 'RNA' else ''} "
                     f"{tsv} does not exist")
        _log(f"merging {nm} metadata from {tsv}")
        side = pd.read_csv(tsv, sep=None, engine="python", index_col=0,
                           dtype=str)
        side.index = side.index.astype(str)
        shared = ad_.obs_names.intersection(side.index)
        if len(shared) == 0:
            sys.exit(
                f"ERROR: no cell ids shared between {nm}.obs_names and {tsv}.\n"
                f"  {nm}.obs_names e.g. {list(ad_.obs_names[:3])}\n"
                f"  {tsv.name} index e.g. {list(side.index[:3])}\n"
                "These are different id spaces. Harmonise the ids -- do NOT "
                "drop the non-matching cells.")
        if len(shared) < ad_.n_obs:
            missing = ad_.n_obs - len(shared)
            sys.exit(
                f"ERROR: {missing:,} of {ad_.n_obs:,} {nm} cells are absent from "
                f"{tsv.name} ({len(shared):,} matched). Refusing to proceed on a "
                "partial merge -- silently dropping cells would bias every "
                "group. Regenerate the TSV for all cells, or subset the .h5ad "
                "deliberately.")
        new = [c for c in side.columns if c not in ad_.obs.columns]
        ad_.obs = ad_.obs.join(side.loc[ad_.obs_names, new])
        _log(f"  merged {len(new)} column(s) for all {ad_.n_obs:,} {nm} cells: "
             f"{new[:8]}{' ...' if len(new) > 8 else ''}")

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

    # Coerce labels to strings, but keep MISSING as missing. A bare
    # .astype(str) turns NaN into the literal string "nan", which then looks
    # like a perfectly good cell type: the unlabelled cells get pooled into a
    # junk group, metacelled, and handed to SCENIC+ as if they were a lineage.
    # Empty strings arrive the same way from a TSV with blank fields.
    MISSING = {"nan", "NaN", "NA", "None", "none", "null", "", "unassigned",
               "Unassigned", "unknown", "Unknown", "NULL", "<NA>"}

    def _labels(adata, nm):
        s = adata.obs[args.group_key]
        arr = s.astype(str).str.strip().to_numpy()
        bad = np.isin(arr, list(MISSING)) | s.isna().to_numpy()
        if bad.any():
            _log(f"{nm}: {bad.sum():,} of {len(arr):,} cells have no "
                 f"'{args.group_key}' label -- EXCLUDED "
                 f"(values seen: {sorted(set(arr[bad]))[:5]})")
            arr = arr.copy()
            arr[bad] = "\x00UNLABELLED"          # cannot collide with a real label
        return arr, int(bad.sum())

    gr, n_bad_r = _labels(rna, "RNA")
    ga, n_bad_a = _labels(atac, "ATAC")
    UNL = "\x00UNLABELLED"
    groups = sorted((set(gr) & set(ga)) - {UNL})
    only_r = sorted(set(gr) - set(ga) - {UNL})
    only_a = sorted(set(ga) - set(gr) - {UNL})
    if n_bad_r or n_bad_a:
        _log(f"unlabelled cells excluded from pairing: "
             f"RNA {n_bad_r:,}, ATAC {n_bad_a:,}")
    if only_r:
        _log(f"RNA-only groups (dropped): {only_r}")
    if only_a:
        _log(f"ATAC-only groups (dropped): {only_a}")

    # ---------------------------------------------------------------------
    # PHASE 1 -- plan only. Choose anchors and metacell memberships for every
    # group WITHOUT aggregating, so the exact output shape is known before any
    # dense memory is touched. Accumulating per-group blocks and np.vstack-ing
    # them at the end would hold two full copies of the dense output at once
    # (91 GB vs 51 GB on this reference); preallocating avoids that.
    # Memberships are tiny: 25k metacells x k x 8 bytes ~ 20 MB.
    # ---------------------------------------------------------------------
    plan, names, diags, dropped = [], [], {}, []
    for g in groups:
        ir = np.flatnonzero(gr == g)
        ia = np.flatnonzero(ga == g)
        if len(ir) < args.min_cells_per_group or len(ia) < args.min_cells_per_group:
            dropped.append({"group": g, "n_rna": int(len(ir)), "n_atac": int(len(ia))})
            continue
        n_cells_eff = min(len(ir), len(ia))
        base_mc = n_cells_eff / max(k_r, k_a)
        n_mc = max(1, int(base_mc * args.oversample))
        eff_ov = args.oversample
        if args.min_metacells_per_group and n_mc < args.min_metacells_per_group:
            capped = int(base_mc * args.max_oversample)
            n_mc = max(1, min(args.min_metacells_per_group, capped))
            eff_ov = n_mc / max(base_mc, 1e-9)
            _log(f"group {g!r}: small population, raising oversample to "
                 f"{eff_ov:.1f}x to reach {n_mc} metacells")
        Zr_g, Za_g = Zr[ir], Za[ia]
        anchors = maximin_anchors(np.vstack([Zr_g, Za_g]), n_mc, rng)
        A = np.vstack([Zr_g, Za_g])[anchors]
        mem_r = knn_rows(A, Zr_g, min(k_r, len(ir)))
        mem_a = knn_rows(A, Za_g, min(k_a, len(ia)))

        plan.append((g, ir, ia, mem_r, mem_a))
        names += [f"{g}_mc{i}" for i in range(A.shape[0])]

        d_r = 1.0 - (A * Zr_g[mem_r[:, 0]]).sum(1)
        d_a = 1.0 - (A * Za_g[mem_a[:, 0]]).sum(1)
        diags[g] = {"n_rna_cells": int(len(ir)), "n_atac_cells": int(len(ia)),
                    "n_metacells": int(A.shape[0]),
                    "effective_oversample": round(float(eff_ov), 2),
                    # Independent-observation equivalent: how many NON-overlapping
                    # metacells this population could support. Use THIS, not
                    # n_metacells, when judging how much a link is really supported.
                    "independent_metacell_equiv": int(max(1, base_mc)),
                    "mean_nn_dist_rna": float(np.mean(d_r)),
                    "mean_nn_dist_atac": float(np.mean(d_a)),
                    "median_crossmodal_gap": float(np.median(np.abs(d_r - d_a))),
                    "rna_cell_reuse_mean": float(mem_r.size / max(len(ir), 1)),
                    "atac_cell_reuse_mean": float(mem_a.size / max(len(ia), 1))}
        _log(f"group {g!r}: RNA {len(ir)}, ATAC {len(ia)} -> {A.shape[0]} metacells")

    if not plan:
        sys.exit("ERROR: every group was dropped; lower --min-cells-per-group")

    # ---------------------------------------------------------------------
    # PHASE 2 -- allocate once, fill in place.
    # ---------------------------------------------------------------------
    n_mc_total = sum(p[3].shape[0] for p in plan)
    gb = n_mc_total * (X_rna.shape[1] + X_atac.shape[1]) * 4 / 1024**3
    _log(f"allocating {n_mc_total:,} x {X_rna.shape[1]:,} (RNA) and "
         f"{n_mc_total:,} x {X_atac.shape[1]:,} (ATAC) float32 = {gb:.1f} GB")
    if args.max_dense_gb and gb > args.max_dense_gb:
        sys.exit(f"ERROR: the metacell matrices need {gb:.1f} GB, above "
                 f"--max-dense-gb {args.max_dense_gb:.1f}. Lower --oversample "
                 f"(halving it roughly halves this), or raise the limit if the "
                 f"job really has the memory.")
    if args.dry_run:
        print()
        print(f"DRY RUN -- nothing aggregated, nothing written.")
        print(f"  groups kept      : {len(plan)}")
        print(f"  groups dropped   : {len(dropped)}"
              + (f" {[d['group'] for d in dropped]}" if dropped else ""))
        print(f"  metacells        : {n_mc_total:,}")
        print(f"  dense output     : {gb:.1f} GB")
        print(f"  + sparse inputs  : "
              f"{(X_rna.data.nbytes + X_rna.indices.nbytes + X_atac.data.nbytes + X_atac.indices.nbytes)/1024**3:.1f} GB")
        print(f"  -> request --mem at least {int(gb*1.35 + 14)}G")
        print()
        print("  per group: name, rna_cells, atac_cells, metacells, indep_equiv")
        for g, ir, ia, mem_r, _ in plan:
            base = min(len(ir), len(ia)) / max(k_r, k_a)
            print(f"    {g:<26} {len(ir):>7,} {len(ia):>7,} "
                  f"{mem_r.shape[0]:>6,} {int(max(1, base)):>6}")
        return

    M_rna = np.empty((n_mc_total, X_rna.shape[1]), dtype=np.float32)
    M_atac = np.empty((n_mc_total, X_atac.shape[1]), dtype=np.float32)
    row = 0
    for g, ir, ia, mem_r, mem_a in plan:
        n = mem_r.shape[0]
        M_rna[row:row + n] = aggregate_sparse(X_rna[ir], mem_r)
        M_atac[row:row + n] = aggregate_sparse(X_atac[ia], mem_a)
        row += n
        _log(f"aggregated {g!r} ({row:,}/{n_mc_total:,} metacells)")
    assert row == n_mc_total, (row, n_mc_total)
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
        "oversample": args.oversample,
        "min_metacells_per_group": args.min_metacells_per_group,
        "caveat": ("Metacells pair cells by GLUE-latent proximity, not by "
                   "co-measurement. Accessibility is a mean of raw counts over k "
                   "cells, NOT topic-imputed: no cross-cell smoothing beyond the "
                   "averaging itself. Region-to-gene links are manifold-level "
                   "covariation, not single-nucleus co-occurrence."),
        "oversampling_caveat": (
            "With oversample > 1 metacells SHARE cells and are not independent "
            "observations. Simulation: oversampling improves true-vs-decoy link "
            "RANKING (AUROC 0.79 -> 0.98 at high integration noise, k=50) while "
            "median rho on true links stays flat -- but under a pure null the "
            "naive p<0.05 rate rises from 0.125 to 0.175 and p<0.001 from 0.000 "
            "to 0.050 as oversample goes 1 -> 16. Treat SCENIC+ region-to-gene "
            "p-values as ranking scores, not calibrated significance, and judge "
            "support using independent_metacell_equiv in the diagnostics rather "
            "than n_metacells.")}
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
