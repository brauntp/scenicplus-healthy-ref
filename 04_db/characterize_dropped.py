#!/usr/bin/env python3
"""
Are the peaks the cisTarget database drops real regulatory elements, or noise?

WHY THIS EXISTS
---------------
peak_overlap_audit.py says WHAT FRACTION of your peaks a database can score. It
cannot say whether losing them matters. An 18.8% loss is cheap if the dropped
peaks are weak, ubiquitous, technical peaks; it is expensive if they are the
cell-type-specific distal elements an eGRN analysis exists to find.

SCREEN is a cCRE catalogue, so peaks with ZERO overlap are by definition regions
ENCODE never nominated as candidate regulatory elements. That is not evidence
they are noise -- ENCODE's cCRE set is built from a fixed panel of biosamples and
misses elements specific to cell types outside it. Bone marrow progenitor
subsets are exactly the kind of population that can fall in that gap.

So this script asks your data, not the annotation:

  1. PREVALENCE  -- in how many cells is each peak accessible? Technical peaks
     are accessible in few cells; real elements in many.
  2. STRENGTH    -- mean accessibility. Same logic, continuous.
  3. SPECIFICITY -- max-group mean / mean over groups. A cell-type-restricted
     enhancer scores high; a housekeeping promoter or a noise peak scores low.

The verdict turns on 3: if dropped peaks are MORE cell-type-specific than kept
peaks, the database is discarding the informative tail and a custom build is
justified. If they are less specific AND lower prevalence, the loss is cheap.

MEMORY
------
Streams the sparse matrix in row blocks and never materialises it. Peak RSS is a
few GB on a 400k-peak object regardless of cell count -- login-node safe, though
a batch job is politer for the full reference.

Usage
-----
    python 04_db/characterize_dropped.py \\
        --atac        "$REF/atac.h5ad" \\
        --per-peak    audit_screen.per_peak.csv \\
        --obs-tsv     "$LABELS" \\
        --group-key   predicted_CellType_Broad \\
        --out         dropped_characterization
"""
import argparse
import json
import sys
from pathlib import Path

try:
    import numpy as np
    import pandas as pd
    import h5py
except ImportError as e:                                         # pragma: no cover
    sys.exit(f"ERROR: needs numpy, pandas and h5py ({e}).\n"
             f"       interpreter: {sys.executable}\n"
             "       conda activate scplus-pairing")

BLOCK = 2000


def _index_key(grp):
    k = grp.attrs.get("_index", "_index")
    return k.decode() if isinstance(k, bytes) else str(k)


def _strings(dset):
    return np.array([v.decode() if isinstance(v, bytes) else str(v)
                     for v in dset[:]])


def _read_categorical(grp, key):
    g = grp[key]
    if isinstance(g, h5py.Group) and "categories" in g:
        cats = _strings(g["categories"])
        codes = g["codes"][:]
        return np.where(codes >= 0, cats[np.clip(codes, 0, None)], "NA")
    return _strings(g)


def stream_stats(path, groups, n_var):
    """Per-peak sum, nonzero count, and per-group sum -- one pass, blocked."""
    with h5py.File(path, "r") as f:
        X = f["X"]
        if not isinstance(X, h5py.Group) or "indptr" not in X:
            sys.exit("ERROR: X is not CSR/CSC sparse; this script streams sparse "
                     "input only. A dense object needs a different approach.")
        enc = X.attrs.get("encoding-type", b"")
        enc = enc.decode() if isinstance(enc, bytes) else str(enc)
        if "csr" not in enc.lower():
            sys.exit(f"ERROR: X encoding is '{enc}'; expected csr_matrix "
                     "(cells x peaks). A csc object would need a column pass.")
        indptr = X["indptr"][:]
        n_obs = len(indptr) - 1
        uniq = np.unique(groups)
        gi = {g: i for i, g in enumerate(uniq)}
        gcode = np.array([gi[g] for g in groups], dtype=np.int32)

        tot = np.zeros(n_var, dtype=np.float64)
        nz = np.zeros(n_var, dtype=np.int64)
        gsum = np.zeros((len(uniq), n_var), dtype=np.float64)
        gn = np.bincount(gcode, minlength=len(uniq)).astype(np.float64)

        for i0 in range(0, n_obs, BLOCK):
            i1 = min(i0 + BLOCK, n_obs)
            lo, hi = indptr[i0], indptr[i1]
            if hi == lo:
                continue
            cols = X["indices"][lo:hi]
            vals = X["data"][lo:hi].astype(np.float64)
            rows = np.repeat(np.arange(i0, i1),
                            np.diff(indptr[i0:i1 + 1]))
            np.add.at(tot, cols, vals)
            np.add.at(nz, cols, 1)
            # per-group accumulation, one group at a time to keep it vectorised
            gr = gcode[rows]
            for gcode_i in np.unique(gr):
                m = gr == gcode_i
                np.add.at(gsum[gcode_i], cols[m], vals[m])
            if (i0 // BLOCK) % 10 == 0:
                print(f"  ... {i1:,}/{n_obs:,} cells", flush=True)
    return tot, nz, gsum, gn, uniq, n_obs


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--atac", type=Path, required=True)
    ap.add_argument("--per-peak", type=Path, required=True,
                    help="audit_screen.per_peak.csv from peak_overlap_audit.py")
    ap.add_argument("--obs-tsv", type=Path, default=None,
                    help="sidecar labels, if the object has no group column")
    ap.add_argument("--group-key", default="predicted_CellType_Broad")
    ap.add_argument("--out", type=Path, required=True, help="output prefix")
    ap.add_argument("--no-plot", action="store_true")
    args = ap.parse_args()

    for p in (args.atac, args.per_peak):
        if not p.exists():
            sys.exit(f"ERROR: missing input: {p}")

    # ---- representability flag, keyed by region id --------------------------
    pp = pd.read_csv(args.per_peak)
    cols = {c.lower(): c for c in pp.columns}
    flag_col = next((cols[k] for k in ("representable", "is_representable",
                                       "passes", "kept") if k in cols), None)
    if flag_col is None:
        sys.exit(f"ERROR: no representability column in {args.per_peak}. "
                 f"Columns: {list(pp.columns)}")
    id_col = next((cols[k] for k in ("region", "region_id", "peak", "peak_id",
                                     "name") if k in cols), None)
    if id_col is None:
        if {"chrom", "start", "end"} <= set(cols):
            pp["__id"] = (pp[cols["chrom"]].astype(str) + ":"
                          + pp[cols["start"]].astype(str) + "-"
                          + pp[cols["end"]].astype(str))
            id_col = "__id"
        else:
            sys.exit(f"ERROR: cannot key {args.per_peak} by region id. "
                     f"Columns: {list(pp.columns)}")
    flag = pp[flag_col]
    if flag.dtype == object:
        flag = flag.astype(str).str.lower().isin(("true", "yes", "1", "t"))
    rep_map = dict(zip(pp[id_col].astype(str), flag.astype(bool)))

    # ---- peaks and labels from the object ----------------------------------
    with h5py.File(args.atac, "r") as f:
        var = f["var"]
        peaks = _strings(var[_index_key(var)])
        obs = f["obs"]
        obs_names = _strings(obs[_index_key(obs)])
        groups = (_read_categorical(obs, args.group_key)
                  if args.group_key in obs else None)

    if groups is None:
        if args.obs_tsv is None:
            sys.exit(f"ERROR: '{args.group_key}' not in obs and no --obs-tsv "
                     "given. The labels live in the sidecar TSV for this "
                     "reference; pass it.")
        side = pd.read_csv(args.obs_tsv, sep="\t")
        if args.group_key not in side.columns:
            sys.exit(f"ERROR: '{args.group_key}' not in {args.obs_tsv} "
                     f"(columns: {list(side.columns)[:8]})")
        if len(side) != len(obs_names):
            sys.exit(f"ERROR: {args.obs_tsv} has {len(side):,} rows but the "
                     f"object has {len(obs_names):,} cells -- refusing to "
                     "assume positional alignment.")
        groups = side[args.group_key].astype(str).to_numpy()
        print(f"[labels] taken from {args.obs_tsv} by row order "
              f"({len(groups):,} rows match cell count)")

    rep = np.array([rep_map.get(p, False) for p in peaks])
    missing = sum(1 for p in peaks if p not in rep_map)
    if missing:
        print(f"WARNING: {missing:,} peaks absent from the audit CSV, counted "
              f"as dropped", file=sys.stderr)
    print(f"[peaks] {len(peaks):,} total | representable {rep.sum():,} "
          f"({rep.mean():.2%}) | dropped {(~rep).sum():,}")

    # ---- one streaming pass -------------------------------------------------
    print("[scan] streaming the sparse matrix in row blocks")
    tot, nz, gsum, gn, uniq, n_obs = stream_stats(args.atac, groups, len(peaks))

    prevalence = nz / n_obs
    strength = tot / n_obs
    gmean = gsum / gn[:, None]
    with np.errstate(divide="ignore", invalid="ignore"):
        overall = gmean.mean(0)
        specificity = np.where(overall > 0, gmean.max(0) / overall, np.nan)

    def summarize(mask, label):
        return dict(label=label, n=int(mask.sum()),
                    median_prevalence=float(np.median(prevalence[mask])),
                    median_strength=float(np.median(strength[mask])),
                    median_specificity=float(np.nanmedian(specificity[mask])))

    kept, drop = summarize(rep, "representable"), summarize(~rep, "dropped")
    from scipy.stats import mannwhitneyu
    res = {}
    for nm, arr in (("prevalence", prevalence), ("strength", strength),
                    ("specificity", specificity)):
        a, b = arr[rep], arr[~rep]
        ok_a, ok_b = a[np.isfinite(a)], b[np.isfinite(b)]
        u = mannwhitneyu(ok_b, ok_a, alternative="two-sided")
        # rank-biserial: >0 means dropped peaks score HIGHER
        rb = 2 * u.statistic / (len(ok_a) * len(ok_b)) - 1
        res[nm] = dict(p=float(u.pvalue), effect_dropped_higher=float(rb))

    print("\n" + "=" * 74)
    print("DROPPED-PEAK CHARACTERISATION")
    print("=" * 74)
    print(f"{'metric':<16}{'representable':>16}{'dropped':>12}"
          f"{'effect':>10}{'p':>12}")
    for nm, k, d in (("prevalence", kept["median_prevalence"],
                      drop["median_prevalence"]),
                     ("strength", kept["median_strength"],
                      drop["median_strength"]),
                     ("specificity", kept["median_specificity"],
                      drop["median_specificity"])):
        e = res[nm]["effect_dropped_higher"]
        print(f"{nm:<16}{k:>16.4f}{d:>12.4f}{e:>+10.3f}"
              f"{res[nm]['p']:>12.2e}")
    print("\neffect is rank-biserial: positive means DROPPED peaks score higher")

    spec_e = res["specificity"]["effect_dropped_higher"]
    prev_e = res["prevalence"]["effect_dropped_higher"]
    print("-" * 74)
    if spec_e > 0.1:
        verdict = "CUSTOM DATABASE JUSTIFIED"
        why = ("dropped peaks are MORE cell-type-specific than retained ones, "
               "so the database is discarding exactly the restricted elements "
               "an eGRN analysis is looking for")
    elif prev_e < -0.1 and spec_e < 0.05:
        verdict = "PRECOMPUTED DATABASE IS FINE"
        why = ("dropped peaks are both less prevalent and no more "
               "cell-type-specific -- consistent with weak or technical peaks, "
               "so the loss is cheap")
    else:
        verdict = "GENUINELY BORDERLINE"
        why = ("dropped peaks are not clearly worse than retained ones. The "
               "loss is unbiased but real; decide on cost, and report the "
               "dropped fraction either way")
    print(f"VERDICT: {verdict}")
    for line in (why[i:i + 70] for i in range(0, len(why), 70)):
        print(f"  {line}")
    print("-" * 74)

    out = {"n_peaks": int(len(peaks)), "n_cells": int(n_obs),
           "n_groups": int(len(uniq)),
           "representable": kept, "dropped": drop, "tests": res,
           "verdict": verdict, "rationale": why}
    Path(f"{args.out}.summary.json").write_text(json.dumps(out, indent=2))
    pd.DataFrame({"region": peaks, "representable": rep,
                  "prevalence": prevalence, "strength": strength,
                  "specificity": specificity}).to_csv(
        f"{args.out}.per_peak.csv", index=False)
    print(f"  summary : {args.out}.summary.json")
    print(f"  per-peak: {args.out}.per_peak.csv")

    if not args.no_plot:
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            plt.rcParams.update({"font.size": 7, "axes.titlesize": 8,
                                 "axes.labelsize": 7, "legend.fontsize": 6,
                                 "xtick.labelsize": 6, "ytick.labelsize": 6})
            fig, axes = plt.subplots(1, 3, figsize=(7.4, 2.6))
            for ax, (nm, arr, xl) in zip(axes, (
                    ("prevalence", prevalence, "fraction of cells accessible"),
                    ("strength", strength, "mean accessibility"),
                    ("specificity", specificity,
                     "max-group mean / mean over groups"))):
                a, b = arr[rep], arr[~rep]
                a, b = a[np.isfinite(a)], b[np.isfinite(b)]
                hi = np.nanpercentile(np.concatenate([a, b]), 99)
                bins = np.linspace(0, hi, 60)
                ax.hist(a, bins=bins, density=True, histtype="step", lw=1.4,
                        color="#3b7dd8", label="in database")
                ax.hist(b, bins=bins, density=True, histtype="step", lw=1.4,
                        color="#c44e52", label="dropped")
                ax.set_xlabel(xl)
                ax.set_title(nm, loc="left")
                for sp in ("top", "right"):
                    ax.spines[sp].set_visible(False)
            axes[0].set_ylabel("density")
            axes[0].legend(frameon=False)
            fig.tight_layout()
            fig.savefig(f"{args.out}.png", dpi=200)
            print(f"  figure  : {args.out}.png")
        except Exception as e:                                # pragma: no cover
            print(f"  (plot skipped: {e})")


if __name__ == "__main__":
    main()
