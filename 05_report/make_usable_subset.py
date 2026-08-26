#!/usr/bin/env python3
"""
Restrict the per-cell-type correlation table to the USABLE cell-type groups.

WHY THIS EXISTS
---------------
`celltype_rho.py` scores every group with >= --min-metacells (default 30) and
writes `celltype_rho.tsv` -- 23 groups on this reference. But only 13 of those
groups carry enough metacells for the attribution to agree with an independent
route: measured against the regulon-membership route, per-group agreement was
~36% for groups with >= 400 metacells and ~1% below (Spearman rho = +0.89
between group size and agreement).

`best_group` / `specificity_z` / `best_rho` in the 23-group table are computed
by competing ALL scored groups, so an underpowered group can win a link on
noise. This script recomputes those three columns over the usable groups only,
which is the table every downstream deliverable (browser tracks, the AML
handoff, the notebook entry) actually uses.

THIS STEP WAS MISSING FROM THE REPO. The 13-group table was originally derived
in an interactive session and published without a script behind it, so
`slurm/browser_tracks.sbatch` defaulted to a path that had never existed on the
cluster and stage 3 failed with "not found:
05_report/celltype_rho_usable.csv.gz". A published table with no generator is
the same defect class as a figure with no generator: it cannot be reproduced,
re-derived after a fix, or checked.

WHAT IS AND IS NOT RECOMPUTED
-----------------------------
Recomputed over the usable groups: `best_group`, `best_abs_z`, `second_abs_z`,
`specificity_z`, `best_rho`, `n_groups_scored`.

Carried through unchanged: `region`, `target`, `Distance`, `importance_global`,
`rho_global`, and each kept group's own `rho__<group>` / `z__<group>` columns --
a within-group correlation does not depend on which other groups exist.

Dropped: the `rho__` / `z__` columns of excluded groups. They are not deleted
data -- `celltype_rho.tsv` keeps them, and the group-size threshold analysis was
derived FROM the full table.

Usage:
    python 05_report/make_usable_subset.py \\
        --in 05_report/celltype_rho.tsv \\
        --floors 05_report/celltype_rho.detection_floors.csv \\
        --out 05_report/celltype_rho_usable.csv.gz \\
        --min-metacells 400
"""
from __future__ import annotations

import argparse
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


def unwrap_distance(s: pd.Series) -> pd.Series:
    """Upstream ships Distance as '[-126154]'; pd.to_numeric on that returns
    ALL-NaN silently. celltype_rho.py already unwraps, but accept both so this
    script is safe against an older table."""
    if s.dtype.kind in "iuf":
        return s
    return (s.astype(str).str.strip("[]").replace("", np.nan)
            .astype("Float64").astype("Int64"))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--in", dest="inp", type=Path, required=True,
                    help="celltype_rho.tsv (all scored groups)")
    ap.add_argument("--floors", type=Path, required=True,
                    help="celltype_rho.detection_floors.csv, written beside it")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--min-metacells", type=int, default=400,
                    help="usability threshold; 400 is where route agreement "
                         "fell from ~36%% to ~1%%")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    for p in (args.inp, args.floors):
        if not p.exists():
            sys.exit(
                f"ERROR: not found: {p}\n"
                "       Run the correlation job first:\n"
                "         sbatch --partition=batch slurm/celltype_rho.sbatch\n"
                "       It writes celltype_rho.tsv and "
                "celltype_rho.detection_floors.csv.")

    print("=" * 74)
    print("usable-group subset")
    print("=" * 74)

    fl = pd.read_csv(args.floors)
    if not {"group", "metacells"} <= set(fl.columns):
        sys.exit(f"ERROR: {args.floors} lacks 'group'/'metacells'")
    keep = sorted(fl.loc[fl.metacells >= args.min_metacells, "group"])
    drop = sorted(fl.loc[fl.metacells < args.min_metacells, "group"])
    if not keep:
        sys.exit(f"ERROR: no group has >= {args.min_metacells} metacells.\n"
                 f"       largest: {fl.metacells.max():,}")

    # Detect the separator by PARSING, not by extension. `Path.suffix` on
    # 'celltype_rho.tsv.gz' is '.gz', so an extension rule silently chose ','
    # for a TAB file: the header came back as one column, no rho__ was found,
    # and the group-mismatch guard below then blamed the floors file. A
    # confidently wrong diagnosis is worse than no guard at all.
    head, sep = None, None
    for cand in ("\t", ","):
        h = pd.read_csv(args.inp, sep=cand, nrows=0).columns.tolist()
        if any(c.startswith("rho__") for c in h):
            head, sep = h, cand
            break
    if head is None:
        sys.exit(f"ERROR: no 'rho__<group>' column found in {args.inp.name} "
                 "under either tab or comma separation.\n"
                 "       Is this really the output of celltype_rho.py?")

    rho_all = [c for c in head if c.startswith("rho__")]
    scored = [c.removeprefix("rho__") for c in rho_all]
    missing = [g for g in keep if g not in scored]
    if missing:
        sys.exit(f"ERROR: group(s) pass the threshold but are absent from "
                 f"{args.inp.name}: {missing}\n"
                 f"       table has {len(scored)} scored group(s): {scored}\n"
                 f"       The floors file and the table must come from the "
                 f"SAME run.")

    carry = [c for c in ("region", "target", "Distance",
                         "importance_global", "rho_global") if c in head]
    rho_keep = [f"rho__{g}" for g in keep]
    z_keep = [f"z__{g}" for g in keep if f"z__{g}" in head]

    print(f"  input      : {args.inp.name}  ({len(scored)} scored groups)")
    print(f"  threshold  : >= {args.min_metacells} metacells")
    print(f"  keep       : {len(keep)}  {keep}")
    print(f"  drop       : {len(drop)}  {drop}")
    print(f"  columns    : {len(carry)} carried + {len(rho_keep)} rho + "
          f"{len(z_keep)} z + 6 recomputed")

    if args.dry_run:
        print()
        print("--dry-run: nothing written.")
        print("=" * 74)
        return

    df = pd.read_csv(args.inp, sep=sep, usecols=carry + rho_keep + z_keep)
    if "Distance" in df.columns:
        df["Distance"] = unwrap_distance(df["Distance"])
    print(f"  read       : {len(df):,} pairs")

    # Recompute the contrast over the kept groups only. Fisher z/SE, NOT raw
    # rho: raw |rho| is not comparable across group sizes (95th percentile under
    # pure noise was 0.054 at n=1200 but 0.202 at n=90), so a raw contrast
    # favours the smallest group.
    n_by_group = fl.set_index("group").metacells
    if z_keep and len(z_keep) == len(rho_keep):
        Z = df[[f"z__{g}" for g in keep]].to_numpy(float)
    else:
        # Older table without z__ columns: derive them.
        M = df[rho_keep].to_numpy(float)
        nvec = np.array([float(n_by_group[g]) for g in keep])
        Z = np.arctanh(np.clip(M, -0.999999, 0.999999)) * np.sqrt(
            np.maximum(nvec - 3.0, 1.0))
        for g, z in zip(keep, Z.T):
            df[f"z__{g}"] = z
        print("  note       : z__ columns absent from the input; derived them")

    M = df[rho_keep].to_numpy(float)
    filled = np.where(np.isnan(Z), -np.inf, np.abs(Z))
    srt = np.sort(filled, axis=1)
    bi = np.argmax(filled, axis=1)
    # Rows where NO kept group has a defined rho carry no evidence at all, yet
    # np.argmax still returns index 0 -- which puts a real cell-type name in
    # best_group as though a call had been made. The published 23-group table
    # has 5,355 such rows and they were the ONLY source of best_group
    # disagreement when this script was checked against it (every numeric column
    # matched to 0.0). Mark them instead of letting column order decide.
    n_def = np.isfinite(M).sum(1)
    df["best_group"] = np.where(n_def > 0, [keep[i] for i in bi], "none")
    df["best_abs_z"] = np.where(np.isfinite(srt[:, -1]), srt[:, -1], np.nan)
    df["second_abs_z"] = np.where(np.isfinite(srt[:, -2]), srt[:, -2], np.nan)
    df["specificity_z"] = df.best_abs_z - df.second_abs_z
    df["best_rho"] = M[np.arange(len(M)), bi]
    df["n_groups_scored"] = np.isfinite(M).sum(1)

    order = (carry + ["best_group", "best_rho", "specificity_z",
                      "n_groups_scored", "best_abs_z", "second_abs_z"]
             + rho_keep + [f"z__{g}" for g in keep])
    df = df[[c for c in order if c in df.columns]]
    df = df.sort_values("specificity_z", ascending=False)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.out, index=False)

    n_pass = int(((df.specificity_z >= 3) & (df.best_rho.abs() >= 0.15)).sum())
    print(f"  wrote      : {args.out}  ({len(df):,} pairs x {len(keep)} groups)")
    print()
    print(f"  links with specificity_z >= 3 and |best_rho| >= 0.15: {n_pass:,}")
    print(f"  Distance dtype: {df.Distance.dtype if 'Distance' in df else 'absent'}")
    print()
    print("  RANK ON specificity_z, not on raw rho differences. An empty")
    print("  rho__<group> cell means zero variance in that group (a peak or gene")
    print("  constant there), NOT zero correlation.")
    print("=" * 74)


if __name__ == "__main__":
    main()
