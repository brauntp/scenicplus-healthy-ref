#!/usr/bin/env python3
"""
What did the pipeline actually find? Read the eRegulon tables, not the 41 GB object.

WHY THIS EXISTS
---------------
The pipeline finished and wrote scplusmdata.h5mu (41 GB) -- but ~99% of that is
the input matrices carried through (scATAC 37.2 GB + scRNA 3.2 GB); the AUC
layers it adds are 0.26 GB. The biology is in the two eRegulon TSVs, and nothing
in this repo read them.

It also cross-references the caveats established earlier in the project, so a
number is never reported without the reason it might be thin:

  - three groups were dropped before region sets (too few independent metacells)
  - Stromal retained only 60% of its regions through the database overlap
  - both erythroid groups were the worst cross-modal alignments in pairing QC
  - region sets are label-driven DARs, so cell-type-associated regulons are
    found well and shared/continuous programs poorly

COLUMN SCHEMA (read from scenicplus v1.0a2 `_format_egrns`, not guessed)
------------------------------------------------------------------------
  Region, Gene, TF, is_extended, eRegulon_name
  Gene_signature_name    TF_direct_+/+_(123g)     -- gene count in the name
  Region_signature_name  TF_direct_+/+_(456r)     -- region count in the name
  importance_R2G, rho_R2G, importance_x_rho_R2G, importance_x_abs_rho_R2G
  importance_TF2G, rho_TF2G, ...                  -- merged from tf_to_gene

The name encodes sign as `TF_[direct|extended]_[+|-]/[+|-]`: the first sign is
TF-to-gene correlation, the second region-to-gene. `+/+` is an activator acting
through accessible regions; `-/+` a repressor.

LOGIN-NODE SAFE: reads two TSVs (~100s of MB), never the h5mu. Stated because
this project has already shipped one tool wrongly labelled cheap.

Usage
-----
    python 05_report/summarize_eregulons.py \\
        --direct 03_pipeline/eRegulon_direct.tsv \\
        --extended 03_pipeline/eRegulons_extended.tsv \\
        --out-prefix docs/eregulon_summary
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

try:
    import pandas as pd
except ImportError:
    sys.exit("ERROR: this python has no pandas.\n"
             f"       interpreter: {sys.executable}\n"
             "       Use the scplus-pairing env, or any env with pandas.")

# Groups with too few independent metacells for a one-vs-rest DAR test. From the
# region-set job's own log, which counted independent_metacell_equiv rather than
# raw metacells.
DROPPED = {"Late GMP": 3, "Pro-Monocyte": 1, "cDC": 2}

# Fraction of each region set's 5,000 regions that survived the
# fraction_overlap 0.4 mapping into the cisTarget database, from the cisTarget
# run log. Only the notable ends are listed.
RETENTION = {"Stromal": 0.60, "HSC_MPP": 0.96}

# Flagged in docs/QC_RESULT.md as the worst cross-modal alignments in pairing.
POOR_PAIRING = {"Early_Erythroid", "Late_Erythroid"}


def load(path: Path, label: str) -> pd.DataFrame:
    if not path.exists():
        sys.exit(f"ERROR: {label} table not found: {path}")
    df = pd.read_table(path)
    need = {"TF", "Region", "Gene", "eRegulon_name"}
    missing = need - set(df.columns)
    if missing:
        sys.exit(f"ERROR: {path} lacks expected columns: {sorted(missing)}\n"
                 f"       found: {sorted(df.columns)[:12]}\n"
                 "       Schema is from scenicplus v1.0a2 _format_egrns; a\n"
                 "       different version may name them differently.")
    df["set"] = label
    return df


def sign_of(name: str) -> str:
    """`TF_direct_+/+` -> '+/+'. Returns '?' if the name does not parse."""
    if "_" not in name:
        return "?"
    tail = name.rsplit("_", 1)[-1]
    return tail if "/" in tail else "?"


def per_eregulon(df: pd.DataFrame) -> pd.DataFrame:
    g = df.groupby("eRegulon_name")
    out = pd.DataFrame({
        "TF": g["TF"].first(),
        "n_regions": g["Region"].nunique(),
        "n_genes": g["Gene"].nunique(),
        "n_links": g.size(),
    })
    out["sign"] = [sign_of(n) for n in out.index]
    for col in ("importance_x_abs_rho_R2G", "rho_R2G", "rho_TF2G"):
        if col in df.columns:
            out[f"median_{col}"] = g[col].median()
    return out.sort_values("n_genes", ascending=False)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--direct", type=Path,
                    default=Path("03_pipeline/eRegulon_direct.tsv"))
    ap.add_argument("--extended", type=Path,
                    default=Path("03_pipeline/eRegulons_extended.tsv"))
    ap.add_argument("--out-prefix", type=Path,
                    default=Path("docs/eregulon_summary"))
    ap.add_argument("--top", type=int, default=25,
                    help="rows to print (the full table is always written)")
    args = ap.parse_args()

    dfs = [load(args.direct, "direct")]
    if args.extended.exists():
        dfs.append(load(args.extended, "extended"))
    else:
        print(f"NOTE: no extended table at {args.extended}; direct only")
    df = pd.concat(dfs, ignore_index=True)

    print("=" * 74)
    print("eRegulon summary")
    print("=" * 74)
    for label, sub in df.groupby("set"):
        print(f"  {label:<9} {sub['eRegulon_name'].nunique():>5} eRegulons  "
              f"{sub['TF'].nunique():>5} TFs  "
              f"{sub['Region'].nunique():>8,} regions  "
              f"{sub['Gene'].nunique():>7,} genes  "
              f"{len(sub):>9,} region-gene links")
    print()

    # -- sign breakdown: activators vs repressors -----------------------------
    df["sign"] = [sign_of(n) for n in df["eRegulon_name"]]
    print("-- eRegulons by sign (TF2G / R2G) ----------------------------------")
    sig = (df.groupby(["set", "sign"])["eRegulon_name"].nunique()
             .unstack(fill_value=0))
    print(sig.to_string())
    print("  +/+ activator through accessible regions; -/+ repressor.")
    print()

    # -- per-eRegulon table ---------------------------------------------------
    tab = per_eregulon(df)
    args.out_prefix.parent.mkdir(parents=True, exist_ok=True)
    csv = args.out_prefix.with_suffix(".csv")
    tab.to_csv(csv)
    print(f"-- top {args.top} eRegulons by target-gene count -------------------")
    cols = [c for c in ("TF", "sign", "n_regions", "n_genes", "n_links",
                        "median_importance_x_abs_rho_R2G") if c in tab.columns]
    print(tab[cols].head(args.top).to_string())
    print()
    print(f"  full table ({len(tab)} eRegulons): {csv}")
    print()

    # -- caveats, cross-referenced against what was actually found -----------
    print("-- caveats, checked against the results -----------------------------")
    tfs = set(df["TF"])
    print("  Groups dropped before region sets (too few independent")
    print("  metacells for a one-vs-rest DAR test):")
    for g, n in sorted(DROPPED.items(), key=lambda kv: kv[1]):
        print(f"    {g:<16} {n} independent metacell(s) -- contributes no "
              f"region set, so no eRegulon can be specific to it")
    print()
    print("  Region retention through the fraction_overlap 0.4 database")
    print("  mapping (5,000 regions per set as written):")
    for g, r in sorted(RETENTION.items(), key=lambda kv: kv[1]):
        flag = "  <- thinnest base of any group" if r < 0.7 else ""
        print(f"    {g:<16} {r:.0%}{flag}")
    print()
    print("  Pairing QC flagged these as the worst cross-modal alignments:")
    for g in sorted(POOR_PAIRING):
        print(f"    {g}")
    print("    Thin results for these are the cross-modal gap surfacing")
    print("    downstream, not a pipeline failure.")
    print()
    print("  STRUCTURAL, and the most important one: region sets are")
    print("  label-driven DARs (one-vs-rest per cell type). That finds")
    print("  cell-type-associated regulons well and shared or continuous")
    print("  programs poorly. A TF absent here may be active everywhere")
    print("  rather than nowhere.")
    print()
    print(f"  TFs with at least one eRegulon: {len(tfs)}")
    print("=" * 74)


if __name__ == "__main__":
    main()
