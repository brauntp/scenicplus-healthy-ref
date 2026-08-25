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

# Curated human TF repertoire size (Lambert et al., Cell 2018). Used only as the
# denominator for the chance baseline on the lineage control.
N_HUMAN_TF = 1639


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


# TFs that share a motif with a family member, or that reach DNA only through a
# partner. cisTarget scores MOTIFS, so it cannot separate these: each expressed
# family member gets its own eRegulon from the same motif hits. Listing them is
# not a criticism of the pipeline -- it is the interpretive limit of motif
# enrichment, and reporting counts without it invites double-counting.
SHARED_MOTIF = {
    "RUNX": ["RUNX1", "RUNX2", "RUNX3"],
    "ELF (ETS)": ["ELF1", "ELF2", "ELF4"],
    "GABP (ETS)": ["GABPA", "GABPB1", "GABPB2"],
    "ETS broad": ["ETS1", "ETS2", "ELK1", "ELK3", "ELK4", "ERF", "ETV3"],
    "E-protein": ["TCF3", "TCF4", "TCF12"],
    "IRF": ["IRF1", "IRF2", "IRF8", "IRF9"],
}

# Present in a network only via a partner -- they make no direct DNA contact, so
# a motif-based method cannot have evidenced them independently.
NO_DIRECT_DNA_CONTACT = {
    "CBFB": "RUNX1/2/3 obligate partner",
    "GABPB1": "GABPA obligate partner",
    "GABPB2": "GABPA obligate partner",
    "EP300": "histone acetyltransferase co-activator",
    "CREBBP": "histone acetyltransferase co-activator",
    "TBL1XR1": "NCoR/SMRT corepressor subunit",
    "BPTF": "NURF remodelling subunit",
    "SMARCA4": "SWI/SNF ATPase",
}


# Canonical haematopoietic regulators, by compartment. A POSITIVE CONTROL: if
# these are absent the network is suspect regardless of how many eRegulons it
# contains. Chosen from established lineage biology, not from this run's output,
# so the check can fail.
LINEAGE_MARKERS = {
    "HSC / progenitor": ["RUNX1", "GATA2", "TAL1", "MYB", "ETV6", "MECOM", "HLF"],
    "erythroid":        ["GATA1", "KLF1", "TAL1", "NFE2", "GATA2"],
    "megakaryocyte":    ["GATA1", "FLI1", "RUNX1", "MEIS1"],
    "myeloid / mono":   ["SPI1", "CEBPA", "CEBPB", "CEBPE", "IRF8", "KLF4"],
    "granulocyte":      ["CEBPE", "GFI1", "SPI1"],
    "B lymphoid":       ["PAX5", "EBF1", "TCF3", "POU2AF1", "IRF4", "SPIB"],
    "T / NK":           ["TCF7", "GATA3", "LEF1", "RUNX3", "TBX21", "EOMES"],
    "pDC":              ["TCF4", "IRF8", "SPIB"],
}


def lineage_control(tfs: set[str]) -> None:
    """Did the network recover the TFs a haematopoietic reference must contain?"""
    print("-- positive control: canonical lineage regulators ---------------------")
    print("   These are chosen from lineage biology, not from this output, so")
    print("   the check can fail. Absence is the informative direction.")
    print()
    tot_hit = tot_all = 0
    for comp, members in LINEAGE_MARKERS.items():
        hit = [m for m in members if m in tfs]
        miss = [m for m in members if m not in tfs]
        tot_hit += len(hit); tot_all += len(members)
        print(f"    {comp:<18} {len(hit)}/{len(members)}")
        if hit:
            print(f"      found  : {', '.join(hit)}")
        if miss:
            print(f"      ABSENT : {', '.join(miss)}")
    print()
    print(f"    overall {tot_hit}/{tot_all} = {tot_hit/tot_all:.0%}")
    print()
    # A raw hit rate means nothing without a chance baseline: a network naming
    # enough TFs will hit markers by coincidence. Compare against drawing the
    # same NUMBER of TFs at random from the human repertoire.
    distinct = {m for ms in LINEAGE_MARKERS.values() for m in ms}
    n_hit = len(distinct & tfs)
    k = len(distinct)
    p_rand = min(1.0, len(tfs) / N_HUMAN_TF)
    exp = p_rand * k
    print(f"    CHANCE BASELINE. {len(tfs)} of ~{N_HUMAN_TF:,} human TFs "
          f"(Lambert 2018)")
    print(f"    have an eRegulon, so a random TF is recovered with "
          f"p = {p_rand:.3f}.")
    print(f"    distinct markers  : {k}")
    print(f"    expected by chance: {exp:.1f} ({p_rand:.0%})")
    print(f"    observed          : {n_hit} ({n_hit/k:.0%})")
    if n_hit > exp:
        try:
            from math import comb
            tail = sum(comb(k, i) * p_rand**i * (1 - p_rand)**(k - i)
                       for i in range(n_hit, k + 1))
            print(f"    P(>= {n_hit} of {k} by chance) = {tail:.2e}")
        except (OverflowError, ValueError):
            pass
    print("    Enrichment over chance is the result; the raw fraction is not.")
    print("      An absent marker has three innocent explanations before")
    print("      'the pipeline failed': its cell type was dropped for too few")
    print("      independent metacells; it acts through a motif shared with a")
    print("      recovered paralogue; or it is broadly active and so invisible")
    print("      to one-vs-rest DAR region sets.")
    print()


def motif_sharing_report(tfs: set[str]) -> None:
    """Which recovered TFs cannot be told apart, and which bind only via partners."""
    print("-- motif sharing: eRegulons that are not independent -----------------")
    any_hit = False
    for fam, members in SHARED_MOTIF.items():
        present = [m for m in members if m in tfs]
        if len(present) > 1:
            any_hit = True
            print(f"    {fam:<12} {', '.join(present)}")
    if not any_hit:
        print("    no families with more than one member recovered")
    print("      Members of a family share a motif, so their eRegulons draw on")
    print("      the SAME motif hits. Do not treat them as independent evidence,")
    print("      and do not sum their target counts.")
    print()
    cof = {t: why for t, why in NO_DIRECT_DNA_CONTACT.items() if t in tfs}
    if cof:
        print("    Recovered but making no direct DNA contact:")
        for t, why in sorted(cof.items()):
            print(f"      {t:<9} {why}")
        print("      These appear because a motif was annotated to them or to a")
        print("      partner. Their presence is not motif evidence for them.")
    print()


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
    print("  The two signs are INDEPENDENT: first is TF-to-gene correlation,")
    print("  second is region-to-gene. All four occur, and the previous version")
    print("  of this legend explained only two of them:")
    print("    +/+  TF rises with targets, regions open with them  -> activator")
    print("    -/-  TF falls with targets, regions close with them -> coherent")
    print("         repressor")
    print("    +/-  TF rises with targets but the regions CLOSE. Mixed. Can be a")
    print("         real repressive-looping mechanism, or the r2g correlation")
    print("         being driven by a different factor in the same locus.")
    print("    -/+  TF falls with targets while regions open. Same caution.")
    _mixed = sig.get("+/-", pd.Series(dtype=int)).sum() + \
             sig.get("-/+", pd.Series(dtype=int)).sum()
    _tot = int(sig.to_numpy().sum())
    if _tot:
        print(f"  MIXED-SIGN here: {_mixed}/{_tot} = {_mixed/_tot:.0%}. Treat those")
        print("  as hypotheses rather than called activators or repressors.")
    print()

    # -- per-eRegulon table ---------------------------------------------------
    tab = per_eregulon(df)
    args.out_prefix.parent.mkdir(parents=True, exist_ok=True)
    csv = args.out_prefix.with_suffix(".csv")
    tab.to_csv(csv)
    print(f"-- top {args.top} eRegulons by target-gene count -------------------")
    print("   CAUTION: ranking by target-gene count structurally favours broadly")
    print("   active and promoter-proximal factors over lineage-restricted ones.")
    print("   A big eRegulon is not a more important one. See the note below the")
    print("   table.")
    cols = [c for c in ("TF", "sign", "n_regions", "n_genes", "n_links",
                        "median_importance_x_abs_rho_R2G") if c in tab.columns]
    print(tab[cols].head(args.top).to_string())
    print()
    print(f"  full table ({len(tab)} eRegulons): {csv}")
    print()

    lineage_control(set(df["TF"]))
    motif_sharing_report(set(df["TF"]))

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
