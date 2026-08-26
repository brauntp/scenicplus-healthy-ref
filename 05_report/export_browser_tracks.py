#!/usr/bin/env python3
"""
Enhancer-promoter link tracks for a genome browser, per cell type.

WHY THIS EXISTS
---------------
`celltype_rho_usable.csv.gz` is 1.59M rows x 13 groups. A browser cannot render
that, and it carries no coordinates for the gene end of each link -- only a
signed `Distance`. This writes per-cell-type arc tracks that a browser can load.

THE TSS ANCHOR, AND WHY IT IS RECOVERED RATHER THAN LOOKED UP
-------------------------------------------------------------
An arc needs two endpoints: the peak (which we have) and the gene's TSS (which
we do not). The obvious route -- look the gene up in a fresh
`genome_annotation.tsv` -- does NOT reproduce the pipeline's own geometry:

  * `get_search_space` (use_gene_boundaries=False by default) builds a promoter
    per TRANSCRIPT, not per gene. hg38 knownGene has a median of 4 distinct TSS
    per gene and up to 145, and 86.7% of genes have more than one.
  * Measured on 4,000 real links against a freshly built annotation: SOME
    transcript's TSS reproduces the stored |Distance| for only 44.9% of links,
    and the min-over-transcripts rule for 29.2%. Interval-to-point distance is
    ruled out entirely (1.0%) -- the distance is peak-MIDPOINT to TSS.

So a fresh annotation cannot tell you which transcript anchored a given link,
and picking the wrong one silently misplaces the arc by a median 1.25 kb.

What IS exact: `Distance` was computed from the peak midpoint, so the TSS sits
at `midpoint +/- |Distance|`. That leaves two candidates; this script picks the
one nearer the gene's annotated TSS set. The result reproduces the stored
|Distance| EXACTLY by construction (verified: 100% of rows), and lands on an
annotated TSS exactly for 44.9% of links, within 5 kb for 61.9%, median error
1.25 kb.

Do not "fix" this by looking the TSS up directly -- that trades an exact
reconstruction of the pipeline's geometry for a plausible-looking wrong one.

`Distance`'s SIGN is strand-relative (upstream/downstream in gene orientation),
which is why the sign alone cannot orient the arc and the candidate test is
needed.

OUTPUTS
-------
Per cell-type group, two formats because browsers differ:

  <group>.interact       UCSC interact (bigInteract source). Arcs with score and
                         colour. Load via a track line, or convert with
                         `bedToBigBed -as=interact.as -type=bed5+13`.
  <group>.bedpe          IGV / generic paired-end. Simpler, no colour semantics.

Plus `<group>.peaks.bed` -- just the peak ends, useful as a companion feature
track and as the interval set for the bigWig quantification.

SCORE AND COLOUR
----------------
`score` is 0-1000 from `specificity_z` (the validated ranking statistic), NOT
from rho. Colour encodes the sign of `best_rho`: blue for positive coupling
(candidate enhancer), orange for negative (candidate silencer). A browser user
therefore sees "how cell-type-specific" as intensity and "which direction" as
hue.

Usage:
    python 05_report/export_browser_tracks.py \\
        --links celltype_rho_usable.csv.gz \\
        --annotation anno/genome_annotation.tsv \\
        --out-dir browser_tracks \\
        --min-specificity-z 3 --min-abs-rho 0.15 --max-per-group 50000
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

# Blue = positive coupling (candidate enhancer); orange = negative (silencer).
# Chosen to survive greyscale printing and the common red-green colour vision
# deficiencies, which a red/green pair would not.
COL_POS = "31,111,180"
COL_NEG = "201,97,31"


def parse_regions(regions: pd.Series) -> pd.DataFrame:
    """chr:start-end -> columns. Fails loudly rather than dropping bad rows."""
    ex = regions.str.extract(r"^(?P<chrom>[^:]+):(?P<start>\d+)-(?P<end>\d+)$")
    bad = ex.chrom.isna()
    if bad.any():
        raise SystemExit(
            f"ERROR: {int(bad.sum())} region ids are not 'chr:start-end', "
            f"e.g. {regions[bad].head(3).tolist()}")
    return pd.DataFrame({
        "chrom": ex.chrom.to_numpy(),
        "start": ex.start.astype(np.int64).to_numpy(),
        "end": ex.end.astype(np.int64).to_numpy()})


def unwrap_distance(s: pd.Series) -> pd.Series:
    """Upstream may ship Distance as '[-126154]'. pd.to_numeric on that returns
    all-NaN SILENTLY, so handle both forms explicitly."""
    if s.dtype.kind in "iuf":
        return s
    return (s.astype(str).str.strip("[]").replace("", np.nan)
            .astype("Float64").astype("Int64"))


def recover_tss(mid: np.ndarray, dist_abs: np.ndarray, genes: np.ndarray,
                tss_by_gene: dict[str, np.ndarray]) -> tuple[np.ndarray, np.ndarray]:
    """TSS = mid +/- |Distance|; choose the candidate nearer the gene's
    annotated TSS set. Returns (tss, error_to_nearest_annotated)."""
    lo, hi = mid - dist_abs, mid + dist_abs
    out = np.where(dist_abs == 0, mid, lo).astype(np.int64)
    err = np.full(len(mid), np.nan)
    for i, g in enumerate(genes):
        arr = tss_by_gene.get(g)
        if arr is None:
            continue
        dl = np.abs(arr - lo[i]).min()
        dh = np.abs(arr - hi[i]).min()
        out[i] = lo[i] if dl <= dh else hi[i]
        err[i] = min(dl, dh)
    return out, err


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--links", type=Path, required=True,
                    help="celltype_rho_usable.csv.gz")
    ap.add_argument("--annotation", type=Path, required=True,
                    help="genome_annotation.tsv from build_genome_annotation.py")
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--min-specificity-z", type=float, default=3.0)
    ap.add_argument("--min-abs-rho", type=float, default=0.15)
    ap.add_argument("--max-per-group", type=int, default=50_000,
                    help="top N by specificity_z per group; 0 = no cap")
    ap.add_argument("--groups", nargs="*", default=None,
                    help="restrict to these groups (default: all present)")
    ap.add_argument("--dry-run", action="store_true",
                    help="report the plan and write nothing")
    args = ap.parse_args()

    for p in (args.links, args.annotation):
        if not p.exists():
            sys.exit(f"ERROR: not found: {p}")

    print("=" * 74)
    print("enhancer-promoter link tracks")
    print("=" * 74)

    ga = pd.read_csv(args.annotation, sep="\t")
    need = {"Gene", "Transcription_Start_Site", "Strand", "Chromosome"}
    if not need <= set(ga.columns):
        sys.exit(f"ERROR: {args.annotation} lacks {sorted(need - set(ga.columns))}")
    tss_by_gene = {g: v.to_numpy(np.int64)
                   for g, v in ga.groupby("Gene").Transcription_Start_Site}
    print(f"  annotation : {len(ga):,} transcripts, {len(tss_by_gene):,} genes")

    usecols = ["region", "target", "Distance", "best_group",
               "specificity_z", "best_rho"]
    lk = pd.read_csv(args.links, usecols=usecols)
    lk["Distance"] = unwrap_distance(lk["Distance"])
    print(f"  links      : {len(lk):,} pairs, {lk.best_group.nunique()} groups")

    keep = ((lk.specificity_z >= args.min_specificity_z)
            & (lk.best_rho.abs() >= args.min_abs_rho)
            & lk.Distance.notna())
    lk = lk[keep].copy()
    print(f"  after specificity_z>={args.min_specificity_z} and "
          f"|best_rho|>={args.min_abs_rho}: {len(lk):,}")
    if not len(lk):
        sys.exit("ERROR: no link passes the filters. Lower --min-specificity-z "
                 "or --min-abs-rho.")

    groups = args.groups or sorted(lk.best_group.unique())
    missing = [g for g in groups if g not in set(lk.best_group)]
    if missing:
        sys.exit(f"ERROR: requested group(s) absent after filtering: {missing}\n"
                 f"       present: {sorted(lk.best_group.unique())}")

    print()
    print(f"{'group':<26}{'links':>10}{'written':>10}")
    plan = {}
    for g in groups:
        sub = lk[lk.best_group == g]
        n_out = len(sub) if not args.max_per_group else min(len(sub), args.max_per_group)
        plan[g] = (len(sub), n_out)
        print(f"{g:<26}{len(sub):>10,}{n_out:>10,}")
    print(f"{'TOTAL':<26}{sum(v[0] for v in plan.values()):>10,}"
          f"{sum(v[1] for v in plan.values()):>10,}")

    if args.dry_run:
        print()
        print("--dry-run: nothing written.")
        print("=" * 74)
        return

    args.out_dir.mkdir(parents=True, exist_ok=True)
    err_all = []
    for g in groups:
        sub = lk[lk.best_group == g]
        if args.max_per_group:
            sub = sub.nlargest(args.max_per_group, "specificity_z")
        sub = sub.reset_index(drop=True)

        co = parse_regions(sub.region)
        mid = ((co.start + co.end) // 2).to_numpy()
        dist_abs = sub.Distance.abs().astype(np.int64).to_numpy()
        tss, err = recover_tss(mid, dist_abs, sub.target.to_numpy(), tss_by_gene)
        err_all.append(err)

        # Score 0-1000 from specificity_z. Scale by the group's own max so each
        # track uses the full dynamic range; the absolute z is in the name field.
        z = sub.specificity_z.to_numpy()
        score = np.clip((z / max(z.max(), 1e-9) * 1000).round(), 0, 1000).astype(int)
        colour = np.where(sub.best_rho.to_numpy() > 0, COL_POS, COL_NEG)

        lo = np.minimum(mid, tss)
        hi = np.maximum(mid, tss)
        safe = g.replace(" ", "_").replace("/", "-")

        # -- UCSC interact ----------------------------------------------------
        # bed5+13: chrom start end name score value exp color
        #          sourceChrom sourceStart sourceEnd sourceName sourceStrand
        #          targetChrom targetStart targetEnd targetName targetStrand
        inter = pd.DataFrame({
            "chrom": co.chrom, "chromStart": lo, "chromEnd": hi + 1,
            "name": [f"{t}|z={zz:.1f}|rho={rr:+.2f}"
                     for t, zz, rr in zip(sub.target, z, sub.best_rho)],
            "score": score, "value": np.round(z, 3), "exp": safe, "color": colour,
            "sourceChrom": co.chrom, "sourceStart": co.start, "sourceEnd": co.end,
            "sourceName": sub.region, "sourceStrand": ".",
            "targetChrom": co.chrom, "targetStart": tss, "targetEnd": tss + 1,
            "targetName": sub.target, "targetStrand": ".",
        }).sort_values(["chrom", "chromStart"])
        pi = args.out_dir / f"{safe}.interact"
        with open(pi, "w") as fh:
            fh.write(f'track type=interact name="{g} E-P links" '
                     f'description="specificity_z>={args.min_specificity_z}, '
                     f'|best_rho|>={args.min_abs_rho}" '
                     'interactDirectional=true maxHeightPixels=200:100:50 '
                     'visibility=full\n')
            inter.to_csv(fh, sep="\t", header=False, index=False)

        # -- BEDPE ------------------------------------------------------------
        bedpe = pd.DataFrame({
            "chrom1": co.chrom, "start1": co.start, "end1": co.end,
            "chrom2": co.chrom, "start2": tss, "end2": tss + 1,
            "name": sub.target, "score": score,
            "strand1": ".", "strand2": ".",
        }).sort_values(["chrom1", "start1"])
        bedpe.to_csv(args.out_dir / f"{safe}.bedpe", sep="\t",
                     header=False, index=False)

        # -- peak ends, as a companion feature track --------------------------
        pk = (pd.DataFrame({"chrom": co.chrom, "start": co.start, "end": co.end,
                            "name": sub.region, "score": score, "strand": "."})
              .drop_duplicates(subset=["chrom", "start", "end"])
              .sort_values(["chrom", "start"]))
        pk.to_csv(args.out_dir / f"{safe}.peaks.bed", sep="\t",
                  header=False, index=False)

        print(f"  {safe:<26} {len(inter):>7,} arcs, {len(pk):>7,} peaks")

    e = np.concatenate(err_all)
    e = e[~np.isnan(e)]
    print()
    print("TSS anchor quality (recovered position vs the gene's annotated TSS set):")
    print(f"  exact          : {(e == 0).mean():.1%}")
    print(f"  within 500 bp  : {(e <= 500).mean():.1%}")
    print(f"  within 5 kb    : {(e <= 5000).mean():.1%}")
    print(f"  median error   : {np.median(e):,.0f} bp")
    print("  The anchor reproduces the stored |Distance| exactly by construction;")
    print("  the error above is against a FRESHLY BUILT annotation, which cannot")
    print("  identify which transcript the pipeline used (see the module docstring).")
    print()
    print(f"  wrote {args.out_dir}")
    print("=" * 74)


if __name__ == "__main__":
    main()
