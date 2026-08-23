#!/usr/bin/env python3
"""
Positive control for the paired object: does it recover known biology?

WHY THIS EXISTS
---------------
`pairing_diagnostics.csv` reports how close RNA and ATAC cells sat in the GLUE
latent space. That is a property of the EMBEDDING, and it can look excellent
while the paired object is still useless -- if labels were misassigned, if the
positional merge was wrong, or if the two modalities were paired at the right
distance but the wrong cells, the cross-modal gap would not notice.

This script asks a question the diagnostics cannot: across metacells, does
marker-gene expression track accessibility at that gene's own locus? For a
lineage-defining gene the answer is emphatically yes in real data -- GATA1
expression and GATA1-proximal accessibility both rise in MEP/erythroid and are
off elsewhere. If that correlation is absent, the pairing failed regardless of
what the gaps say.

Three checks, cheapest first:

  1. STRUCTURE -- shapes, obs_names identity across modalities, group counts
     against the expected 24 groups / 25,323 metacells, NaN/zero-row scan.

  2. MARKER SPECIFICITY (RNA only) -- is each marker's expression highest in
     the group it marks? This tests the LABELS, independently of pairing.
     Failure here means the label merge or the group key is wrong.

  3. LINKAGE (the real control) -- for each marker, correlate its expression
     with mean accessibility of peaks within +/-`--window` of its TSS, across
     metacells. Compared against a null built from DISTANT peaks matched on
     accessibility, so the number means something. If true-locus correlations
     are not shifted above the null, region-to-gene inference will not work on
     this object and there is no point running SCENIC+.

Marker TSS coordinates come from 00_inspect/marker_tss_hg38.tsv, fetched from
Ensembl (GRCh38) rather than typed from memory.

Usage
-----
    python 02_pair/qc_paired.py --h5mu ACC_GEX.h5mu \\
        --markers 00_inspect/marker_tss_hg38.tsv \\
        --group-key predicted_CellType_Broad --out qc_paired

Writes <out>.md, <out>.json and <out>.png. Read the markdown.
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np


def _log(m):
    print(f"[qc] {m}", flush=True)


def parse_regions(names):
    """chr:start-end (or chr-start-end) -> arrays. Unparseable -> masked out."""
    chrom = np.empty(len(names), dtype=object)
    start = np.full(len(names), -1, dtype=np.int64)
    end = np.full(len(names), -1, dtype=np.int64)
    for i, n in enumerate(names):
        s = str(n).replace("-", ":", 1) if ":" not in str(n) else str(n)
        try:
            c, rest = s.split(":", 1)
            a, b = rest.split("-", 1)
            chrom[i] = c
            start[i] = int(a)
            end[i] = int(b)
        except Exception:
            chrom[i] = None
    return chrom, start, end


def spearman(x, y):
    """Rank correlation without scipy, so this runs in the pairing env."""
    n = len(x)
    if n < 3:
        return float("nan")
    rx = np.argsort(np.argsort(x)).astype(np.float64)
    ry = np.argsort(np.argsort(y)).astype(np.float64)
    rx -= rx.mean()
    ry -= ry.mean()
    d = np.sqrt((rx ** 2).sum() * (ry ** 2).sum())
    return float((rx * ry).sum() / d) if d > 0 else float("nan")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--h5mu", type=Path, required=True)
    ap.add_argument("--markers", type=Path,
                    default=Path("00_inspect/marker_tss_hg38.tsv"))
    ap.add_argument("--group-key", default="predicted_CellType_Broad")
    ap.add_argument("--window", type=int, default=150_000,
                    help="bp either side of the TSS; matches SCENIC+'s default "
                         "region-to-gene search space")
    ap.add_argument("--expect-groups", type=int, default=None)
    ap.add_argument("--expect-metacells", type=int, default=None)
    ap.add_argument("--out", type=Path, default=Path("qc_paired"))
    ap.add_argument("--no-plot", action="store_true")
    args = ap.parse_args()

    if not args.h5mu.exists():
        sys.exit(f"ERROR: {args.h5mu} does not exist. Did the pairing job finish?")
    if not args.markers.exists():
        sys.exit(f"ERROR: marker table {args.markers} not found (it ships in the repo)")

    import mudata
    import pandas as pd

    # backed=True keeps the matrices on disk. The ATAC block is ~37 GB dense on
    # this reference and QC needs only the few hundred marker-window columns, so
    # loading it whole would cost ~90 GB for no benefit.
    _log(f"reading {args.h5mu} (backed)")
    try:
        md = mudata.read(str(args.h5mu), backed="r")
    except (TypeError, ValueError, OSError) as e:
        _log(f"backed read unavailable ({type(e).__name__}); reading in memory")
        _log("NOTE: this materialises every matrix -- submit it "
             "(slurm/qc_paired.sbatch) rather than running on a login node")
        md = mudata.read(str(args.h5mu))
    if not {"scRNA", "scATAC"} <= set(md.mod):
        sys.exit(f"ERROR: expected modalities scRNA and scATAC, found {list(md.mod)}")
    R, A = md["scRNA"], md["scATAC"]
    rep = {"file": str(args.h5mu), "checks": {}}

    # ---------------------------------------------------------------- 1 STRUCTURE
    s = {"n_metacells_rna": int(R.shape[0]), "n_metacells_atac": int(A.shape[0]),
         "n_genes": int(R.shape[1]), "n_regions": int(A.shape[1]),
         "obs_names_identical": list(R.obs_names) == list(A.obs_names)}
    def _dense(X):
        return X.toarray() if hasattr(X, "toarray") else np.asarray(X)

    # RNA is small (~3 GB dense here) so read it whole. ATAC is not: only the
    # marker-window columns are ever needed, selected further down.
    Xr = _dense(R.X[:] if hasattr(R.X, "shape") and not hasattr(R.X, "toarray")
                else R.X)
    s["rna_nan"] = int(np.isnan(Xr).sum())
    s["rna_allzero_metacells"] = int((Xr.sum(1) == 0).sum())
    # ATAC row sums stream in metacell blocks -- no full densification.
    _atac_rowsum = np.zeros(A.shape[0], dtype=np.float64)
    _atac_nan = 0
    _BLK = 2048
    for i0 in range(0, A.shape[0], _BLK):
        blk = _dense(A.X[i0:i0 + _BLK])
        _atac_nan += int(np.isnan(blk).sum())
        _atac_rowsum[i0:i0 + blk.shape[0]] = blk.sum(1)
    s["atac_nan"] = _atac_nan
    s["atac_allzero_metacells"] = int((_atac_rowsum == 0).sum())
    s["atac_scan"] = "streamed in blocks (full scan, never fully densified)"
    if args.group_key in R.obs.columns:
        vc = R.obs[args.group_key].astype(str).value_counts()
        s["n_groups"] = int(len(vc))
        s["metacells_per_group"] = {k: int(v) for k, v in vc.items()}
        junk = [g for g in vc.index if g.strip().lower() in
                {"nan", "", "na", "none", "unassigned", "unknown"}]
        s["junk_groups"] = junk
    else:
        s["n_groups"] = None
        s["group_key_missing"] = True
    if args.expect_groups:
        s["groups_as_expected"] = (s.get("n_groups") == args.expect_groups)
    if args.expect_metacells:
        s["metacells_as_expected"] = (s["n_metacells_rna"] == args.expect_metacells)
    rep["checks"]["structure"] = s
    _log(f"structure: {s['n_metacells_rna']:,} metacells, {s.get('n_groups')} groups, "
         f"obs_names identical={s['obs_names_identical']}")

    mk = pd.read_csv(args.markers, sep="\t")
    genes_present = [g for g in mk.gene.unique() if g in set(R.var_names)]
    _log(f"markers: {len(genes_present)}/{mk.gene.nunique()} present in the object")

    # -------------------------------------------------- 2 MARKER SPECIFICITY (RNA)
    spec_rows = []
    if s.get("n_groups"):
        grp = R.obs[args.group_key].astype(str).to_numpy()
        gi = {g: np.flatnonzero(grp == g) for g in np.unique(grp)}
        # CP10K-like normalisation so groups with deeper metacells do not dominate
        tot = Xr.sum(1, keepdims=True)
        Xn = np.log1p(np.divide(Xr, np.where(tot == 0, 1, tot)) * 1e4)
        gcol = {g: i for i, g in enumerate(R.var_names)}
        # A marker listed for n groups cannot rank first in all of them: GATA1
        # genuinely marks both MEP and Early Erythroid. Score such a gene as
        # correct if the group is within the top n, so a shared marker is not
        # counted as a failure.
        n_listed = mk.groupby("gene").group.nunique().to_dict()
        for _, r in mk.iterrows():
            if r.gene not in gcol or r.group not in gi:
                continue
            v = Xn[:, gcol[r.gene]]
            means = {g: float(v[idx].mean()) for g, idx in gi.items() if len(idx)}
            order = sorted(means, key=lambda g: -means[g])
            rank = order.index(r.group) + 1
            nl = int(n_listed.get(r.gene, 1))
            spec_rows.append({"gene": r.gene, "marks_group": r.group,
                              "rank_of_marked_group": rank,
                              "n_groups_listed_for": nl,
                              "within_listed_top": bool(rank <= nl),
                              "top_group": order[0],
                              "mean_in_marked": round(means[r.group], 4),
                              "max_mean": round(means[order[0]], 4)})
        ranks = np.array([r["rank_of_marked_group"] for r in spec_rows])
        within = np.array([r["within_listed_top"] for r in spec_rows])
        rep["checks"]["marker_specificity"] = {
            "n_tested": int(len(ranks)),
            "frac_top1": float((ranks == 1).mean()) if len(ranks) else None,
            "frac_top3": float((ranks <= 3).mean()) if len(ranks) else None,
            # The headline number: correct within the number of groups the gene
            # is listed for, so shared markers are not penalised.
            "frac_within_listed": float(within.mean()) if len(ranks) else None,
            "median_rank": float(np.median(ranks)) if len(ranks) else None,
            "n_groups": int(s["n_groups"]),
            "detail": spec_rows}
        if len(ranks):
            _log(f"marker specificity: {(ranks==1).mean():.0%} rank-1, "
                 f"{(ranks<=3).mean():.0%} top-3 of {s['n_groups']} groups")

    # ------------------------------------------------------------- 3 LINKAGE
    chrom, rs, re_ = parse_regions(list(A.var_names))
    ok = np.array([c is not None for c in chrom])
    _log(f"regions parsed: {ok.sum():,}/{len(ok):,}")
    mid = np.where(ok, (rs + re_) / 2.0, np.nan)
    # Accessibility strata for the null: match decoy peaks on mean signal, so a
    # correlation is not attributable to peaks simply being more accessible.
    # Column means, streamed. Needed for the accessibility-matched null, and it
    # is the last thing that would otherwise force a full densification.
    _csum = np.zeros(A.shape[1], dtype=np.float64)
    for i0 in range(0, A.shape[0], _BLK):
        _csum += _dense(A.X[i0:i0 + _BLK]).sum(0)
    acc_mean = (_csum / A.shape[0]).astype(np.float32)

    # Only the columns QC actually touches get materialised: peaks within the
    # window of some marker TSS, plus a pool of accessibility-matched decoys.
    # That is a few thousand of ~394k columns, so this is MBs not tens of GBs.
    _need = np.zeros(A.shape[1], dtype=bool)
    for gene in genes_present:
        _r = mk[mk.gene == gene].iloc[0]
        _need |= (ok & (chrom == _r.chrom) & (np.abs(mid - _r.tss) <= args.window))
    _n_true = int(_need.sum())
    _rng0 = np.random.default_rng(12345)
    _pool = np.flatnonzero(ok & ~_need)
    if len(_pool):
        _decoy = _rng0.choice(_pool, min(len(_pool), max(20_000, 40 * _n_true)),
                              replace=False)
        _need[_decoy] = True
    _cols = np.flatnonzero(_need)
    _log(f"materialising {len(_cols):,} of {A.shape[1]:,} ATAC columns "
         f"({_n_true:,} in marker windows + decoy pool) = "
         f"{A.shape[0]*len(_cols)*4/1024**3:.2f} GB")
    _remap = np.full(A.shape[1], -1, dtype=np.int64)
    _remap[_cols] = np.arange(len(_cols))
    Xa_sub = np.empty((A.shape[0], len(_cols)), dtype=np.float32)
    for i0 in range(0, A.shape[0], _BLK):
        Xa_sub[i0:i0 + _BLK] = _dense(A.X[i0:i0 + _BLK])[:, _cols]
    rng = np.random.default_rng(0)
    gcol = {g: i for i, g in enumerate(R.var_names)}
    tot = Xr.sum(1, keepdims=True)
    Xn = np.log1p(np.divide(Xr, np.where(tot == 0, 1, tot)) * 1e4)

    link_rows = []
    for gene in genes_present:
        row = mk[mk.gene == gene].iloc[0]
        near = ok & (chrom == row.chrom) & (np.abs(mid - row.tss) <= args.window)
        n_near = int(near.sum())
        if n_near == 0:
            link_rows.append({"gene": gene, "n_peaks_in_window": 0,
                              "rho_true": None, "rho_null_median": None})
            continue
        expr = Xn[:, gcol[gene]]
        acc_true = Xa_sub[:, _remap[np.flatnonzero(near)]].mean(1)
        rho_true = spearman(acc_true, expr)
        # null: same number of peaks, same accessibility stratum, other chromosomes
        lo, hi = np.quantile(acc_mean[near], [0.1, 0.9])
        # restrict decoys to columns we materialised
        far = (ok & (chrom != row.chrom) & (acc_mean >= lo)
               & (acc_mean <= hi) & (_remap >= 0))
        far_idx = np.flatnonzero(far)
        nulls = []
        if len(far_idx) >= n_near:
            for _ in range(30):
                pick = rng.choice(far_idx, n_near, replace=False)
                nulls.append(spearman(Xa_sub[:, _remap[pick]].mean(1), expr))
        nulls = np.array([x for x in nulls if np.isfinite(x)])
        link_rows.append({
            "gene": gene, "chrom": row.chrom, "tss": int(row.tss),
            "marks_group": row.group, "n_peaks_in_window": n_near,
            "rho_true": None if not np.isfinite(rho_true) else round(float(rho_true), 4),
            "rho_null_median": None if not len(nulls) else round(float(np.median(nulls)), 4),
            "rho_null_p95": None if not len(nulls) else round(float(np.quantile(nulls, 0.95)), 4),
            "beats_null_p95": bool(len(nulls) and np.isfinite(rho_true)
                                   and rho_true > np.quantile(nulls, 0.95))})

    tested = [r for r in link_rows if r.get("rho_true") is not None]
    rt = np.array([r["rho_true"] for r in tested]) if tested else np.array([])
    rn = np.array([r["rho_null_median"] for r in tested
                   if r["rho_null_median"] is not None])
    rep["checks"]["linkage"] = {
        "window_bp": args.window,
        "n_markers_tested": len(tested),
        "n_markers_no_peaks": sum(1 for r in link_rows if r["n_peaks_in_window"] == 0),
        "median_rho_true": None if not len(rt) else round(float(np.median(rt)), 4),
        "median_rho_null": None if not len(rn) else round(float(np.median(rn)), 4),
        "frac_beating_null_p95": None if not tested else
            round(float(np.mean([r["beats_null_p95"] for r in tested])), 4),
        "detail": link_rows}
    if len(rt):
        _log(f"linkage: median rho_true {np.median(rt):.3f} vs null "
             f"{np.median(rn):.3f}; "
             f"{np.mean([r['beats_null_p95'] for r in tested]):.0%} beat null p95")

    # ------------------------------------------------------------------ verdict
    v, why = "PASS", []
    st = rep["checks"]["structure"]
    if not st["obs_names_identical"]:
        v, _ = "FAIL", why.append("RNA and ATAC obs_names differ -- not a paired object")
    if st["rna_nan"] or st["atac_nan"]:
        v = "FAIL"; why.append("NaNs in a metacell matrix")
    if st.get("junk_groups"):
        v = "FAIL"; why.append(f"junk group(s) present: {st['junk_groups']}")
    ms = rep["checks"].get("marker_specificity", {})
    _spec = ms.get("frac_within_listed", ms.get("frac_top3"))
    if _spec is not None and _spec < 0.5:
        v = "FAIL" if v != "FAIL" else v
        why.append(f"only {_spec:.0%} of markers are top-ranked in a group they "
                   "mark -- suspect the labels, not the pairing")
    lk = rep["checks"]["linkage"]
    if lk["frac_beating_null_p95"] is not None:
        # Judge against CHANCE, not a fixed fraction. Beating one's own null 95th
        # percentile happens 5% of the time under the null, so the question is
        # whether the observed rate is enriched over 0.05 -- a binomial test.
        #
        # A fixed threshold was wrong here: it was calibrated on a simulation
        # whose decoy median rho was ~0, but real metacells are grouped BY CELL
        # TYPE, and peaks and genes share lineage programs, so a distal peak
        # correlates with an unrelated gene at rho ~0.2. Verified by simulation
        # that this inflates the absolute rho of BOTH true and decoy sets while
        # leaving the per-gene beat rate calibrated at ~5-9%: the absolute
        # numbers do not transfer between datasets, but the enrichment does.
        n_t = lk["n_markers_tested"]
        frac = lk["frac_beating_null_p95"]
        n_hit = int(round(frac * n_t))
        try:
            from scipy.stats import binomtest
            p_enrich = float(binomtest(n_hit, n_t, 0.05,
                                       alternative="greater").pvalue)
        except Exception:                                      # pragma: no cover
            p_enrich = None
        lk["n_beating_null"] = n_hit
        lk["enrichment_over_chance"] = round(frac / 0.05, 2)
        lk["binomial_p_vs_chance"] = p_enrich
        if p_enrich is not None and p_enrich < 1e-4 and frac >= 0.15:
            pass                                    # decisive locus-specific signal
        elif p_enrich is not None and p_enrich < 0.05:
            if v == "PASS":
                v = "WARN"
            why.append(
                f"{n_hit}/{n_t} marker loci beat their own null p95 "
                f"({frac:.0%}, {frac/0.05:.1f}x chance, p={p_enrich:.1e}) -- "
                "significant but modest; distal links deserve more scepticism "
                "than proximal ones")
        else:
            v = "FAIL"
            why.append(
                f"only {n_hit}/{n_t} marker loci beat their own null p95 "
                f"({frac:.0%}, {frac/0.05:.1f}x chance"
                + (f", p={p_enrich:.2g}" if p_enrich is not None else "")
                + ") -- not distinguishable from chance, so region-to-gene "
                "inference will not work on this object")
        # The absolute medians are informative but must not be read as effect
        # sizes: a high decoy median means shared cell-type structure, not a
        # broken pairing.
        if (lk["median_rho_null"] is not None
                and lk["median_rho_null"] > 0.10):
            why.append(
                f"note: decoy median rho is {lk['median_rho_null']} -- peaks and "
                "genes share strong lineage structure, so absolute rho is "
                "inflated for true AND decoy loci alike. Judge by the "
                "enrichment above, not by the raw correlation.")
    rep["verdict"] = v
    rep["reasons"] = why

    # ------------------------------------------------------------------- outputs
    L = [f"# Paired-object QC: `{args.h5mu.name}`", "",
         f"**Verdict: {v}**", ""]
    for w in why:
        L.append(f"- {w}")
    if not why:
        L.append("- structure, marker specificity and locus linkage all as expected")
    L += ["", "## 1. Structure", "",
          f"- metacells: {st['n_metacells_rna']:,} (RNA) / {st['n_metacells_atac']:,} (ATAC)",
          f"- features: {st['n_genes']:,} genes x {st['n_regions']:,} regions",
          f"- obs_names identical across modalities: **{st['obs_names_identical']}**",
          f"- groups: {st.get('n_groups')}"]
    if st.get("groups_as_expected") is not None:
        L.append(f"- matches expected group count: **{st['groups_as_expected']}**")
    if st.get("metacells_as_expected") is not None:
        L.append(f"- matches expected metacell count: **{st['metacells_as_expected']}**")
    L.append(f"- all-zero metacells: {st['rna_allzero_metacells']} RNA, "
             f"{st['atac_allzero_metacells']} ATAC")
    if ms:
        L += ["", "## 2. Marker specificity (tests the LABELS)", "",
              f"Of {ms['n_tested']} marker-group pairs across {ms['n_groups']} groups: "
              f"**{ms['frac_within_listed']:.0%} correct** (the marked group ranks "
              f"within the number of groups that marker is listed for -- GATA1 is "
              f"listed for MEP and Early Erythroid, so top-2 counts as correct). "
              f"Strict rank-1: {ms['frac_top1']:.0%}; median rank "
              f"{ms['median_rank']:.0f}.", "",
              "Markers that miss even that allowance:", ""]
        off = [r for r in ms["detail"] if not r["within_listed_top"]]
        L += ([f"- `{r['gene']}` marks *{r['marks_group']}* (listed for "
               f"{r['n_groups_listed_for']}) but ranks {r['rank_of_marked_group']}; "
               f"highest in *{r['top_group']}*" for r in off[:15]]
              or ["- none"])
    L += ["", "## 3. Locus linkage (tests the PAIRING)", "",
          f"Expression vs mean accessibility within ±{lk['window_bp']//1000} kb of "
          "each marker's TSS, across metacells, against a null of "
          "accessibility-matched peaks on other chromosomes.", ""]
    if lk["median_rho_true"] is not None:
        L += [f"- median rho at the true locus: **{lk['median_rho_true']}**",
              f"- median rho for matched decoys: {lk['median_rho_null']}"
              + ("  *(inflated by shared lineage structure -- see the note above)*"
                 if (lk['median_rho_null'] or 0) > 0.10 else ""),
              f"- markers beating their own null 95th percentile: "
              f"**{lk.get('n_beating_null', '?')}/{lk['n_markers_tested']} "
              f"= {lk['frac_beating_null_p95']:.0%}**"
              + (f", i.e. **{lk['enrichment_over_chance']}x chance** "
                 f"(binomial p={lk['binomial_p_vs_chance']:.2g})"
                 if lk.get('binomial_p_vs_chance') is not None else ""),
              "",
              "  The 5% expected under the null is the reference, not 50%: each "
              "gene is tested against its OWN decoys, so the enrichment over "
              "chance is the quantity that transfers between datasets. Absolute "
              "rho does not.",
              "",
              f"- markers with no peak in window: {lk['n_markers_no_peaks']}", ""]
        top = sorted([r for r in lk["detail"] if r.get("rho_true") is not None],
                     key=lambda r: -r["rho_true"])[:10]
        L += ["Strongest loci:", "",
              "| gene | marks | peaks in window | rho | null p95 |",
              "|---|---|---|---|---|"]
        L += [f"| {r['gene']} | {r['marks_group']} | {r['n_peaks_in_window']} | "
              f"{r['rho_true']} | {r['rho_null_p95']} |" for r in top]
    else:
        L.append("- no marker locus could be tested; check region name format")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    Path(f"{args.out}.md").write_text("\n".join(L))
    with open(f"{args.out}.json", "w") as fh:
        json.dump(rep, fh, indent=2, default=str)

    if not args.no_plot and len(rt):
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            fig, ax = plt.subplots(figsize=(5.0, 3.2))
            ax.hist(rn, bins=20, alpha=0.65, label="accessibility-matched decoys",
                    color="#9aa0a6")
            ax.hist(rt, bins=20, alpha=0.8, label="true marker locus", color="#08306b")
            ax.set_xlabel("Spearman rho, expression vs local accessibility")
            ax.set_ylabel("markers")
            ax.set_title("Locus linkage in the paired object", loc="left")
            ax.legend(frameon=False, fontsize=7)
            for sp in ("top", "right"):
                ax.spines[sp].set_visible(False)
            fig.tight_layout()
            fig.savefig(f"{args.out}.png", dpi=200)
            _log(f"wrote {args.out}.png")
        except Exception as e:                                  # pragma: no cover
            _log(f"plot skipped: {type(e).__name__}: {e}")

    print()
    print("\n".join(L))
    print(f"\nwrote {args.out}.md / .json")
    sys.exit(0 if v != "FAIL" else 2)


if __name__ == "__main__":
    main()
