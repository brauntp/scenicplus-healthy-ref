#!/usr/bin/env python3
"""
Prevalence-stratified re-analysis of dropped vs retained peaks.

WHY THIS EXISTS
---------------
characterize_dropped.py compared the specificity of dropped and retained peaks
directly and concluded "custom database justified" because dropped peaks scored
higher (4.80 vs 4.16). That conclusion was unsound.

`specificity = max-group mean / mean over groups` is NOT prevalence-invariant.
A peak accessible in 230 cells spread over 24 groups of wildly different sizes
will concentrate in some group by chance alone. Simulating peaks with ZERO true
cell-type preference at this reference's real group sizes:

    prevalence 0.0030 (retained median) -> apparent specificity 2.55
    prevalence 0.0014 (dropped  median) -> apparent specificity 3.27

i.e. the null alone predicts dropped peaks scoring 1.28x higher, while the
observed ratio is only 1.15x. The specificity gap is entirely a side effect of
dropped peaks being sparser -- and the sparseness is the larger, opposite-signed
effect (rank-biserial -0.49 on prevalence).

THE FIX
-------
Compare only peaks of comparable prevalence. Within each prevalence stratum,
dropped and retained peaks face the same sampling-noise floor, so a specificity
difference there is real. Reported as:

  * per-stratum medians and rank-biserial effect
  * a pooled effect (weighted by stratum size)
  * the null-expected effect from a permutation of the labels WITHIN strata,
    which is the honest reference

This reads the per-peak CSV that characterize_dropped.py already wrote, so it
needs no second pass over the matrix.

Usage
-----
    python 04_db/reanalyze_dropped.py \\
        --per-peak dropped_characterization.per_peak.csv \\
        --out      dropped_stratified
"""
import argparse
import json
import sys
from pathlib import Path

try:
    import numpy as np
    import pandas as pd
    from scipy.stats import mannwhitneyu
except ImportError as e:                                         # pragma: no cover
    sys.exit(f"ERROR: needs numpy, pandas, scipy ({e})")


def rank_biserial(a, b):
    """>0 means b (dropped) scores higher than a (retained)."""
    a, b = a[np.isfinite(a)], b[np.isfinite(b)]
    if len(a) < 5 or len(b) < 5:
        return np.nan, np.nan
    u = mannwhitneyu(b, a, alternative="two-sided")
    return 2 * u.statistic / (len(a) * len(b)) - 1, float(u.pvalue)


def match_on_prevalence(prev, spec, rep, tol=0.02, window=40, seed=0):
    """1:1 nearest-neighbour matching on prevalence, without replacement.

    Quantile stratification is not enough: even inside a decile the dropped
    peaks sit at the low end and the retained at the high end, because the
    marginal prevalence distributions differ. On a fixture where the ONLY
    difference was prevalence, ten strata left a pooled effect of +0.047 that
    still cleared the permutation null -- a false positive. Exact matching
    reduces the same fixture to +0.009 (p=0.24) while leaving a genuine
    specificity difference at +0.673.

    Each dropped peak is paired with the nearest unused retained peak whose
    prevalence is within `tol` relatively; unmatched peaks are discarded.
    """
    r = np.random.default_rng(seed)
    ri, di = np.flatnonzero(rep), np.flatnonzero(~rep)
    order = np.argsort(prev[ri])
    ri_s, pv_s = ri[order], prev[ri][order]
    used = np.zeros(len(ri_s), dtype=bool)
    pairs = []
    for j in r.permutation(di):            # random order avoids greedy bias
        k = np.searchsorted(pv_s, prev[j])
        best, bestd = -1, np.inf
        for c in range(max(0, k - window), min(len(pv_s), k + window)):
            if used[c]:
                continue
            d = abs(pv_s[c] - prev[j])
            if d < bestd:
                best, bestd = c, d
        if best < 0 or (prev[j] > 0 and bestd / prev[j] > tol):
            continue
        used[best] = True
        pairs.append((ri_s[best], j))
    if len(pairs) < 200:
        return None
    ia = [a for a, _ in pairs]
    ib = [b for _, b in pairs]
    e, p = rank_biserial(spec[ia], spec[ib])
    return dict(n_pairs=len(pairs), effect=float(e), p=float(p),
                prevalence_ratio=float(np.median(prev[ib])
                                       / np.median(prev[ia])),
                median_retained=float(np.nanmedian(spec[ia])),
                median_dropped=float(np.nanmedian(spec[ib])))


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--per-peak", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--n-strata", type=int, default=10,
                    help="prevalence quantile bins (default 10)")
    ap.add_argument("--min-per-cell", type=int, default=50,
                    help="skip strata with fewer than this in either class")
    ap.add_argument("--no-plot", action="store_true")
    args = ap.parse_args()

    if not args.per_peak.exists():
        sys.exit(f"ERROR: {args.per_peak} not found")
    df = pd.read_csv(args.per_peak)
    need = {"representable", "prevalence", "specificity"}
    if not need <= set(df.columns):
        sys.exit(f"ERROR: need columns {sorted(need)}; have {list(df.columns)}")
    rep = df["representable"].astype(bool).to_numpy()
    prev = df["prevalence"].to_numpy()
    spec = df["specificity"].to_numpy()

    print(f"[input] {len(df):,} peaks | retained {rep.sum():,} | "
          f"dropped {(~rep).sum():,}")
    e_naive, p_naive = rank_biserial(spec[rep], spec[~rep])
    e_prev, p_prev = rank_biserial(prev[rep], prev[~rep])
    print(f"[naive] specificity effect {e_naive:+.3f} (p={p_naive:.2e})  "
          f"-- CONFOUNDED by prevalence effect {e_prev:+.3f}")

    # Quantile strata on prevalence, pooled across both classes so the bin
    # edges do not depend on which class is being compared.
    qs = np.quantile(prev, np.linspace(0, 1, args.n_strata + 1))
    qs = np.unique(qs)
    if len(qs) < 3:
        sys.exit("ERROR: prevalence has too few distinct values to stratify")
    binid = np.clip(np.digitize(prev, qs[1:-1]), 0, len(qs) - 2)

    rows, eff, wts = [], [], []
    for b in range(len(qs) - 1):
        m = binid == b
        a, d = spec[m & rep], spec[m & ~rep]
        if min(len(a), len(d)) < args.min_per_cell:
            rows.append(dict(stratum=b, prev_lo=float(qs[b]),
                             prev_hi=float(qs[b + 1]),
                             n_retained=int(len(a)), n_dropped=int(len(d)),
                             median_retained=None, median_dropped=None,
                             effect=None, p=None, note="too few"))
            continue
        e, p = rank_biserial(a, d)
        rows.append(dict(stratum=b, prev_lo=float(qs[b]),
                         prev_hi=float(qs[b + 1]),
                         n_retained=int(len(a)), n_dropped=int(len(d)),
                         median_retained=float(np.nanmedian(a)),
                         median_dropped=float(np.nanmedian(d)),
                         effect=float(e), p=float(p), note=""))
        eff.append(e)
        wts.append(min(len(a), len(d)))

    if not eff:
        sys.exit("ERROR: no stratum had enough peaks in both classes")
    eff, wts = np.array(eff), np.array(wts, dtype=float)
    pooled = float(np.average(eff, weights=wts))

    # Null reference: shuffle the class label WITHIN each stratum. This holds
    # prevalence fixed by construction, so it isolates sampling noise.
    rngs = np.random.default_rng(0)
    null = []
    for _ in range(200):
        perm = rep.copy()
        for b in range(len(qs) - 1):
            idx = np.flatnonzero(binid == b)
            if len(idx) > 1:
                perm[idx] = rngs.permutation(perm[idx])
        e_b, w_b = [], []
        for b in range(len(qs) - 1):
            m = binid == b
            a, d = spec[m & perm], spec[m & ~perm]
            if min(len(a), len(d)) < args.min_per_cell:
                continue
            e, _ = rank_biserial(a, d)
            if np.isfinite(e):
                e_b.append(e)
                w_b.append(min(len(a), len(d)))
        if e_b:
            null.append(np.average(e_b, weights=w_b))
    null = np.array(null)
    lo, hi = np.percentile(null, [2.5, 97.5])

    print()
    print("=" * 78)
    print("PREVALENCE-STRATIFIED SPECIFICITY")
    print("=" * 78)
    print(f"{'stratum':>8}{'prevalence range':>26}{'n ret':>9}{'n drop':>9}"
          f"{'med ret':>10}{'med drop':>10}{'effect':>9}")
    for r in rows:
        rng_s = f"{r['prev_lo']:.5f}-{r['prev_hi']:.5f}"
        if r["effect"] is None:
            print(f"{r['stratum']:>8}{rng_s:>26}{r['n_retained']:>9,}"
                  f"{r['n_dropped']:>9,}{'--':>10}{'--':>10}"
                  f"{'skipped':>9}")
            continue
        print(f"{r['stratum']:>8}{rng_s:>26}{r['n_retained']:>9,}"
              f"{r['n_dropped']:>9,}{r['median_retained']:>10.3f}"
              f"{r['median_dropped']:>10.3f}{r['effect']:>+9.3f}")
    print("-" * 78)
    print(f"pooled effect (size-weighted) : {pooled:+.4f}")
    print(f"null 95% interval (label perm): [{lo:+.4f}, {hi:+.4f}]")
    print(f"naive, unstratified           : {e_naive:+.4f}")

    # --- the primary test: exact matching on prevalence --------------------
    mt = match_on_prevalence(prev, spec, rep)
    print()
    if mt is None:
        print("matched test: too few pairs -- falling back to stratification")
    else:
        print("=" * 78)
        print("PREVALENCE-MATCHED SPECIFICITY  (the primary test)")
        print("=" * 78)
        print(f"  matched pairs                 : {mt['n_pairs']:,}")
        print(f"  prevalence ratio, drop/ret    : {mt['prevalence_ratio']:.4f}"
              f"   (1.000 = confound removed)")
        print(f"  median specificity, retained  : {mt['median_retained']:.3f}")
        print(f"  median specificity, dropped   : {mt['median_dropped']:.3f}")
        print(f"  effect (rank-biserial)        : {mt['effect']:+.4f}"
              f"   p = {mt['p']:.2e}")

    drive = mt["effect"] if mt else pooled
    drive_p = mt["p"] if mt else None
    if mt is not None and (drive_p is None or drive_p >= 0.01
                           or abs(drive) < 0.05):
        verdict = "PRECOMPUTED DATABASE IS ADEQUATE"
        why = ("at matched prevalence there is no meaningful specificity "
               "difference. The naive gap was the prevalence confound: dropped "
               "peaks are sparser, and sparser peaks look more specific by "
               "chance. The loss is weak peaks, not restricted enhancers")
    elif drive > 0:
        verdict = "CUSTOM DATABASE JUSTIFIED"
        why = ("at matched prevalence, dropped peaks are still more "
               "cell-type-specific than retained ones -- the database is "
               "discarding restricted elements, not merely sparse ones")
    else:
        verdict = "PRECOMPUTED DATABASE IS FINE"
        why = ("at matched prevalence, dropped peaks are LESS "
               "cell-type-specific than retained ones; the naive difference "
               "was an artifact of their lower prevalence")
    print("-" * 78)
    print(f"VERDICT: {verdict}")
    import textwrap
    for line in textwrap.wrap(why, 74):
        print(f"  {line}")
    print("=" * 78)

    out = dict(n_peaks=int(len(df)), naive_specificity_effect=float(e_naive),
               prevalence_effect=float(e_prev), pooled_stratified_effect=pooled,
               null_ci=[float(lo), float(hi)], matched=mt, strata=rows,
               verdict=verdict, rationale=why)
    Path(f"{args.out}.summary.json").write_text(json.dumps(out, indent=2))
    pd.DataFrame(rows).to_csv(f"{args.out}.strata.csv", index=False)
    print(f"  summary: {args.out}.summary.json")
    print(f"  strata : {args.out}.strata.csv")


if __name__ == "__main__":
    main()
