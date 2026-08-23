#!/usr/bin/env python3
"""
Recompute every quantitative claim in the docs from its source, and check the
claims are internally consistent.

WHY THIS EXISTS
---------------
A substring check is not verification. DATABASE_DECISION.md once asserted that
SCREEN regions are "150-350 bp (median 272), all narrower than 0.4 x 500 =
200 bp" -- arithmetically impossible, since a median of 272 cannot come from a
population entirely below 200. My check at the time confirmed that the strings
"150", "272" and "0.4 x 500 = 200" each appeared in the file, which they did.
Every individual number was right and the sentence was still false, because the
check never asked whether the claims agreed with each other.

Two classes of check, and the second is the one that was missing:

  RECOMPUTED  -- the claim is re-derived from the source artifact (a CSV, the
                 region catalog) and compared to what the doc says.
  CONSISTENCY -- the claim is checked against OTHER claims in the same document.
                 A stated median must lie inside a stated range; parts must sum
                 to their stated whole; a percentage must match its own
                 numerator and denominator.

Usage
-----
    python docs/verify_claims.py                      # all docs it knows about
    python docs/verify_claims.py --doc docs/DATABASE_DECISION.md
    python docs/verify_claims.py --regions screen_db_regions.csv   # + recompute

Exit status is nonzero if any check fails, so this can gate a commit.
"""
import argparse
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent


def _num(s):
    return float(str(s).replace(",", "").replace("\u2212", "-").strip())


class Checker:
    def __init__(self):
        self.rows = []

    def check(self, kind, name, ok, detail=""):
        self.rows.append((kind, name, bool(ok), detail))

    def report(self):
        bad = [r for r in self.rows if not r[2]]
        width = max(len(r[1]) for r in self.rows) if self.rows else 10
        for kind, name, ok, detail in self.rows:
            mark = "ok  " if ok else "FAIL"
            print(f"  [{mark}] {kind:<11} {name:<{width}}  {detail}")
        print(f"\n{len(self.rows) - len(bad)}/{len(self.rows)} checks passed")
        return 1 if bad else 0


# ---------------------------------------------------------------- consistency

def median_inside_range(text, ck):
    """Any 'A-B bp (median M)' must satisfy A <= M <= B.

    This is the check that would have caught the width claim.
    """
    pat = re.compile(r"(\d[\d,]*)\s*[-\u2013]\s*(\d[\d,]*)\s*bp\s*"
                     r"\(median\s+(\d[\d,]*)\)")
    found = 0
    for m in pat.finditer(text):
        lo, hi, med = (_num(m.group(i)) for i in (1, 2, 3))
        found += 1
        ck.check("consistency", f"median in range {lo:.0f}-{hi:.0f}",
                 lo <= med <= hi, f"median {med:.0f}")
    if not found:
        ck.check("consistency", "range/median patterns", True, "none present")


def universal_vs_median(text, ck):
    """Catch 'all X are <T>' / 'all narrower than T' when a stated median >= T.

    The exact failure mode from DATABASE_DECISION.md: a universal quantifier
    about a population whose own reported median contradicts it.
    """
    meds = [_num(m) for m in re.findall(r"median\s+(\d[\d,]*)", text)]
    # Note the character class: an earlier version used [^.\n], which silently
    # failed on the very sentence this check exists for, because "0.4 x 500"
    # contains a period. Allow any non-newline character and let the non-greedy
    # quantifier find the number that precedes the unit.
    hits = list(re.finditer(
        r"\ball\b[^\n]{0,60}?(?:narrower|smaller|less|below|under|<=|\u2264)"
        r"[^\n]{0,50}?(\d[\d,]*)\s*bp", text, re.I))
    if not hits:
        ck.check("consistency", "universal claims", True,
                 "no 'all ... below N bp' claims")
        return
    for m in hits:
        thresh = _num(m.group(1))
        conflict = [x for x in meds if x >= thresh]
        ck.check("consistency", f"'all < {thresh:.0f} bp' vs medians",
                 not conflict,
                 f"conflicting medians: {conflict}" if conflict
                 else "no stated median contradicts it")


def parts_sum_to_whole(text, ck):
    """The dropped-peak decomposition must add up, and match the total."""
    def grab(pat):
        m = re.search(pat, text)
        return _num(m.group(1)) if m else None

    total = grab(r"silently dropped\s*\|\s*\*\*([\d,]+)")
    zero = grab(r"zero overlap with any cCRE\s*\|\s*([\d,]+)")
    sub = grab(r"overlaps but below threshold\s*\|\s*([\d,]+)")
    rep = grab(r"representable in DB\s*\|\s*\*\*([\d,]+)")
    audited = grab(r"(\d[\d,]*)\s+consensus peaks")
    if None in (total, zero, sub):
        ck.check("consistency", "dropped parts sum", False,
                 "could not parse the decomposition")
    else:
        ck.check("consistency", "dropped parts sum", zero + sub == total,
                 f"{zero:.0f}+{sub:.0f}={zero+sub:.0f} vs {total:.0f}")
    if None not in (rep, total, audited):
        ck.check("consistency", "rep + dropped = audited",
                 rep + total == audited,
                 f"{rep:.0f}+{total:.0f}={rep+total:.0f} vs {audited:.0f}")


def percentages_match(text, ck):
    """'N (P%)' and 'N (**P% of the loss**)' must be arithmetically right."""
    total = re.search(r"silently dropped\s*\|\s*\*\*([\d,]+)", text)
    total = _num(total.group(1)) if total else None
    for label, pat in (("82% of loss",
                        r"zero overlap with any cCRE\s*\|\s*([\d,]+)\s*"
                        r"\(\*\*(\d+)% of the loss\*\*\)"),
                       ("18% of loss",
                        r"overlaps but below threshold\s*\|\s*([\d,]+)\s*"
                        r"\((\d+)% of the loss\)")):
        m = re.search(pat, text)
        if not m or total is None:
            ck.check("consistency", label, False, "not parsed")
            continue
        n, pct = _num(m.group(1)), _num(m.group(2))
        ck.check("consistency", label, abs(n / total * 100 - pct) < 0.6,
                 f"{n:.0f}/{total:.0f} = {n/total*100:.1f}% vs stated {pct:.0f}%")


def overlap_routes(text, ck):
    """The two overlap requirements must follow from the widths and 0.4."""
    m = re.search(r"0\.4\s*[x\u00d7]\s*(\d[\d,]*)\s*[-\u2013]\s*(\d[\d,]*)\s*=\s*"
                  r"\*\*(\d[\d,]*)\s*[-\u2013]\s*(\d[\d,]*)\s*bp\*\*", text)
    if m:
        lo, hi, rlo, rhi = (_num(m.group(i)) for i in (1, 2, 3, 4))
        import math
        ck.check("consistency", "query-route arithmetic",
                 math.ceil(0.4 * lo) == rlo and math.ceil(0.4 * hi) == rhi,
                 f"0.4x({lo:.0f},{hi:.0f}) -> ({rlo:.0f},{rhi:.0f})")
    else:
        ck.check("consistency", "query-route arithmetic", False, "not parsed")
    m = re.search(r"0\.4\s*[x\u00d7]\s*(\d[\d,]*)\s*=\s*\*\*(\d[\d,]*)\s*bp\*\*",
                  text)
    if m:
        p, r = _num(m.group(1)), _num(m.group(2))
        ck.check("consistency", "target-route arithmetic", 0.4 * p == r,
                 f"0.4x{p:.0f}={0.4*p:.0f} vs {r:.0f}")
    else:
        ck.check("consistency", "target-route arithmetic", False, "not parsed")


# ---------------------------------------------------------------- recomputed

def recompute_from_catalog(text, path, ck):
    """Re-derive the width claims from the region catalog itself."""
    try:
        import numpy as np
        import pandas as pd
    except ImportError:
        ck.check("recomputed", "catalog", False, "numpy/pandas unavailable")
        return
    p = Path(path)
    if not p.exists():
        ck.check("recomputed", "catalog", True, f"{path} absent -- skipped")
        return
    df = (pd.read_parquet(p) if p.suffix == ".parquet"
          else pd.read_csv(p))
    w = (df["end"] - df["start"]).to_numpy()
    stated = re.search(r"(\d[\d,]*)\s*[-\u2013]\s*(\d[\d,]*)\s*bp\s*"
                       r"\(median\s+(\d[\d,]*)\)", text)
    if stated:
        lo, hi, med = (_num(stated.group(i)) for i in (1, 2, 3))
        ck.check("recomputed", "width min", w.min() == lo, f"{w.min()} vs {lo:.0f}")
        ck.check("recomputed", "width max", w.max() == hi, f"{w.max()} vs {hi:.0f}")
        ck.check("recomputed", "width median",
                 int(np.median(w)) == med, f"{int(np.median(w))} vs {med:.0f}")
    m = re.search(r"median\s+(\d+)\s*bp\s+of\s+overlap|median\s+(\d+)\s*bp,\s*not",
                  text)
    need_q = np.ceil(0.4 * w)
    if m:
        stated_med = _num(next(g for g in m.groups() if g))
        ck.check("recomputed", "median required overlap",
                 int(np.median(need_q)) == stated_med,
                 f"{int(np.median(need_q))} vs {stated_med:.0f}")
    m = re.search(r"\*\*(\d+)%\*\*\s*of regions", text)
    if m:
        stated_pct = _num(m.group(1))
        actual = (need_q < 0.4 * 500).mean() * 100
        ck.check("recomputed", "query route easier for",
                 abs(actual - stated_pct) < 0.5,
                 f"{actual:.1f}% vs {stated_pct:.0f}%")
    m = re.search(r"([\d.]+)% of SCREEN regions are themselves\s*\u2264\s*(\d+)\s*bp",
                  text)
    if m:
        pct, thr = _num(m.group(1)), _num(m.group(2))
        actual = (w <= thr).mean() * 100
        ck.check("recomputed", f"fraction <= {thr:.0f} bp",
                 abs(actual - pct) < 0.15, f"{actual:.1f}% vs {pct:.1f}%")
    ck.check("recomputed", "catalog rows",
             f"{len(df):,}" in text, f"{len(df):,}")


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--doc", type=Path,
                    default=ROOT / "docs" / "DATABASE_DECISION.md")
    ap.add_argument("--regions", default="screen_db_regions.csv",
                    help="region catalog, for the recomputed checks "
                         "(skipped if absent)")
    args = ap.parse_args()

    if not args.doc.exists():
        sys.exit(f"ERROR: {args.doc} not found")
    text = args.doc.read_text()
    ck = Checker()
    print(f"verifying {args.doc}\n")
    median_inside_range(text, ck)
    universal_vs_median(text, ck)
    parts_sum_to_whole(text, ck)
    percentages_match(text, ck)
    overlap_routes(text, ck)
    recompute_from_catalog(text, args.regions, ck)
    sys.exit(ck.report())


if __name__ == "__main__":
    main()
