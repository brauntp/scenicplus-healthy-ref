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

    def skip(self, kind, name, detail=""):
        """A claim that could NOT be checked here -- neither pass nor fail.

        A missing input is not evidence the claim holds, so it must not report
        green; it is also not a defect, so it must not fail the gate. Both were
        previously conflated: an absent catalog scored as a pass (silently
        green) while an absent numpy scored as a failure (crying wolf). Either
        way the analyst learns the wrong thing from the summary line.
        """
        self.rows.append((kind, name, None, detail))

    def report(self):
        bad = [r for r in self.rows if r[2] is False]
        skipped = [r for r in self.rows if r[2] is None]
        ran = len(self.rows) - len(skipped)
        width = max(len(r[1]) for r in self.rows) if self.rows else 10
        for kind, name, ok, detail in self.rows:
            mark = "skip" if ok is None else ("ok  " if ok else "FAIL")
            print(f"  [{mark}] {kind:<11} {name:<{width}}  {detail}")
        print(f"\n{ran - len(bad)}/{ran} checks passed", end="")
        if skipped:
            print(f", {len(skipped)} NOT CHECKED (unverified, not confirmed):")
            for _, name, _, detail in skipped:
                print(f"    - {name}: {detail}")
        else:
            print()
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
        ck.skip("recomputed", "catalog widths", "numpy/pandas unavailable")
        return
    p = Path(path)
    if not p.exists():
        ck.skip("recomputed", "catalog widths", f"{path} absent")
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



# ------------------------------------------------------- README pairing claims
# The README's central quantitative argument -- that GLUE-anchored pairing beats
# stock label-random pairing -- is a table of numbers copied out of
# pairing_sensitivity.csv. Copied numbers drift from their source. These checks
# re-derive every one of them from the CSV. Pure stdlib, so they run in any
# environment, including one without numpy/pandas.

def _sensitivity_rows(csv_path):
    import csv as _csv
    with open(csv_path, newline="") as fh:
        return list(_csv.DictReader(fh))


def readme_pairing_table(text, csv_path, ck):
    """The two table rows must match noise=0.15, k=25 in the CSV exactly."""
    if not Path(csv_path).exists():
        ck.skip("recomputed", "pairing table", f"{csv_path} not found")
        return
    rows = _sensitivity_rows(csv_path)
    for pair, label in (("GLUE", "GLUE-anchored"), ("stock", "SCENIC+ stock")):
        r = [x for x in rows
             if x["pair"] == pair and _num(x["noise"]) == 0.15 and _num(x["k"]) == 25]
        if len(r) != 1:
            ck.check("recomputed", f"{pair} row present", False,
                     f"{len(r)} rows at noise=0.15,k=25")
            continue
        r = r[0]
        # Each stated figure must appear in the README's table row for this pair.
        line = next((l for l in text.splitlines()
                     if l.startswith("|") and label in l), "")
        ok = bool(line)
        detail = []
        for field, fmt in (("median_rho_true", None),
                           ("median_rho_decoy", None),
                           ("auroc", None)):
            v = _num(r[field])
            # README renders -0.001 with a unicode minus and bolds some cells;
            # normalise both sides to a bare number string before comparing.
            norm = line.replace("\u2212", "-").replace("*", "")
            cand = {f"{v:g}", f"{v:.3f}", f"{v:.2f}", f"{v:.1f}"}
            hit = any(c in norm for c in cand)
            ok = ok and hit
            detail.append(f"{field}={v:g}{'' if hit else ' MISSING'}")
        ck.check("recomputed", f"README row: {pair}", ok, ", ".join(detail))

    # n_metacells = 320 for both rows, as the README's caption states.
    ns = {r["n_metacells"] for r in rows
          if _num(r["noise"]) == 0.15 and _num(r["k"]) == 25}
    ck.check("recomputed", "README n_metacells",
             ns == {"320"} and "320 metacells" in text,
             f"csv={sorted(ns)}, README says 320")


def readme_auroc_degradation(text, csv_path, ck):
    """'AUROC >= 0.987 at noise <= 0.6, 0.89 at 1.0, 0.63 at 2.5 (best k)'."""
    if not Path(csv_path).exists():
        ck.skip("recomputed", "auroc degradation", f"{csv_path} not found")
        return
    rows = [r for r in _sensitivity_rows(csv_path) if r["pair"] == "GLUE"]
    best = {}
    for r in rows:
        n = _num(r["noise"])
        best[n] = max(best.get(n, 0.0), _num(r["auroc"]))

    m = re.search(r"AUROC\s*[\u2265>]=?\s*([\d.]+)\s*\n?\s*while.*?"
                  r"at or below\s*([\d.]+)", text, re.S)
    if m:
        thr, upto = _num(m.group(1)), _num(m.group(2))
        low = [v for n, v in best.items() if n <= upto]
        ck.check("recomputed", "auroc floor below noise",
                 low and min(low) >= thr,
                 f"min best-AUROC at noise<={upto:g} is {min(low):.3f} >= {thr}")
    else:
        ck.check("consistency", "auroc floor below noise", False,
                 "claim sentence not found in README")

    for noise, pat in ((1.0, r"(?:falls to|0\.89 at)\s*([\d.]+)?\s*at noise 1\.0"),
                       (2.5, r"reaches\s*([\d.]+)\s*at 2\.5")):
        mm = re.search(pat, text)
        stated = _num(mm.group(1)) if (mm and mm.group(1)) else None
        if stated is None and noise == 1.0:
            mm2 = re.search(r"falls to\s*([\d.]+)", text)
            stated = _num(mm2.group(1)) if mm2 else None
        if stated is None:
            ck.check("consistency", f"auroc at noise {noise:g}", False, "not stated")
            continue
        actual = best.get(noise)
        ok = actual is not None and abs(round(actual, 2) - stated) < 0.011
        ck.check("recomputed", f"auroc at noise {noise:g}", ok,
                 f"csv best {actual:.3f} vs README {stated}")


def readme_k_choice(text, csv_path, ck):
    """'k=50 has the highest median rho at every noise level' + 798/320/160."""
    if not Path(csv_path).exists():
        ck.skip("recomputed", "k=50 claim", f"{csv_path} not found")
        return
    rows = [r for r in _sensitivity_rows(csv_path) if r["pair"] == "GLUE"]
    by_noise = {}
    for r in rows:
        by_noise.setdefault(_num(r["noise"]), []).append(
            (_num(r["median_rho_true"]), int(_num(r["k"]))))
    losers = [n for n, v in by_noise.items() if max(v)[1] != 50]
    ck.check("recomputed", "k=50 best at every noise", not losers,
             "all noise levels" if not losers else f"fails at {sorted(losers)}")

    counts = {int(_num(r["k"])): r["n_metacells"] for r in rows}
    stated = re.search(r"k=10/25/50\s*[\u2192-]+\s*(\d+)/(\d+)/(\d+)", text)
    if stated:
        want = [stated.group(i) for i in (1, 2, 3)]
        got = [counts.get(k) for k in (10, 25, 50)]
        ck.check("recomputed", "metacell counts", want == got,
                 f"README {'/'.join(want)} vs csv {'/'.join(map(str, got))}")
    else:
        ck.check("consistency", "metacell counts", False, "claim not found")


# --------------------------------------------------- Snakefile structural facts
# This repo makes load-bearing claims about a file it does not contain. Fetch it
# with docs/fetch_snakefile.sh, then these re-derive each claim from the source.

def snakefile_claims(sf_path, ck):
    sf = Path(sf_path)
    if not sf.exists():
        ck.skip("recomputed", "snakefile claims",
                "not fetched -- run: bash docs/fetch_snakefile.sh")
        return
    src = sf.read_text()

    ck.check("recomputed", "no conda: directives", src.count("conda:") == 0,
             f"{src.count('conda:')} found (env must be pre-activated)")

    keys = set(re.findall(r'config\["(\w+)"\]\["(\w+)"\]', src))
    lookups = len(re.findall(r'config\["\w+"\]\["\w+"\]', src))
    ck.check("recomputed", "config surface",
             len(keys) == 67 and lookups == 151,
             f"{len(keys)} keys over {lookups} lookups")

    # The wrapper's REQUIRED_KEYS must equal the Snakefile's key set exactly:
    # a missing key means a config passes validation and dies inside snakemake.
    rp = (ROOT / "03_pipeline" / "run_pipeline.sh").read_text()
    block = re.search(r"REQUIRED_KEYS=\((.*?)\n\)", rp, re.S)
    if block:
        wrapper = {tuple(l.strip().split("."))
                   for l in block.group(1).splitlines() if l.strip()}
        ck.check("recomputed", "wrapper key list exact",
                 wrapper == keys,
                 f"missing={sorted(keys - wrapper)} extra={sorted(wrapper - keys)}")
    else:
        ck.check("consistency", "wrapper key list exact", False,
                 "REQUIRED_KEYS block not found in run_pipeline.sh")

    # DEM is gated by a conditional, but BOTH branches declare its outputs, so
    # it always runs. This is why --skip-scores breaks the DAG.
    #
    # Count the HTML output, not dem_result_fname: the latter legitimately
    # appears a third time as prepare_menr's INPUT, which an earlier version of
    # this check miscounted as a third branch.
    dem_html = len(re.findall(r'output_fname_dem_html=config\["output_data"\]', src))
    ck.check("recomputed", "DEM outputs unconditional", dem_html == 2,
             f"declared as output in {dem_html}/2 conditional branches")

    # ...and DEM's result is consumed downstream, so it gates the whole DAG:
    # 2 output declarations + exactly 1 input declaration (prepare_menr).
    dem_all = len(re.findall(r'dem_result_fname=config\["output_data"\]', src))
    ck.check("recomputed", "DEM gates prepare_menr", dem_all == dem_html + 1,
             f"{dem_all} refs = {dem_html} outputs + {dem_all - dem_html} input")

    ctx = len(re.findall(r'=config\["input_data"\]\["ctx_db_fname"\]', src))
    ck.check("recomputed", "ctx_db read by 3 rules", ctx == 3,
             f"{ctx} rules (cistarget + both eGRN)")


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--doc", type=Path,
                    default=ROOT / "docs" / "DATABASE_DECISION.md")
    ap.add_argument("--regions", default="screen_db_regions.csv",
                    help="region catalog, for the recomputed checks "
                         "(skipped if absent)")
    ap.add_argument("--readme", type=Path, default=ROOT / "README.md")
    ap.add_argument("--sensitivity", type=Path,
                    default=ROOT / "docs" / "pairing_sensitivity.csv")
    ap.add_argument("--snakefile", type=Path,
                    default=ROOT / "docs" / "reference" / "Snakefile.v1.0a2",
                    help="fetch with docs/fetch_snakefile.sh")
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

    if args.readme.exists():
        print(f"\nverifying {args.readme}\n")
        rtext = args.readme.read_text()
        readme_pairing_table(rtext, args.sensitivity, ck)
        readme_auroc_degradation(rtext, args.sensitivity, ck)
        readme_k_choice(rtext, args.sensitivity, ck)

    print(f"\nverifying claims about the pinned Snakefile\n")
    snakefile_claims(args.snakefile, ck)

    sys.exit(ck.report())


if __name__ == "__main__":
    main()
