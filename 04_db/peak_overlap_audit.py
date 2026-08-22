#!/usr/bin/env python3
"""
peak_overlap_audit.py -- how much of your ATAC peak set is representable in a
cisTarget database?

WHY THIS SCRIPT EXISTS
----------------------
SCENIC+ / pycistarget do NOT run motif enrichment on your peaks directly. They
map your peaks onto the regions that exist in the cisTarget database, and then
score the *database* regions. Any peak of yours that fails to map is dropped --
silently, with no warning, no log line, and no entry in the results. If you hand
pycistarget a peak set whose regions are only 60% representable in the
precomputed database, you lose 40% of your regulatory landscape and nothing in
the output tells you so.

This script computes that loss number before you commit, so you can decide
whether to build a custom database (expensive) or use the precomputed one
(free).

THE RULE, TAKEN FROM THE SOURCE (NOT FROM MEMORY)
-------------------------------------------------
pycistarget/src/pycistarget/utils.py :: target_to_query(), the function
cisTargetDatabase.load_db() calls to build its peak->DB mapping:

    join_pr = target_pr.join(query_pr, report_overlap=True)
    join_pr.Overlap_query  = join_pr.Overlap / (join_pr.End_b - join_pr.Start_b)
    join_pr.Overlap_target = join_pr.Overlap / (join_pr.End   - join_pr.Start)
    join_pr = join_pr[(join_pr.Overlap_query  > fraction_overlap) |
                      (join_pr.Overlap_target > fraction_overlap)]

Three details that matter and that a from-memory implementation gets wrong:

  1. TARGET is your peak set, QUERY is the database.
     In cisTargetDatabase.load_db():
         target_to_query(region_sets[x], list(db_regions), fraction_overlap=...)
     i.e. target=your regions, query=db_regions. So:
         Overlap_target = overlap / len(your peak)
         Overlap_query  = overlap / len(DB region)

  2. The test is a reciprocal OR, not an AND, and not one-sided.
     A peak is kept if EITHER fraction clears the threshold. This is much more
     permissive than an AND and it is the reason small DB regions still work:
     a 150 bp DB region sitting entirely inside your 501 bp peak gives
     Overlap_target = 150/501 = 0.299 (FAILS 0.4) but Overlap_query = 1.000
     (PASSES). Implementing this as an AND, or as peak-side-only, will
     dramatically over-report your loss.

  3. The comparison is STRICTLY GREATER THAN (`>`), not `>=`.
     At the default 0.4 an overlap fraction of exactly 0.4 is REJECTED. This
     script reproduces `>` exactly.

Coordinates: BED is 0-based half-open. pycistarget's region_names_to_coordinates()
splits "chr:start-end" and uses start/end as PyRanges Start/End with no
adjustment, so DB region IDs are also treated as 0-based half-open and the two
are directly comparable. Overlap = max(0, min(end1,end2) - max(start1,start2)).

ALGORITHM
---------
Exact interval intersection via sort + numpy.searchsorted, per chromosome:

  * DB regions for a chromosome are sorted by start; `starts`/`ends` are numpy
    arrays. W = max(ends - starts) is the widest DB region on that chromosome.
  * For a query peak [ps, pe), all DB regions that can possibly overlap satisfy
    start < pe  AND  end > ps. Since end <= start + W, any region with
    start <= ps - W has end <= ps and cannot overlap. Therefore the candidate
    window is exactly
        lo = searchsorted(starts, ps - W, side="left")
        hi = searchsorted(starts, pe,     side="left")
    and scanning [lo, hi) is guaranteed to miss no overlap. This is an exact
    bound, not a heuristic -- no interval is skipped.
  * Within the window, overlaps are computed vectorised with numpy.

This is O(n log n) for the sort plus O(k) per peak where k is the local DB
density in a W-wide neighbourhood. For the hg38 SCREEN database (1.84M regions,
max width 350 bp) k is small, and 300k peaks audit in well under a minute with
no pyranges dependency.

Correctness is not taken on faith: --self-test runs this implementation against
an O(n*m) brute-force on randomised interval sets and aborts on any mismatch.
If pyranges is installed, --self-test additionally cross-checks the pass/fail
verdict against a real pycistarget-style pyranges join.

USAGE
-----
  ./peak_overlap_audit.py \
      --peaks consensus_peaks.bed \
      --db    hg38_screen_v10_clust.regions_vs_motifs.rankings.feather \
      --out-prefix audit/screen_vs_mypeaks

  ./peak_overlap_audit.py --self-test        # verify the algorithm, no inputs

Only the feather SCHEMA is read (plus, for motifs_vs_regions layouts, a single
column). The 35 GB rankings file is never loaded into memory.

Outputs:
  <out-prefix>.per_peak.csv   one row per input peak, best overlap + verdict
  <out-prefix>.summary.json   aggregate stats, per-chromosome and per-type loss
  stdout                      the verdict: % of peaks SILENTLY DROPPED
"""

from __future__ import annotations

import argparse
import json
import os
import re
import sys
from typing import Dict, List, Optional, Tuple

import numpy as np

# --------------------------------------------------------------------------- #
# Region ID parsing
# --------------------------------------------------------------------------- #

# Accept "chr1:100-200" and also "chr1:100-200,+" / trailing junk is rejected.
# Deliberately permissive on the contig token (scaffolds, no-"chr" builds, etc.)
# so we can *detect* naming mismatches rather than crash on them.
_REGION_RE = re.compile(r"^(?P<chrom>[^:]+):(?P<start>\d+)-(?P<end>\d+)$")


def parse_region_id(region_id: str) -> Optional[Tuple[str, int, int]]:
    """Parse 'chr:start-end' -> (chrom, start, end), or None if unparseable."""
    m = _REGION_RE.match(region_id.strip())
    if not m:
        return None
    start, end = int(m.group("start")), int(m.group("end"))
    if end <= start:
        return None
    return m.group("chrom"), start, end


# --------------------------------------------------------------------------- #
# Reading DB region IDs out of a cisTarget feather
# --------------------------------------------------------------------------- #

def read_db_region_ids(feather_path: str, verbose: bool = True) -> List[str]:
    """
    Extract the region IDs from a cisTarget feather database.

    The layout is NOT the same for every cisTarget feather, so we inspect the
    schema and dispatch. From create_cisTarget_databases/cistarget_db.py:

        DatabaseTypes.create_db_filename() ->
            f"{db_prefix}.{column_kind}_vs_{row_kind}.{scores_or_rankings}.feather"

        CisTargetDatabase.write_db() ->
            self.df[self.db_type.row_kind] = self.df.index.to_series()

    so the *row* kind becomes an extra named column holding what were the
    dataframe's row labels, and the *column* kind supplies all the other column
    names. Two region-based layouts therefore exist:

      regions_vs_motifs.*.feather  (column_kind=regions, row_kind=motifs)
          -> schema = [<region IDs> ..., "motifs"]
          -> region IDs are the COLUMN NAMES.  Schema read only; no data read.
             This is what pycistarget consumes and what SCENIC+ expects.

      motifs_vs_regions.*.feather  (column_kind=motifs, row_kind=regions)
          -> schema = [<motif IDs> ..., "regions"]
          -> region IDs are the VALUES of the single "regions" column.
             We read that one column only.

    Verified empirically against the precomputed hg38 SCREEN database
    (hg38_screen_v10_clust.regions_vs_motifs.scores.feather): the Arrow schema
    carries 1,837,305 fields = 1,837,304 "chr:start-end" names + one field
    literally named "motifs", over 5,876 record-batch rows.

    pycistarget also tolerates a "prefix__region" form and strips the prefix
    (motif_enrichment_cistarget.py: `if '__' in db_regions[0]`). We mirror that.
    """
    try:
        import pyarrow.feather as feather
        import pyarrow.ipc as ipc
    except ImportError:
        sys.exit(
            "ERROR: pyarrow is required to read the cisTarget feather schema.\n"
            "       conda install -c conda-forge 'pyarrow>=7.0.0'\n"
            "       (pyarrow>=7.0.0 is what create_cisTarget_databases itself pins.)"
        )

    if not os.path.exists(feather_path):
        sys.exit(f"ERROR: cisTarget database not found: {feather_path}")

    # Read the schema WITHOUT reading the body. Feather v2 == Arrow IPC file.
    try:
        with open(feather_path, "rb") as fh:
            schema = ipc.open_file(fh).schema
        names = list(schema.names)
    except Exception as exc:  # noqa: BLE001 - fall back to the feather reader
        try:
            schema = feather.read_table(feather_path, columns=[]).schema
            names = list(schema.names)
        except Exception:
            sys.exit(
                f"ERROR: could not read an Arrow/Feather schema from {feather_path}\n"
                f"       underlying error: {exc}\n"
                "       Is this really a cisTarget .feather database? Feather v1\n"
                "       databases must be converted first (see\n"
                "       convert_cistarget_databases_v1_to_v2.py in\n"
                "       create_cisTarget_databases)."
            )

    if not names:
        sys.exit(f"ERROR: {feather_path} has an empty schema.")

    if verbose:
        print(f"[db] {os.path.basename(feather_path)}")
        print(f"[db] schema fields: {len(names):,}")

    # -- dispatch on layout ------------------------------------------------- #
    if "regions" in names:
        # motifs_vs_regions layout: region IDs live in the "regions" column.
        if verbose:
            print("[db] layout: *_vs_regions  (region IDs are VALUES of the "
                  "'regions' column; reading that one column)")
        col = feather.read_table(feather_path, columns=["regions"])["regions"]
        region_ids = [str(v) for v in col.to_pylist()]

    elif "motifs" in names or "tracks" in names:
        # regions_vs_motifs / regions_vs_tracks layout: regions are column names.
        idx_col = "motifs" if "motifs" in names else "tracks"
        if verbose:
            print(f"[db] layout: regions_vs_{idx_col}  (region IDs are COLUMN "
                  f"NAMES; index column '{idx_col}' excluded)")
        region_ids = [n for n in names if n != idx_col]

    elif "genes" in names:
        sys.exit(
            f"ERROR: {feather_path} is a GENE-based cisTarget database "
            "(schema contains a 'genes' column).\n"
            "       An ATAC peak set cannot be audited against a gene-based\n"
            "       database. You need a *region-based* database, e.g.\n"
            "       hg38_screen_v10_clust.regions_vs_motifs.rankings.feather"
        )

    else:
        # No recognised index column. Guess from the shape of the names.
        parseable = sum(1 for n in names[:1000] if parse_region_id(n))
        if parseable < 500:
            sys.exit(
                f"ERROR: cannot identify the region IDs in {feather_path}.\n"
                f"       No 'motifs'/'tracks'/'regions'/'genes' index column in the\n"
                f"       schema, and the column names do not look like\n"
                f"       'chr:start-end' ({parseable}/1000 parsed).\n"
                f"       First few fields: {names[:5]}"
            )
        if verbose:
            print("[db] layout: no index column found; treating all column "
                  "names as region IDs")
        region_ids = list(names)

    # -- pycistarget's "prefix__region" handling ----------------------------- #
    if region_ids and "__" in region_ids[0]:
        prefix = region_ids[0].split("__")[0]
        if verbose:
            print(f"[db] stripping pycistarget-style prefix '{prefix}__' from "
                  "region IDs")
        region_ids = [r.split("__", 1)[1] if "__" in r else r for r in region_ids]

    if verbose:
        print(f"[db] region IDs: {len(region_ids):,}")
    return region_ids


# --------------------------------------------------------------------------- #
# BED reading
# --------------------------------------------------------------------------- #

def read_peaks_bed(bed_path: str, verbose: bool = True):
    """
    Read a BED file of consensus peaks.

    Returns (chroms, starts, ends, names, scores, n_skipped) where names/scores
    are None when the BED has fewer than 4/5 columns. Handles track/browser/#
    header lines and gzip.
    """
    if not os.path.exists(bed_path):
        sys.exit(f"ERROR: peak BED not found: {bed_path}")

    opener = open
    if bed_path.endswith(".gz"):
        import gzip
        opener = gzip.open

    chroms: List[str] = []
    starts: List[int] = []
    ends: List[int] = []
    names: List[str] = []
    scores: List[str] = []
    n_cols_seen = 0
    n_skipped = 0

    with opener(bed_path, "rt") as fh:
        for lineno, line in enumerate(fh, 1):
            if not line.strip():
                continue
            if line.startswith(("#", "track", "browser")):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 3:
                f = line.split()
            if len(f) < 3:
                n_skipped += 1
                continue
            try:
                s, e = int(f[1]), int(f[2])
            except ValueError:
                n_skipped += 1
                if n_skipped <= 3:
                    print(f"WARNING: line {lineno}: non-integer coordinates, "
                          f"skipping: {line[:80]!r}", file=sys.stderr)
                continue
            if e <= s:
                n_skipped += 1
                continue
            n_cols_seen = max(n_cols_seen, len(f))
            chroms.append(f[0])
            starts.append(s)
            ends.append(e)
            names.append(f[3] if len(f) > 3 else "")
            scores.append(f[4] if len(f) > 4 else "")

    if not chroms:
        sys.exit(f"ERROR: no usable BED records parsed from {bed_path}")

    has_name = n_cols_seen >= 4 and any(n not in ("", ".") for n in names)
    has_score = n_cols_seen >= 5 and any(s not in ("", ".") for s in scores)

    if verbose:
        w = np.asarray(ends, dtype=np.int64) - np.asarray(starts, dtype=np.int64)
        uniq_w = np.unique(w)
        width_desc = (f"fixed-width {int(uniq_w[0])} bp" if uniq_w.size == 1
                      else f"variable width: min {w.min()} median "
                           f"{int(np.median(w))} max {w.max()}")
        print(f"[peaks] {os.path.basename(bed_path)}")
        print(f"[peaks] {len(chroms):,} peaks, {width_desc}")
        print(f"[peaks] columns detected: {n_cols_seen} "
              f"(name={'yes' if has_name else 'no'}, "
              f"score={'yes' if has_score else 'no'})")
        if n_skipped:
            print(f"[peaks] skipped {n_skipped:,} malformed/zero-length lines")

    return (chroms, starts, ends,
            names if has_name else None,
            scores if has_score else None,
            n_skipped)


# --------------------------------------------------------------------------- #
# The overlap engine
# --------------------------------------------------------------------------- #

class ChromIndex:
    """Sorted interval index for one chromosome, with an exact scan bound."""

    __slots__ = ("starts", "ends", "max_width")

    def __init__(self, starts: np.ndarray, ends: np.ndarray):
        order = np.argsort(starts, kind="stable")
        self.starts = starts[order]
        self.ends = ends[order]
        self.max_width = int((self.ends - self.starts).max()) if starts.size else 0


def read_db_regions_catalog(path: str, verbose: bool = True) -> List[str]:
    """
    Read DB region IDs from a pre-extracted catalog instead of the feather.

    Accepts, by extension:
      .parquet            -- columns chrom,start,end (as written by fetch_db_regions.py)
      .csv / .tsv / .txt  -- either those same columns, or one chr:start-end per line
      .bed                -- chrom, start, end in the first three columns

    Motivation: the precomputed hg38 SCREEN rankings feather is ~33 GiB, but the
    audit only needs its region IDs. fetch_db_regions.py pulls those from the
    Arrow footer over an HTTP range request (~121 MB) and writes this catalog, so
    the representability question can be answered before committing the download.
    """
    if not os.path.exists(path):
        sys.exit(f"ERROR: region catalog not found: {path}")
    ext = os.path.splitext(path)[1].lower()

    def _from_frame(df):
        cols = {c.lower(): c for c in df.columns}
        need = [cols.get(k) for k in ("chrom", "start", "end")]
        if all(need):
            c, s, e = need
            return [f"{a}:{int(b)}-{int(d)}"
                    for a, b, d in zip(df[c], df[s], df[e])]
        return None

    ids = None
    if ext == ".parquet":
        try:
            import pandas as _pd
        except ImportError:
            sys.exit("ERROR: pandas is required to read a .parquet catalog.")
        ids = _from_frame(_pd.read_parquet(path))
        if ids is None:
            sys.exit(f"ERROR: {path} has no chrom/start/end columns.")
    else:
        try:
            import pandas as _pd
            sep = "\t" if ext in (".tsv", ".bed") else None
            df = _pd.read_csv(path, sep=sep, engine="python")
            ids = _from_frame(df)
            if ids is None and ext == ".bed" and df.shape[1] >= 3:
                ids = [f"{r[0]}:{int(r[1])}-{int(r[2])}"
                       for r in df.itertuples(index=False)]
        except Exception:                                     # noqa: BLE001
            ids = None
        if ids is None:                                      # one ID per line
            with open(path) as fh:
                ids = [ln.strip() for ln in fh if ln.strip() and ":" in ln]
    if not ids:
        sys.exit(f"ERROR: no region IDs parsed from {path}")
    if verbose:
        print(f"[db] region catalog {os.path.basename(path)}")
        print(f"[db] {len(ids):,} region IDs (no feather required)")
        print(f"[db] e.g. {', '.join(ids[:3])}")
    return ids


def build_db_index(region_ids: List[str], verbose: bool = True):
    """Group DB region IDs by chromosome into sorted numpy indexes."""
    by_chrom: Dict[str, Tuple[List[int], List[int]]] = {}
    n_bad = 0
    for rid in region_ids:
        p = parse_region_id(rid)
        if p is None:
            n_bad += 1
            continue
        c, s, e = p
        if c not in by_chrom:
            by_chrom[c] = ([], [])
        by_chrom[c][0].append(s)
        by_chrom[c][1].append(e)

    index = {
        c: ChromIndex(np.asarray(ss, dtype=np.int64), np.asarray(ee, dtype=np.int64))
        for c, (ss, ee) in by_chrom.items()
    }
    if verbose:
        n_ok = sum(i.starts.size for i in index.values())
        print(f"[db] indexed {n_ok:,} regions across {len(index):,} contigs")
        if n_bad:
            print(f"[db] WARNING: {n_bad:,} DB region IDs were not parseable as "
                  "'chr:start-end' and were ignored")
        if index:
            mw = max(i.max_width for i in index.values())
            print(f"[db] widest DB region: {mw:,} bp "
                  "(sets the searchsorted scan bound)")
    return index, n_bad


def audit_peaks(peak_chroms: List[str],
                peak_starts: List[int],
                peak_ends: List[int],
                db_index: Dict[str, ChromIndex],
                fraction_overlap: float):
    """
    For every peak, find its best overlap with the DB and apply pycistarget's
    reciprocal-OR, strictly-greater-than rule.

    Returns a dict of numpy arrays, one entry per peak:
      best_ov_bp        largest single-region overlap in bp
      best_frac_target  best overlap/len(peak)      over all DB regions
      best_frac_query   best overlap/len(DB region) over all DB regions
      passes            True if ANY single DB region satisfies
                        (frac_query > X) OR (frac_target > X)
      n_overlapping     number of DB regions with >=1 bp overlap
      n_passing         number of DB regions satisfying the rule
      best_db_start/end the DB region achieving best_ov_bp (-1 if none)
      chrom_in_db       whether the peak's contig exists in the DB at all

    NOTE on `passes`: the rule is evaluated PER CANDIDATE REGION and then OR-ed
    across candidates, exactly as pycistarget does (it filters rows of the join).
    It is NOT max(frac_target) vs max(frac_query) evaluated independently --
    those maxima can come from different DB regions, which would be wrong.
    """
    n = len(peak_chroms)
    best_ov = np.zeros(n, dtype=np.int64)
    best_ft = np.zeros(n, dtype=np.float64)
    best_fq = np.zeros(n, dtype=np.float64)
    passes = np.zeros(n, dtype=bool)
    n_ov = np.zeros(n, dtype=np.int32)
    n_pass = np.zeros(n, dtype=np.int32)
    best_s = np.full(n, -1, dtype=np.int64)
    best_e = np.full(n, -1, dtype=np.int64)
    chrom_ok = np.zeros(n, dtype=bool)

    pc = np.asarray(peak_chroms, dtype=object)
    ps = np.asarray(peak_starts, dtype=np.int64)
    pe = np.asarray(peak_ends, dtype=np.int64)
    plen = (pe - ps).astype(np.float64)

    # Process one chromosome at a time so the searchsorted arrays stay hot.
    for chrom in np.unique(pc):
        rows = np.flatnonzero(pc == chrom)
        idx = db_index.get(str(chrom))
        if idx is None or idx.starts.size == 0:
            continue  # contig absent from DB -> every peak here is lost
        chrom_ok[rows] = True

        starts, ends, W = idx.starts, idx.ends, idx.max_width
        lo_all = np.searchsorted(starts, ps[rows] - W, side="left")
        hi_all = np.searchsorted(starts, pe[rows], side="left")

        for k, r in enumerate(rows):
            lo, hi = int(lo_all[k]), int(hi_all[k])
            if hi <= lo:
                continue
            cs = starts[lo:hi]
            ce = ends[lo:hi]
            ov = np.minimum(ce, pe[r]) - np.maximum(cs, ps[r])
            m = ov > 0
            if not m.any():
                continue
            ov = ov[m].astype(np.float64)
            cs, ce = cs[m], ce[m]

            ft = ov / plen[r]                      # overlap / len(user peak)
            fq = ov / (ce - cs).astype(np.float64)  # overlap / len(DB region)

            # pycistarget: strictly greater-than, reciprocal OR, per region.
            keep = (fq > fraction_overlap) | (ft > fraction_overlap)

            n_ov[r] = ov.size
            n_pass[r] = int(keep.sum())
            passes[r] = bool(keep.any())
            j = int(np.argmax(ov))
            best_ov[r] = int(ov[j])
            best_ft[r] = float(ft.max())
            best_fq[r] = float(fq.max())
            best_s[r] = int(cs[j])
            best_e[r] = int(ce[j])

    return {
        "best_ov_bp": best_ov,
        "best_frac_target": best_ft,
        "best_frac_query": best_fq,
        "passes": passes,
        "n_overlapping": n_ov,
        "n_passing": n_pass,
        "best_db_start": best_s,
        "best_db_end": best_e,
        "chrom_in_db": chrom_ok,
    }


# --------------------------------------------------------------------------- #
# Self-test: the audit engine vs brute force (and vs pyranges if available)
# --------------------------------------------------------------------------- #

def _brute_force(pchrom, pstart, pend, db_by_chrom, frac):
    """O(n*m) reference implementation of the same rule."""
    out_pass, out_ft, out_fq, out_ov = [], [], [], []
    for c, s, e in zip(pchrom, pstart, pend):
        bp, bft, bfq, bov = False, 0.0, 0.0, 0
        for ds, de in db_by_chrom.get(c, []):
            ov = min(de, e) - max(ds, s)
            if ov <= 0:
                continue
            ft = ov / (e - s)
            fq = ov / (de - ds)
            if fq > frac or ft > frac:
                bp = True
            bft, bfq = max(bft, ft), max(bfq, fq)
            bov = max(bov, ov)
        out_pass.append(bp)
        out_ft.append(bft)
        out_fq.append(bfq)
        out_ov.append(bov)
    return np.array(out_pass), np.array(out_ft), np.array(out_fq), np.array(out_ov)


def self_test(seed: int = 0) -> int:
    """Randomised differential test. Returns process exit code."""
    rng = np.random.default_rng(seed)
    n_fail = 0
    print("=== self-test: searchsorted engine vs brute force ===")

    for trial in range(200):
        frac = float(rng.choice([0.0, 0.2, 0.4, 0.5, 0.8, 1.0]))
        chroms = ["chrA", "chrB", "chrC"]

        # DB regions: mixed widths, including pathologically wide ones so the
        # max_width scan bound is actually exercised.
        db_ids, db_by_chrom = [], {}
        for c in chroms:
            m = int(rng.integers(0, 60))
            ss = rng.integers(0, 3000, size=m)
            ww = rng.integers(1, 400, size=m)
            if m and rng.random() < 0.3:
                ww[0] = int(rng.integers(400, 4000))  # one huge region
            ivs = sorted({(int(s), int(s + w)) for s, w in zip(ss, ww)})
            db_by_chrom[c] = ivs
            db_ids += [f"{c}:{s}-{e}" for s, e in ivs]

        # Peaks: include fixed-width 501 (the real use case) plus random widths,
        # plus a contig that does not exist in the DB.
        npk = int(rng.integers(1, 80))
        pchrom = [str(x) for x in rng.choice(chroms + ["chrZ"], size=npk)]
        pstart = [int(x) for x in rng.integers(0, 3200, size=npk)]
        pend = [
            s + (501 if rng.random() < 0.5 else int(rng.integers(1, 600)))
            for s in pstart
        ]

        idx, _ = build_db_index(db_ids, verbose=False)
        got = audit_peaks(pchrom, pstart, pend, idx, frac)
        exp_pass, exp_ft, exp_fq, exp_ov = _brute_force(
            pchrom, pstart, pend, db_by_chrom, frac)

        for label, a, b in (
            ("passes", got["passes"], exp_pass),
            ("best_ov_bp", got["best_ov_bp"], exp_ov),
        ):
            if not np.array_equal(a, b):
                n_fail += 1
                print(f"  FAIL trial {trial} frac={frac} field={label}")
        for label, a, b in (
            ("best_frac_target", got["best_frac_target"], exp_ft),
            ("best_frac_query", got["best_frac_query"], exp_fq),
        ):
            if not np.allclose(a, b, atol=1e-12):
                n_fail += 1
                print(f"  FAIL trial {trial} frac={frac} field={label}")

    print(f"  200 randomised trials, {n_fail} mismatch(es)")

    # Cross-check against a genuine pyranges join, if pyranges is installed.
    try:
        import pandas as pd
        import pyranges as pr
    except ImportError:
        print("=== pyranges cross-check: SKIPPED (pyranges not installed) ===")
    else:
        print("=== pyranges cross-check (pycistarget-style join) ===")
        pr_fail = 0
        for trial in range(25):
            frac = 0.4
            ivs = sorted({(int(s), int(s) + int(w)) for s, w in
                          zip(rng.integers(0, 5000, 150),
                              rng.integers(50, 400, 150))})
            db_ids = [f"chr1:{s}-{e}" for s, e in ivs]
            pstart = [int(x) for x in rng.integers(0, 5000, 120)]
            pend = [s + 501 for s in pstart]

            idx, _ = build_db_index(db_ids, verbose=False)
            got = audit_peaks(["chr1"] * len(pstart), pstart, pend, idx, frac)

            tgt = pr.PyRanges(pd.DataFrame({
                "Chromosome": ["chr1"] * len(pstart),
                "Start": pstart, "End": pend}))
            qry = pr.PyRanges(pd.DataFrame({
                "Chromosome": ["chr1"] * len(ivs),
                "Start": [s for s, _ in ivs], "End": [e for _, e in ivs]}))
            j = tgt.join(qry, report_overlap=True)
            if len(j) == 0:
                kept = set()
            else:
                df = j.df
                oq = df.Overlap / (df.End_b - df.Start_b)
                ot = df.Overlap / (df.End - df.Start)
                df = df[(oq > frac) | (ot > frac)]
                kept = set(zip(df.Start, df.End))
            mine = {(s, e) for s, e, p in zip(pstart, pend, got["passes"]) if p}
            if mine != kept:
                pr_fail += 1
                print(f"  FAIL trial {trial}: "
                      f"{len(mine - kept)} extra, {len(kept - mine)} missing")
        print(f"  25 trials, {pr_fail} mismatch(es)")
        n_fail += pr_fail

    print("SELF-TEST PASSED" if n_fail == 0 else f"SELF-TEST FAILED ({n_fail})")
    return 0 if n_fail == 0 else 1


# --------------------------------------------------------------------------- #
# Reporting
# --------------------------------------------------------------------------- #

def _hist(values: np.ndarray, edges: List[float]) -> Dict[str, int]:
    out = {}
    for i in range(len(edges) - 1):
        lo, hi = edges[i], edges[i + 1]
        last = i == len(edges) - 2
        m = (values >= lo) & (values <= hi) if last else (values >= lo) & (values < hi)
        out[f"[{lo:.1f},{hi:.1f}{']' if last else ')'}"] = int(m.sum())
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description="Audit what fraction of an ATAC peak set is representable "
                    "in a cisTarget database.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="See the module docstring for the exact overlap rule and its "
               "provenance in the pycistarget source.")
    ap.add_argument("--peaks", help="Consensus peak BED (may be .gz). "
                                    "BED3+ ; col4=name, col5=score if present.")
    ap.add_argument("--db", help="cisTarget *.feather (regions_vs_motifs "
                                 "rankings or scores).")
    ap.add_argument("--db-regions", help="Region catalog instead of --db: a "
                    "parquet/csv with chrom,start,end columns, or a text/BED "
                    "file of chr:start-end IDs. Use this to audit WITHOUT the "
                    "multi-GB feather on disk (see fetch_db_regions.py).")
    ap.add_argument("--out-prefix", default="peak_overlap_audit",
                    help="Output prefix for .per_peak.csv and .summary.json "
                         "(default: %(default)s)")
    ap.add_argument("--fraction-overlap", type=float, default=0.4,
                    help="Overlap threshold. MUST match SCENIC+'s "
                         "fraction_overlap_w_ctx_database / pycistarget's "
                         "fraction_overlap. (default: %(default)s)")
    ap.add_argument("--peak-type-column", choices=["name", "score", "none"],
                    default="name",
                    help="Which BED column to break the loss down by "
                         "(default: %(default)s)")
    ap.add_argument("--max-peak-types", type=int, default=50,
                    help="Skip the per-type breakdown above this many distinct "
                         "values, i.e. when the column is a unique peak ID "
                         "rather than a category (default: %(default)s)")
    ap.add_argument("--no-csv", action="store_true",
                    help="Skip the per-peak CSV (summary JSON only).")
    ap.add_argument("--self-test", action="store_true",
                    help="Validate the interval algorithm and exit.")
    args = ap.parse_args(argv)

    if args.self_test:
        return self_test()

    if not args.peaks:
        ap.error("missing required argument: --peaks (or pass --self-test)")
    if bool(args.db) == bool(args.db_regions):
        ap.error("give exactly one of --db (the feather) or --db-regions "
                 "(a pre-extracted region catalog)")
    if not 0.0 <= args.fraction_overlap < 1.0:
        ap.error("--fraction-overlap must be in [0, 1)")

    out_dir = os.path.dirname(os.path.abspath(args.out_prefix))
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    print("=" * 74)
    print("cisTarget peak representability audit")
    print("=" * 74)

    if args.db:
        region_ids = read_db_region_ids(args.db)
    else:
        region_ids = read_db_regions_catalog(args.db_regions)
    db_index, n_unparseable_db = build_db_index(region_ids)

    pchrom, pstart, pend, pnames, pscores, n_skipped = read_peaks_bed(args.peaks)
    n = len(pchrom)

    # -- contig naming sanity check (chr1 vs 1 is the classic silent killer) -- #
    peak_contigs = set(pchrom)
    db_contigs = set(db_index)
    shared = peak_contigs & db_contigs
    print(f"[check] peak contigs: {len(peak_contigs)}, "
          f"DB contigs: {len(db_contigs)}, shared: {len(shared)}")
    if not shared:
        print("\n" + "!" * 74)
        print("FATAL: no contig name is shared between your peaks and the "
              "database.")
        print(f"  peak contigs (first 5): {sorted(peak_contigs)[:5]}")
        print(f"  DB contigs  (first 5): {sorted(db_contigs)[:5]}")
        print("  This is almost always a 'chr1' vs '1' naming mismatch, or the")
        print("  wrong genome build. pycistarget would raise ValueError here.")
        print("!" * 74)
        return 2
    only_peaks = sorted(peak_contigs - db_contigs)
    if only_peaks:
        n_on_missing = sum(1 for c in pchrom if c in set(only_peaks))
        print(f"[check] {len(only_peaks)} peak contig(s) absent from the DB "
              f"({n_on_missing:,} peaks, {100 * n_on_missing / n:.2f}%): "
              f"{only_peaks[:8]}{' ...' if len(only_peaks) > 8 else ''}")

    print(f"\n[run] threshold: fraction_overlap > {args.fraction_overlap} "
          "(strict >, reciprocal OR)")
    res = audit_peaks(pchrom, pstart, pend, db_index, args.fraction_overlap)

    passes = res["passes"]
    n_pass = int(passes.sum())
    n_drop = n - n_pass
    no_ov = res["n_overlapping"] == 0
    n_no_ov = int(no_ov.sum())
    # Peaks that touch a DB region but not by enough: the subtle losses.
    n_marginal = int((~passes & ~no_ov).sum())

    best_frac = np.maximum(res["best_frac_target"], res["best_frac_query"])

    # -- per-chromosome breakdown ------------------------------------------- #
    pc = np.asarray(pchrom, dtype=object)
    per_chrom = {}
    for c in sorted(set(pchrom)):
        m = pc == c
        tot = int(m.sum())
        kept = int(passes[m].sum())
        per_chrom[str(c)] = {
            "n_peaks": tot,
            "n_representable": kept,
            "n_dropped": tot - kept,
            "frac_dropped": round((tot - kept) / tot, 6) if tot else 0.0,
            "n_zero_overlap": int(no_ov[m].sum()),
            "in_db": bool(str(c) in db_contigs),
        }

    # -- per-peak-type breakdown -------------------------------------------- #
    per_type = None
    type_source = None
    col = None
    if args.peak_type_column == "name" and pnames is not None:
        col, type_source = pnames, "BED name column (4)"
    elif args.peak_type_column == "score" and pscores is not None:
        col, type_source = pscores, "BED score column (5)"
    if col is not None:
        vals = np.asarray(col, dtype=object)
        distinct = set(col)
        if len(distinct) > args.max_peak_types:
            print(f"[type] {type_source} has {len(distinct):,} distinct values "
                  f"(> --max-peak-types {args.max_peak_types}); it looks like a "
                  "per-peak ID, not a category -- skipping type breakdown")
            type_source = (f"{type_source} -- skipped, "
                           f"{len(distinct)} distinct values")
        else:
            per_type = {}
            for v in sorted(distinct):
                m = vals == v
                tot = int(m.sum())
                kept = int(passes[m].sum())
                per_type[str(v)] = {
                    "n_peaks": tot,
                    "n_representable": kept,
                    "n_dropped": tot - kept,
                    "frac_dropped": round((tot - kept) / tot, 6) if tot else 0.0,
                }

    # -- how many DB regions actually get recruited ------------------------- #
    # pycistarget scores DB regions, not your peaks. Many peaks can collapse
    # onto one DB region, so the effective resolution can be lower than the
    # pass rate alone suggests.
    recruited = set()
    for s, e, ok in zip(res["best_db_start"], res["best_db_end"], passes):
        if ok and s >= 0:
            recruited.add((int(s), int(e)))

    # -- summary JSON -------------------------------------------------------- #
    summary = {
        "inputs": {
            "peaks_bed": os.path.abspath(args.peaks),
            "cistarget_db": os.path.abspath(args.db) if args.db else None,
            "db_regions_catalog": (os.path.abspath(args.db_regions)
                                   if args.db_regions else None),
            "fraction_overlap": args.fraction_overlap,
            "overlap_rule": ("(overlap/len_db_region > X) OR "
                             "(overlap/len_user_peak > X), strict >, "
                             "evaluated per candidate DB region "
                             "(pycistarget.utils.target_to_query)"),
        },
        "peak_set": {
            "n_peaks": n,
            "n_malformed_lines_skipped": n_skipped,
            "n_contigs": len(peak_contigs),
            "width_min": int((np.asarray(pend) - np.asarray(pstart)).min()),
            "width_median": float(np.median(np.asarray(pend) - np.asarray(pstart))),
            "width_max": int((np.asarray(pend) - np.asarray(pstart)).max()),
        },
        "database": {
            "n_regions": len(region_ids),
            "n_regions_indexed": int(sum(i.starts.size for i in db_index.values())),
            "n_region_ids_unparseable": n_unparseable_db,
            "n_contigs": len(db_contigs),
            "region_width_min": int(min(
                int((i.ends - i.starts).min()) for i in db_index.values())),
            "region_width_max": int(max(i.max_width for i in db_index.values())),
        },
        "verdict": {
            "n_representable": n_pass,
            "n_dropped": n_drop,
            "frac_representable": round(n_pass / n, 6),
            "frac_silently_dropped": round(n_drop / n, 6),
            "n_dropped_zero_overlap": n_no_ov,
            "n_dropped_insufficient_overlap": n_marginal,
            "n_unique_db_regions_recruited": len(recruited),
            "peak_to_db_region_collapse_ratio": (
                round(n_pass / len(recruited), 4) if recruited else None),
        },
        "best_overlap_fraction_distribution": _hist(
            best_frac, [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]),
        "best_overlap_fraction_quantiles": {
            q: round(float(np.quantile(best_frac, v)), 6)
            for q, v in (("p01", .01), ("p05", .05), ("p25", .25), ("p50", .50),
                         ("p75", .75), ("p95", .95), ("p99", .99))
        },
        "per_chromosome": per_chrom,
        "per_peak_type": per_type,
        "per_peak_type_source": type_source,
        "contigs_in_peaks_not_in_db": only_peaks,
    }

    json_path = f"{args.out_prefix}.summary.json"
    with open(json_path, "w") as fh:
        json.dump(summary, fh, indent=2)

    # -- per-peak CSV -------------------------------------------------------- #
    csv_path = f"{args.out_prefix}.per_peak.csv"
    if not args.no_csv:
        with open(csv_path, "w") as fh:
            hdr = ["peak_id", "chrom", "start", "end", "width"]
            if pnames is not None:
                hdr.append("name")
            if pscores is not None:
                hdr.append("score")
            hdr += ["best_overlap_bp", "best_frac_of_peak", "best_frac_of_db_region",
                    "best_frac_either", "n_db_regions_overlapping",
                    "n_db_regions_passing", "best_db_region", "chrom_in_db",
                    "representable"]
            fh.write(",".join(hdr) + "\n")
            for i in range(n):
                bs, be = int(res["best_db_start"][i]), int(res["best_db_end"][i])
                row = [f"{pchrom[i]}:{pstart[i]}-{pend[i]}", str(pchrom[i]),
                       str(pstart[i]), str(pend[i]), str(pend[i] - pstart[i])]
                if pnames is not None:
                    row.append('"' + pnames[i].replace('"', '""') + '"')
                if pscores is not None:
                    row.append('"' + pscores[i].replace('"', '""') + '"')
                row += [
                    str(int(res["best_ov_bp"][i])),
                    f"{res['best_frac_target'][i]:.6f}",
                    f"{res['best_frac_query'][i]:.6f}",
                    f"{best_frac[i]:.6f}",
                    str(int(res["n_overlapping"][i])),
                    str(int(res["n_passing"][i])),
                    f"{pchrom[i]}:{bs}-{be}" if bs >= 0 else "NA",
                    "1" if res["chrom_in_db"][i] else "0",
                    "1" if passes[i] else "0",
                ]
                fh.write(",".join(row) + "\n")

    # -- stdout verdict ------------------------------------------------------ #
    pct_drop = 100.0 * n_drop / n
    print("\n" + "=" * 74)
    print("VERDICT")
    print("=" * 74)
    print(f"  peaks audited                     : {n:,}")
    print(f"  representable in DB               : {n_pass:,} "
          f"({100.0 * n_pass / n:.2f}%)")
    print(f"  SILENTLY DROPPED                  : {n_drop:,} ({pct_drop:.2f}%)")
    print(f"    of which zero DB overlap        : {n_no_ov:,}")
    print(f"    of which overlap <= threshold   : {n_marginal:,}")
    print(f"  unique DB regions recruited       : {len(recruited):,}")
    if recruited:
        print(f"  peak -> DB region collapse        : "
              f"{n_pass / len(recruited):.2f} peaks per DB region")
    print(f"  median best-overlap fraction      : "
          f"{float(np.median(best_frac)):.3f}")

    worst = sorted(((v["frac_dropped"], k) for k, v in per_chrom.items()
                    if v["n_peaks"] >= 100), reverse=True)[:5]
    if worst:
        print("  worst chromosomes (>=100 peaks)   : " + ", ".join(
            f"{k} {100 * f:.1f}%" for f, k in worst))
    if per_type:
        wt = sorted(((v["frac_dropped"], k) for k, v in per_type.items()),
                    reverse=True)[:5]
        print("  worst peak types                  : " + ", ".join(
            f"{k} {100 * f:.1f}%" for f, k in wt))

    print("-" * 74)
    if pct_drop < 1:
        rec = ("THIS DATABASE COVERS YOUR PEAK SET. Essentially nothing is "
               "dropped.\n  (If you just built this DB from these same peaks, "
               "that is the expected\n  result and confirms the region IDs line "
               "up.)")
    elif pct_drop < 5:
        rec = ("THIS DATABASE IS ADEQUATE. Loss is negligible; building a "
               "custom DB\n  would cost hundreds of CPU-hours to recover <5% of "
               "your peaks.")
    elif pct_drop < 15:
        rec = ("THIS DATABASE IS PROBABLY FINE. Report the dropped "
               "fraction in\n  your methods. Build a custom DB only if the "
               "dropped peaks are\n  enriched in cell types or peak classes "
               "central to your question --\n  check the per-chromosome and "
               "per-type breakdowns above.")
    elif pct_drop < 30:
        rec = ("BORDERLINE -- INSPECT THE PER-TYPE BREAKDOWN. A loss this size "
               "will\n  measurably distort motif enrichment if it is not "
               "uniform. If the\n  dropped peaks are distal/non-cCRE (likely, "
               "since SCREEN is cCRE-based),\n  a custom database is justified.")
    else:
        rec = ("BUILD A CUSTOM DATABASE. You would lose more than 30% of your "
               "peak set\n  with no warning from SCENIC+. Run "
               "build_cistarget_db.sh on your own\n  consensus peaks.")
    print("RECOMMENDATION:\n  " + rec)
    print("-" * 74)
    print(f"  summary JSON : {json_path}")
    if not args.no_csv:
        print(f"  per-peak CSV : {csv_path}")
    print("=" * 74)
    return 0


if __name__ == "__main__":
    sys.exit(main())
