#!/usr/bin/env python3
"""
Attach transferred cell-type labels to ATAC cells, by VERIFIED positional merge.

THE SITUATION
-------------
`atac.h5ad` carries no cell-type annotation -- only `Clusters` (C1..C32) and QC.
The transferred labels live in `atac_metadata_with_transferred_labels.tsv`, which
has **no cell-barcode column at all**: 163,969 rows, 25 columns, none of them an
id. So the only possible join is positional, row i <-> cell i.

A positional merge on 163,969 cells is exactly the kind of assumption that fails
silently and poisons everything downstream. But it is checkable: 23 columns
appear in BOTH the TSV and `atac.h5ad.obs` (Clusters, TSSEnrichment, FRIP,
nFrags, Sample, ...). If TSV row i agrees with obs row i across all of them, the
row order is established as fact rather than hope.

This script does that check FIRST and refuses to emit anything if it fails.

WHAT IT WRITES
--------------
A small TSV keyed by the ATAC `obs_names`, carrying the label columns -- ready
for `--obs-tsv`. The 1.5 GB `atac.h5ad` is never rewritten.

Optionally (`--rna`), it derives a fine -> broad mapping from the RNA object,
which carries BOTH `predicted_CellType` (53) and `predicted_CellType_Broad`
(24), and applies it to the ATAC labels so both modalities can be grouped at
the broad level. Mapping purity is reported per label; anything below 1.0 is
listed so you can see the ambiguity rather than inherit it silently.

h5py + stdlib only.

Usage
-----
    python attach_atac_labels.py \
        --atac atac.h5ad \
        --tsv  atac_metadata_with_transferred_labels.tsv \
        --rna  rna.h5ad \
        --label-col predicted_CellType \
        --out  atac_labels.tsv
"""
from __future__ import annotations

import argparse
import csv
import sys
from collections import Counter, defaultdict
from pathlib import Path

try:
    import h5py
except ImportError:
    sys.exit("ERROR: needs h5py only.  pip install --user h5py")

FLOAT_TOL = 1e-4


def _dec(x):
    return x.decode("utf-8", "replace") if isinstance(x, bytes) else str(x)


def read_obs(path: Path):
    """Every obs column as a python list, plus obs_names. h5py only."""
    cols, names = {}, []
    with h5py.File(path, "r") as f:
        obs = f["obs"]
        idx = _dec(obs.attrs.get("_index", "_index"))
        names = [_dec(v) for v in obs[idx][:]]
        for k in obs.keys():
            if k == idx or k.startswith("__"):
                continue
            node = obs[k]
            if isinstance(node, h5py.Group) and "categories" in node:
                cats = [_dec(c) for c in node["categories"][:]]
                codes = node["codes"][:]
                cols[k] = [cats[c] if 0 <= c < len(cats) else "" for c in codes]
            else:
                try:
                    cols[k] = node[:].tolist()
                except Exception:                                # noqa: BLE001
                    pass
    return names, cols


def read_tsv(path: Path):
    with open(path, newline="") as fh:
        head = fh.readline()
        delim = "\t" if head.count("\t") >= head.count(",") else ","
        fh.seek(0)
        rdr = csv.reader(fh, delimiter=delim)
        header = next(rdr)
        rows = [r for r in rdr if r]
    widths = {len(r) for r in rows}
    if widths and min(widths) == len(header) + 1:
        header = ["__rowname__"] + header
        print("[attach] header had one fewer field than the data rows "
              "(R write.table(row.names=TRUE)); prepended __rowname__")
    cols = {c: [r[i] if i < len(r) else "" for r in rows]
            for i, c in enumerate(header)}
    return header, cols, len(rows)


def agree(a_list, b_list, name):
    """Column agreement, numeric-tolerant. Returns (n_mismatch, first_example)."""
    mism, example = 0, None
    for i, (a, b) in enumerate(zip(a_list, b_list)):
        sa, sb = str(a).strip(), str(b).strip()
        if sa == sb:
            continue
        try:
            fa, fb = float(sa), float(sb)
            if fa == fb or abs(fa - fb) <= FLOAT_TOL * max(1.0, abs(fa), abs(fb)):
                continue
            # ints written as "1" vs 1.0
            if fa == int(fb) or int(fa) == fb:
                continue
        except (ValueError, TypeError):
            pass
        mism += 1
        if example is None:
            example = (i, sa, sb)
    return mism, example


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--atac", type=Path, required=True)
    ap.add_argument("--tsv", type=Path, required=True)
    ap.add_argument("--rna", type=Path, default=None,
                    help="RNA .h5ad, used to derive a fine -> broad label map")
    ap.add_argument("--label-col", default="predicted_CellType")
    ap.add_argument("--rna-fine-col", default="predicted_CellType")
    ap.add_argument("--rna-broad-col", default="predicted_CellType_Broad")
    ap.add_argument("--extra-cols", nargs="*",
                    default=["predicted_CellType_confidence", "Cell_confidence"],
                    help="additional TSV columns to carry across")
    ap.add_argument("--min-witnesses", type=int, default=3,
                    help="refuse unless at least this many shared columns agree")
    ap.add_argument("--max-mismatch-frac", type=float, default=0.0,
                    help="tolerated per-column mismatch fraction (default 0 = exact)")
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()

    print(f"[attach] reading {args.atac}")
    obs_names, obs_cols = read_obs(args.atac)
    print(f"[attach] reading {args.tsv}")
    header, tsv_cols, n_rows = read_tsv(args.tsv)

    print(f"[attach] atac cells: {len(obs_names):,}   tsv rows: {n_rows:,}")
    if n_rows != len(obs_names):
        sys.exit(f"ERROR: row count {n_rows:,} != cell count {len(obs_names):,}. "
                 "A positional merge is impossible; these files do not correspond.")

    # ---- prove the row order ------------------------------------------------
    witnesses = sorted(set(obs_cols) & set(tsv_cols))
    if not witnesses:
        sys.exit("ERROR: no columns shared between atac.h5ad obs and the TSV, so "
                 "positional alignment cannot be verified. Refusing to guess.")
    print(f"[attach] verifying row order against {len(witnesses)} shared columns")
    ok, bad = [], []
    for c in witnesses:
        n_mis, ex = agree(obs_cols[c], tsv_cols[c], c)
        frac = n_mis / max(len(obs_names), 1)
        if frac <= args.max_mismatch_frac:
            ok.append((c, n_mis))
        else:
            bad.append((c, n_mis, frac, ex))
    for c, n in ok[:40]:
        print(f"    OK   {c}: {n} mismatches")
    for c, n, frac, ex in bad:
        print(f"    FAIL {c}: {n:,} mismatches ({frac:.2%}); "
              f"first at row {ex[0]}: obs={ex[1]!r} tsv={ex[2]!r}")
    if bad:
        sys.exit(f"\nERROR: {len(bad)} shared column(s) disagree, so TSV row i is "
                 "NOT cell i. Do not merge positionally. The labels must be joined "
                 "on a cell id, which this TSV does not carry -- regenerate it with "
                 "barcodes, or obtain the mapping from whoever produced it.")
    if len(ok) < args.min_witnesses:
        sys.exit(f"ERROR: only {len(ok)} column(s) could be checked; need "
                 f"{args.min_witnesses}. Too weak to trust a positional merge.")
    print(f"[attach] ROW ORDER VERIFIED: {len(ok)} columns agree on all "
          f"{len(obs_names):,} rows")

    if args.label_col not in tsv_cols:
        sys.exit(f"ERROR: '{args.label_col}' not in the TSV "
                 f"(have: {sorted(tsv_cols)})")
    labels = [v.strip() for v in tsv_cols[args.label_col]]
    cnt = Counter(labels)
    n_blank = cnt.get("", 0)
    print(f"[attach] '{args.label_col}': {len(cnt)} distinct values"
          + (f", {n_blank:,} BLANK (will be written as NA and dropped by "
             "--min-cells-per-group)" if n_blank else ""))

    # ---- optional fine -> broad map, derived from RNA -----------------------
    broad = None
    if args.rna:
        print(f"[attach] reading {args.rna} for the fine -> broad map")
        _, rna_cols = read_obs(args.rna)
        for c in (args.rna_fine_col, args.rna_broad_col):
            if c not in rna_cols:
                sys.exit(f"ERROR: '{c}' not in RNA obs (have: {sorted(rna_cols)})")
        pairs = defaultdict(Counter)
        for f_, b_ in zip(rna_cols[args.rna_fine_col], rna_cols[args.rna_broad_col]):
            pairs[f_][b_] += 1
        fine2broad, impure = {}, []
        for f_, c in pairs.items():
            b_, n = c.most_common(1)[0]
            tot = sum(c.values())
            fine2broad[f_] = b_
            if n < tot:
                impure.append((f_, b_, n / tot, dict(c.most_common(3))))
        print(f"[attach] derived {len(fine2broad)} fine -> broad mappings from RNA")
        if impure:
            print(f"[attach] {len(impure)} fine label(s) map to >1 broad label; "
                  "taking the modal broad label:")
            for f_, b_, pur, top in sorted(impure, key=lambda x: x[2])[:10]:
                print(f"    {f_!r} -> {b_!r} (purity {pur:.2f}) {top}")
        unmapped = sorted({l for l in labels if l and l not in fine2broad})
        if unmapped:
            print(f"[attach] WARNING: {len(unmapped)} ATAC label(s) have no RNA "
                  f"counterpart and will be blank at broad level: {unmapped[:10]}")
        broad = [fine2broad.get(l, "") if l else "" for l in labels]

    # ---- write ---------------------------------------------------------------
    out_cols = ["cell", args.label_col]
    data = {"cell": obs_names, args.label_col: labels}
    if broad is not None:
        out_cols.append(args.rna_broad_col)
        data[args.rna_broad_col] = broad
    for c in args.extra_cols:
        if c in tsv_cols and c != args.label_col:
            out_cols.append(c)
            data[c] = tsv_cols[c]

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(out_cols)
        for i in range(len(obs_names)):
            w.writerow([data[c][i] for c in out_cols])
    print(f"[attach] wrote {args.out}  ({len(obs_names):,} rows, "
          f"{len(out_cols)} columns)")

    if broad is not None:
        bc = Counter(b for b in broad if b)
        print(f"[attach] broad-level ATAC composition ({len(bc)} groups):")
        for b_, n in bc.most_common():
            print(f"    {n:>8,}  {b_}")
    print("\nNext: pass this to the pairing step with")
    print(f"    --obs-tsv {args.out} --group-key "
          f"{args.rna_broad_col if broad is not None else args.label_col}")


if __name__ == "__main__":
    main()
