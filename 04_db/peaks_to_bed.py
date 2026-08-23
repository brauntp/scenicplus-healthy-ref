#!/usr/bin/env python3
"""
Export peak coordinates from an .h5ad (or .h5mu) var_names to BED -- h5py only.

WHY THIS EXISTS
---------------
peak_overlap_audit.py needs a BED of your consensus peaks. Those coordinates are
already in the ATAC object's var_names (`chr:start-end`), so there is no reason
to load a 37 GB matrix -- or to have ArchR installed -- to get them. This reads
the var index dataset and nothing else: seconds on a login node, any file size.

Optionally attaches a per-peak label from an obs/var column so the audit can
break representability down by peak class.

Usage
-----
    # from an ATAC .h5ad
    python 04_db/peaks_to_bed.py --h5ad "$REF/atac.h5ad" --out consensus_peaks.bed

    # from the paired .h5mu (scATAC modality)
    python 04_db/peaks_to_bed.py --h5mu ACC_GEX.h5mu --out consensus_peaks.bed

Then:
    python 04_db/peak_overlap_audit.py --peaks consensus_peaks.bed \\
        --db-regions screen_db_regions.parquet --out-prefix audit_screen
"""
import argparse
import sys
from pathlib import Path

try:
    import h5py
except ImportError:
    sys.exit("ERROR: needs h5py and nothing else.\n"
             f"       interpreter: {sys.executable}\n"
             "       bash 00_inspect/find_inspect_env.sh finds one that has it")


def _index_key(grp):
    k = grp.attrs.get("_index", "_index")
    return k.decode() if isinstance(k, bytes) else str(k)


def _strings(dset):
    return [v.decode() if isinstance(v, bytes) else str(v) for v in dset[:]]


def parse_one(name):
    """chr:start-end, chr-start-end, or chr:start:end -> (chrom, start, end)."""
    s = str(name)
    if ":" in s:
        chrom, rest = s.split(":", 1)
        rest = rest.replace(":", "-")
    else:
        parts = s.rsplit("-", 2)
        if len(parts) != 3:
            return None
        chrom, a, b = parts
        rest = f"{a}-{b}"
    if "-" not in rest:
        return None
    a, b = rest.rsplit("-", 1)
    try:
        return chrom, int(a), int(b)
    except ValueError:
        return None


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--h5ad", type=Path)
    src.add_argument("--h5mu", type=Path)
    ap.add_argument("--modality", default="scATAC",
                    help="which .h5mu modality holds the peaks (default scATAC)")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--name-column", default=None,
                    help="var column to write as BED col4 (a peak class, if you "
                         "have one); default writes the region id")
    args = ap.parse_args()

    path = args.h5ad or args.h5mu
    if not path.exists():
        sys.exit(f"ERROR: {path} does not exist")

    with h5py.File(path, "r") as f:
        root = f if args.h5ad else f["mod"][args.modality]
        if args.h5mu and args.modality not in f.get("mod", {}):
            sys.exit(f"ERROR: modality '{args.modality}' not in "
                     f"{list(f.get('mod', {}).keys())}")
        var = root["var"]
        names = _strings(var[_index_key(var)])
        labels = None
        if args.name_column:
            if args.name_column not in var:
                sys.exit(f"ERROR: --name-column '{args.name_column}' not in var "
                         f"(have: {list(var.keys())})")
            g = var[args.name_column]
            if isinstance(g, h5py.Group) and "categories" in g:
                cats = _strings(g["categories"])
                labels = [cats[c] if c >= 0 else "NA" for c in g["codes"][:]]
            else:
                labels = _strings(g)

    parsed, bad = [], []
    for i, n in enumerate(names):
        p = parse_one(n)
        if p is None:
            bad.append(n)
            continue
        chrom, a, b = p
        if b <= a:
            bad.append(n)
            continue
        parsed.append((chrom, a, b, labels[i] if labels else n))

    if not parsed:
        sys.exit(f"ERROR: no var_names parsed as coordinates. First few: "
                 f"{names[:3]} -- is this a peak matrix?")

    # Sort by coordinate so the BED is well formed for any downstream tool.
    parsed.sort(key=lambda r: (r[0], r[1], r[2]))
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w") as fh:
        for chrom, a, b, nm in parsed:
            fh.write(f"{chrom}\t{a}\t{b}\t{nm}\n")

    widths = [b - a for _, a, b, _ in parsed]
    chroms = sorted({c for c, _, _, _ in parsed})
    print(f"wrote {args.out}: {len(parsed):,} peaks")
    print(f"  widths: min {min(widths)} median {sorted(widths)[len(widths)//2]} "
          f"max {max(widths)}")
    print(f"  {len(chroms)} contigs; first few: {chroms[:6]}")
    if bad:
        print(f"  WARNING: {len(bad):,} var_names did not parse as coordinates "
              f"and were skipped, e.g. {bad[:3]}", file=sys.stderr)
    print()
    print("Next:")
    print(f"  python 04_db/peak_overlap_audit.py --peaks {args.out} \\")
    print(f"      --db-regions screen_db_regions.parquet --out-prefix audit_screen")


if __name__ == "__main__":
    main()
