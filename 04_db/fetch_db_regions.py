#!/usr/bin/env python3
"""
Extract the region catalog from a remote cisTarget feather WITHOUT downloading it.

WHY
---
`peak_overlap_audit.py` needs only the region IDs of a cisTarget database to
decide whether your peak set is representable in it. The precomputed hg38 SCREEN
rankings feather is ~33 GiB; its region IDs live in the Arrow schema, which sits
in the file's footer. Feather v2 is Arrow IPC: the last 10 bytes are
`<int32 footer_length><6-byte "ARROW1">`, so one HTTP range request for the tail
gives the footer length, and a second gives the footer itself (~121 MB for the
SCREEN database -- roughly 276x smaller than the full file).

This means the "custom database or precomputed?" question can be answered before
committing ~46 GiB of download and the disk to hold it.

Requires the server to honour byte ranges (resources.aertslab.org does:
`Accept-Ranges: bytes`). Falls back with a clear message if it does not.

Usage
-----
    # default: hg38 SCREEN region-based rankings
    python fetch_db_regions.py --out screen_db_regions.parquet

    # any other cisTarget feather URL
    python fetch_db_regions.py --url <...feather> --out regions.parquet

    # audit a local feather you already have, no network
    python fetch_db_regions.py --local /path/to/db.feather --out regions.parquet

Then:
    python peak_overlap_audit.py --peaks consensus_peaks.bed \
        --db-regions screen_db_regions.parquet
"""
from __future__ import annotations

import argparse
import os
import re
import struct
import sys
import urllib.error
import urllib.request

DEFAULT_URL = (
    "https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/"
    "screen/mc_v10_clust/region_based/"
    "hg38_screen_v10_clust.regions_vs_motifs.rankings.feather"
)
# Sanity anchor: the published hg38 SCREEN database carries this many regions.
KNOWN_COUNTS = {"hg38_screen_v10_clust": 1_837_304}
REGION_RE = re.compile(rb"(chr[0-9XYM]{1,2}|[0-9]{1,2}|[XYM]):([0-9]{1,10})-([0-9]{1,10})")


def _log(m: str) -> None:
    print(f"[fetch_db_regions] {m}", flush=True)


def _range_get(url: str, start: int, end: int, timeout: int = 1800) -> bytes:
    """Inclusive byte range."""
    req = urllib.request.Request(url, headers={"Range": f"bytes={start}-{end}"})
    try:
        with urllib.request.urlopen(req, timeout=timeout) as r:
            if r.status != 206:
                sys.exit(f"ERROR: server ignored the range request (HTTP {r.status}). "
                         "It does not support partial content; download the feather "
                         "with download_precomputed_db.sh and use --db instead.")
            buf = bytearray()
            while chunk := r.read(8 << 20):
                buf += chunk
            return bytes(buf)
    except urllib.error.HTTPError as e:
        sys.exit(f"ERROR: HTTP {e.code} fetching range from {url}\n{e.read()[:300]!r}")
    except urllib.error.URLError as e:
        sys.exit(f"ERROR: could not reach {url}: {e.reason}\n"
                 "If this host is not on your network allowlist, download the "
                 "feather on a machine that can reach it and use --local.")


def footer_from_url(url: str) -> bytes:
    with urllib.request.urlopen(urllib.request.Request(url, method="HEAD"),
                                timeout=120) as r:
        size = int(r.headers["Content-Length"])
        if r.headers.get("Accept-Ranges", "").lower() != "bytes":
            _log("WARNING: server did not advertise Accept-Ranges; trying anyway")
    _log(f"remote size {size / 1024**3:.1f} GiB")
    tail = _range_get(url, size - 64, size - 1)
    if tail[-6:] != b"ARROW1":
        sys.exit("ERROR: file does not end with the ARROW1 magic; not a Feather v2 "
                 "file. Use --db with the real feather instead.")
    flen = struct.unpack("<i", tail[-10:-6])[0]
    _log(f"footer metadata is {flen / 1024**2:.1f} MB "
         f"({size / max(flen, 1):.0f}x smaller than the file)")
    return _range_get(url, size - 10 - flen, size - 11)


def footer_from_local(path: str) -> bytes:
    size = os.path.getsize(path)
    with open(path, "rb") as fh:
        fh.seek(size - 64)
        tail = fh.read(64)
        if tail[-6:] != b"ARROW1":
            sys.exit("ERROR: not a Feather v2 file (no ARROW1 magic).")
        flen = struct.unpack("<i", tail[-10:-6])[0]
        fh.seek(size - 10 - flen)
        return fh.read(flen)


def parse_regions(raw: bytes):
    """Region IDs appear as schema field names in the footer flatbuffer."""
    hits = REGION_RE.findall(raw)
    seen, out = set(), []
    for c, s, e in hits:
        key = (c, s, e)
        if key in seen:
            continue
        seen.add(key)
        out.append((c.decode(), int(s), int(e)))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    src = ap.add_mutually_exclusive_group()
    src.add_argument("--url", default=DEFAULT_URL, help="cisTarget feather URL")
    src.add_argument("--local", help="path to a local feather (no network)")
    ap.add_argument("--out", default="db_regions.parquet",
                    help="output .parquet (or .csv) region catalog")
    ap.add_argument("--expect", type=int, default=None,
                    help="assert this many regions were parsed")
    args = ap.parse_args()

    raw = footer_from_local(args.local) if args.local else footer_from_url(args.url)
    regions = parse_regions(raw)
    if not regions:
        sys.exit("ERROR: no chr:start-end region IDs found in the footer. This may "
                 "be a GENE-based database, or a layout where regions are stored "
                 "as column VALUES rather than field names -- use --db instead.")
    _log(f"parsed {len(regions):,} unique region IDs")

    name = os.path.basename(args.local or args.url)
    expect = args.expect
    if expect is None:
        for k, v in KNOWN_COUNTS.items():
            if k in name:
                expect = v
                break
    if expect is not None:
        if len(regions) != expect:
            sys.exit(f"ERROR: parsed {len(regions):,} regions but expected "
                     f"{expect:,} for {name}. Refusing to write a catalog that "
                     "does not match the known database size.")
        _log(f"count matches the known size for {name} ({expect:,})")

    try:
        import pandas as pd
    except ImportError:
        sys.exit("ERROR: pandas is required to write the catalog.")
    df = pd.DataFrame(regions, columns=["chrom", "start", "end"])
    df["width"] = df.end - df.start
    if args.out.endswith(".parquet"):
        df.to_parquet(args.out, index=False)
    else:
        df.to_csv(args.out, index=False)
    _log(f"wrote {args.out} ({os.path.getsize(args.out) / 1024**2:.1f} MB)")
    _log(f"region width: median {df.width.median():.0f} bp, "
         f"range {df.width.min()}-{df.width.max()}; "
         f"{df.chrom.nunique()} contigs")
    print()
    print("Next:")
    print(f"  python peak_overlap_audit.py --peaks <consensus_peaks.bed> "
          f"--db-regions {args.out}")


if __name__ == "__main__":
    main()
