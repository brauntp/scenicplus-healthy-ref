#!/usr/bin/env python3
"""
Dependency-light .h5ad inspector -- needs ONLY h5py.

WHY THIS EXISTS
---------------
inspect_anndata.py needs anndata (hence numpy/scipy/pandas). On a cluster the
env that wrote the files may belong to a different user, or may not be findable
at all. But .h5ad is just HDF5 with a documented layout, so everything stage 1
actually needs can be read straight off the file.

This answers exactly six questions, and nothing else:

  1. shapes                -- how many cells and features, per modality
  2. ATAC var_names format -- peaks (chr:start-end) or gene/tile names?
  3. obsm keys + dims      -- is the GLUE latent IN the file, or only in a TSV?
  4. obs columns + levels  -- which column can be key_to_group_by?
  5. cell-id format        -- do ArchR ids and GLUE ids match?
  6. X dtype / integrality -- raw counts or normalized?

Usage
-----
    python inspect_h5ad_lite.py rna.h5ad atac.h5ad --out report_lite
    python inspect_h5ad_lite.py "$REF"/*.h5ad --out report_lite

If h5py is missing:  pip install h5py   (binary wheel, no compiler, seconds)

Emits <out>.json and <out>.md, and prints a summary. Reads only structure,
names and a small sample of values -- never a whole matrix.
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

try:
    import h5py
except ImportError:
    sys.exit("ERROR: this needs h5py and nothing else.\n"
             f"       interpreter: {sys.executable}\n"
             "       pip install h5py    (binary wheel, no compiler, seconds)")

PEAK_RE = re.compile(r"^(chr)?[0-9XYMT]+[:\-_]\d+[\-_]\d+$", re.I)
MAX_LEVELS = 60


# ---------------------------------------------------------------- helpers
def _decode(x):
    if isinstance(x, bytes):
        return x.decode("utf-8", "replace")
    return str(x)


def _read_strings(node, limit=None):
    """Read a string dataset (or a categorical group) into a python list."""
    if isinstance(node, h5py.Group):
        # categorical: {categories, codes}
        if "categories" in node and "codes" in node:
            cats = [_decode(c) for c in node["categories"][:]]
            codes = node["codes"][:limit] if limit else node["codes"][:]
            return [cats[c] if 0 <= c < len(cats) else None for c in codes], cats
        return None, None
    data = node[:limit] if limit else node[:]
    return [_decode(v) for v in data], None


def _x_info(grp):
    """Shape + dtype + integrality for X, dense or sparse, without loading it."""
    info = {}
    if isinstance(grp, h5py.Group):                       # sparse
        enc = grp.attrs.get("encoding-type", b"")
        info["storage"] = _decode(enc) or "sparse"
        shp = grp.attrs.get("shape")
        if shp is not None:
            info["shape"] = [int(v) for v in shp]
        if "data" in grp:
            d = grp["data"]
            info["dtype"] = str(d.dtype)
            info["nnz"] = int(d.shape[0])
            s = d[: min(200_000, d.shape[0])]
            info["all_integer_sample"] = bool((s == s.astype("int64")).all()) \
                if s.size else None
            if s.size:
                info["min"], info["max"] = float(s.min()), float(s.max())
    else:                                                  # dense
        info["storage"] = "dense"
        info["shape"] = [int(v) for v in grp.shape]
        info["dtype"] = str(grp.dtype)
        n = min(2000, grp.shape[0])
        s = grp[:n]
        info["all_integer_sample"] = bool((s == s.astype("int64")).all())
        info["min"], info["max"] = float(s.min()), float(s.max())
    if info.get("all_integer_sample") is True:
        info["looks_like"] = "raw_counts"
    elif info.get("all_integer_sample") is False:
        info["looks_like"] = ("log1p_or_normalized"
                              if info.get("max", 1e9) < 30 else "normalized_nonlog")
    return info


def inspect(path: Path) -> dict:
    rep = {"file": str(path), "size_gb": round(path.stat().st_size / 1024**3, 2)}
    with h5py.File(path, "r") as f:
        rep["root_keys"] = sorted(f.keys())

        # ---- X -----------------------------------------------------------
        if "X" in f:
            rep["X"] = _x_info(f["X"])
        rep["layers"] = sorted(f["layers"].keys()) if "layers" in f else []
        rep["has_raw"] = "raw" in f
        if rep["has_raw"] and "X" in f["raw"]:
            rep["raw_X"] = _x_info(f["raw"]["X"])

        # ---- obs ----------------------------------------------------------
        obs = f.get("obs")
        if obs is not None:
            idx_key = _decode(obs.attrs.get("_index", "_index"))
            names, _ = _read_strings(obs[idx_key], limit=5) if idx_key in obs else (None, None)
            rep["n_obs"] = int(obs[idx_key].shape[0]) if idx_key in obs else None
            rep["obs_index_name"] = idx_key
            rep["obs_names_sample"] = names
            # A larger sample, kept out of the printed report, so TSV comparison
            # can test real membership instead of guessing from id shape.
            if idx_key in obs:
                big, _ = _read_strings(obs[idx_key], limit=3000)
                rep["_obs_id_sample"] = set(big or [])
            if names:
                rep["cell_id_has_hash"] = any("#" in n for n in names)
                rep["cell_id_has_suffix_dash"] = any(re.search(r"-\d+$", n) for n in names)
            cols = {}
            for k in obs.keys():
                if k == idx_key or k.startswith("__"):
                    continue
                node = obs[k]
                entry = {}
                if isinstance(node, h5py.Group) and "categories" in node:
                    cats = [_decode(c) for c in node["categories"][:]]
                    entry["kind"] = "categorical"
                    entry["n_levels"] = len(cats)
                    entry["levels"] = cats[:MAX_LEVELS]
                    if len(cats) > MAX_LEVELS:
                        entry["levels_truncated"] = True
                else:
                    entry["kind"] = str(getattr(node, "dtype", "group"))
                cols[k] = entry
            rep["obs_columns"] = cols

        # ---- var ----------------------------------------------------------
        var = f.get("var")
        if var is not None:
            vidx = _decode(var.attrs.get("_index", "_index"))
            rep["n_var"] = int(var[vidx].shape[0]) if vidx in var else None
            vnames, _ = _read_strings(var[vidx], limit=2000) if vidx in var else (None, None)
            rep["var_names_sample"] = vnames[:5] if vnames else None
            rep["var_columns"] = sorted(k for k in var.keys() if k != vidx)
            if vnames:
                hits = [bool(PEAK_RE.match(v)) for v in vnames]
                frac = sum(hits) / len(hits)
                rep["var_names_peaklike_frac"] = round(frac, 3)
                rep["feature_type_guess"] = (
                    "PEAKS (chr:start-end)" if frac > 0.9
                    else "GENES / tiles / other (NOT a peak matrix)")
                if frac > 0.9:
                    widths = []
                    for v in vnames:
                        m = re.match(r"^(chr)?[0-9XYMT]+[:\-_](\d+)[\-_](\d+)$", v, re.I)
                        if m:
                            widths.append(int(m.group(3)) - int(m.group(2)))
                    if widths:
                        widths.sort()
                        rep["peak_width"] = {
                            "min": widths[0], "median": widths[len(widths)//2],
                            "max": widths[-1],
                            "fixed_width": widths[0] == widths[-1],
                            "n_sampled": len(widths)}
                    rep["chr_prefix"] = vnames[0].lower().startswith("chr")

        # ---- obsm ----------------------------------------------------------
        if "obsm" in f:
            shapes = {}
            for k in f["obsm"].keys():
                node = f["obsm"][k]
                shapes[k] = [int(v) for v in node.shape] if hasattr(node, "shape") else None
            rep["obsm"] = shapes
            rep["glue_obsm_keys"] = [k for k in shapes if "glue" in k.lower()]
        else:
            rep["obsm"] = {}
            rep["glue_obsm_keys"] = []
    return rep


def to_md(reps):
    L = ["# .h5ad structural report (h5py only)", ""]
    for r in reps:
        L += [f"## `{Path(r['file']).name}`  ({r['size_gb']} GB)", ""]
        L.append(f"- **{r.get('n_obs')} cells x {r.get('n_var')} features**")
        x = r.get("X", {})
        if x:
            L.append(f"- X: {x.get('storage')} {x.get('dtype')} "
                     f"-> looks like **{x.get('looks_like','?')}** "
                     f"(min {x.get('min')}, max {x.get('max')})")
        if r.get("has_raw"):
            rx = r.get("raw_X", {})
            L.append(f"- `.raw` present: {rx.get('storage')} {rx.get('dtype')} "
                     f"-> {rx.get('looks_like','?')}")
        if r.get("layers"):
            L.append(f"- layers: {', '.join(r['layers'])}")
        L.append(f"- features look like: **{r.get('feature_type_guess','?')}** "
                 f"(peak-like fraction {r.get('var_names_peaklike_frac')})")
        if r.get("peak_width"):
            pw = r["peak_width"]
            L.append(f"  - peak width: median {pw['median']} bp, "
                     f"range {pw['min']}-{pw['max']}, fixed={pw['fixed_width']}")
            L.append(f"  - 'chr' prefix: {r.get('chr_prefix')}")
        L.append(f"- var_names e.g. `{r.get('var_names_sample')}`")
        L.append(f"- obs_names e.g. `{r.get('obs_names_sample')}` "
                 f"(contains '#': {r.get('cell_id_has_hash')})")
        L.append(f"- **obsm**: {r.get('obsm') or 'NONE'}  "
                 f"-> GLUE-like keys: {r.get('glue_obsm_keys') or 'NONE'}")
        L += ["", "### obs columns", ""]
        for k, v in (r.get("obs_columns") or {}).items():
            if v.get("kind") == "categorical":
                lv = v.get("levels", [])
                L.append(f"- `{k}` -- categorical, {v['n_levels']} levels: "
                         f"{', '.join(map(str, lv[:12]))}"
                         f"{' ...' if v.get('levels_truncated') or len(lv) > 12 else ''}")
            else:
                L.append(f"- `{k}` -- {v.get('kind')}")
        L.append("")
    return "\n".join(L)


def peek_tsv(path: Path, n_ids: int = 5) -> dict:
    """
    Structure of a sidecar TSV, read with the csv module only.

    The question this answers: are the row ids the SAME ids as the .h5ad
    obs_names? If the GLUE embedding TSV uses plain barcodes while the .h5ad
    uses Sample#BARCODE-1, pairing will refuse and the fix is to harmonise the
    ids -- better to learn that now than at stage 5.
    """
    import csv
    rep = {"file": str(path), "size_mb": round(path.stat().st_size / 1024**2, 1)}
    with open(path, newline="") as fh:
        sample = fh.read(64_000)
        fh.seek(0)
        delim = "\t" if sample.count("\t") >= sample.count(",") else ","
        rep["delimiter"] = "tab" if delim == "\t" else "comma"
        rdr = csv.reader(fh, delimiter=delim)
        header = next(rdr, [])
        rep["n_columns"] = len(header)
        rep["header_first"] = header[:8]
        rep["header_last"] = header[-4:] if len(header) > 8 else []
        ids, nrow, widths = [], 0, set()
        id_set = set()
        for row in rdr:
            if row:
                if nrow < n_ids:
                    ids.append(row[0])
                if nrow < 200_000:
                    id_set.add(row[0])
                widths.add(len(row))
            nrow += 1
        rep["n_rows"] = nrow
        rep["first_ids"] = ids
        rep["_id_sample"] = id_set          # for real overlap testing, not printed
        rep["row_widths"] = sorted(widths)[:4]
        # R's write.table(row.names=TRUE) emits a header with ONE FEWER field
        # than each data row. Detected here because it silently shifts every
        # column by one when read with a naive header=0 parser.
        if widths and len(header) == (min(widths) - 1):
            rep["r_rownames_quirk"] = True
            rep["note"] = ("header has one fewer field than the data rows -- this is "
                           "R write.table(row.names=TRUE). The first data field is the "
                           "ROW NAME (likely the cell id); read with index_col=0.")
        else:
            rep["r_rownames_quirk"] = False
        rep["id_has_hash"] = any("#" in i for i in ids)
        rep["id_has_suffix_dash"] = any(re.search(r"-\d+$", i) for i in ids)
        # numeric width = plausible embedding dimensionality
        numeric = 0
        for c in header[1:]:
            try:
                float(c)
                numeric += 1
            except ValueError:
                pass
        rep["header_all_numeric"] = (numeric == len(header) - 1 and numeric > 0)
    return rep


def tsv_md(tsvs, h5reps):
    if not tsvs:
        return ""
    L = ["", "## sidecar TSVs", ""]
    h5_ids = {Path(r["file"]).name: (r.get("obs_names_sample") or [None])[0]
              for r in h5reps}
    for t in tsvs:
        L.append(f"### `{Path(t['file']).name}`  ({t['size_mb']} MB)")
        L.append(f"- {t['n_rows']} rows x {t['n_columns']} columns "
                 f"({t['delimiter']}-separated)")
        L.append(f"- header starts: `{t['header_first']}`"
                 + (f" ... ends `{t['header_last']}`" if t['header_last'] else ""))
        L.append(f"- first ids: `{t['first_ids']}` "
                 f"(contains '#': {t['id_has_hash']})")
        L.append("")
        # Compare ACTUAL ids, not format flags. Two id sets can share a format
        # and still be disjoint (integer row numbers vs opaque hashes), which a
        # flags-only check reports as compatible -- exactly backwards.
        for r in h5reps:
            fn = Path(r["file"]).name
            obs_sample = set(r.get("_obs_id_sample") or [])
            ex = (r.get("obs_names_sample") or [None])[0]
            if not obs_sample:
                continue
            hit = len(obs_sample & t.get("_id_sample", set()))
            frac = hit / max(len(obs_sample), 1)
            if frac >= 0.99:
                verdict = f"**IDS MATCH** ({hit}/{len(obs_sample)} sampled ids found)"
            elif hit:
                verdict = (f"**PARTIAL OVERLAP** -- only {hit}/{len(obs_sample)} "
                           f"sampled ids found ({frac:.0%})")
            else:
                verdict = ("**NO SHARED IDS** -- 0 of "
                           f"{len(obs_sample)} sampled ids appear in this TSV. "
                           "These are different id spaces; they must be mapped, "
                           "not just reformatted.")
            L.append(f"  - vs `{fn}` obs_names (e.g. `{ex}`): {verdict}")
        L.append("")
    return "\n".join(L)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("files", nargs="+", type=Path, help=".h5ad file(s)")
    ap.add_argument("--tsv", type=Path, nargs="*", default=[],
                    help="sidecar TSV(s) to peek at (embeddings / metadata): "
                         "reports the header, row count and first ids, so you can "
                         "see whether the ids match the .h5ad obs_names")
    ap.add_argument("--out", type=Path, default=Path("report_lite"))
    args = ap.parse_args()

    reps = []
    for p in args.files:
        if not p.exists():
            print(f"WARNING: missing {p}", file=sys.stderr)
            continue
        print(f"reading {p} ...", flush=True)
        reps.append(inspect(p))
    if not reps:
        sys.exit("ERROR: no readable input files")

    tsvs = [peek_tsv(p) for p in args.tsv if p.exists()]
    for p in args.tsv:
        if not p.exists():
            print(f"WARNING: missing {p}", file=sys.stderr)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    def _strip(d):
        return {k: v for k, v in d.items() if not k.startswith("_")}
    with open(f"{args.out}.json", "w") as fh:
        json.dump({"h5ad": [_strip(r) for r in reps],
                   "tsv": [_strip(t) for t in tsvs]}, fh, indent=2, default=str)
    md = to_md(reps) + tsv_md(tsvs, reps)
    with open(f"{args.out}.md", "w") as fh:
        fh.write(md)
    print()
    print(md)
    print(f"\nwrote {args.out}.json and {args.out}.md")


if __name__ == "__main__":
    main()
