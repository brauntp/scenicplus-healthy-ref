#!/usr/bin/env python3
"""
Dump and compare cell-type label vocabularies across .h5ad files and TSVs.

WHY
---
Pairing is confined WITHIN a shared grouping key, so RNA and ATAC must carry
the same label column with the SAME level names. Three things can go wrong, and
all three are silent until something refuses hours later:

  1. the label is not in the .h5ad at all (it lives in a sidecar TSV)
  2. the level names differ ("HSC MPP" vs "HSC/MPP" vs "HSC_MPP")
  3. the TSV's row ids are not the object's obs_names, so a merge yields NaN

This reads every categorical obs column, every plausible label column in the
TSVs, and reports the vocabularies side by side with an exact overlap count.

h5py + stdlib only.

Usage
-----
    python compare_labels.py --h5ad rna.h5ad atac.h5ad \
                             --tsv atac_metadata_with_transferred_labels.tsv \
                             --out label_report
"""
from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from pathlib import Path

try:
    import h5py
except ImportError:
    sys.exit("ERROR: needs h5py only.  pip install --user h5py")

# Columns whose names suggest a cell-type / lineage annotation.
LABEL_HINT = re.compile(
    r"(cell.?type|celltype|lineage|annotation|label|cluster|ident|phenotype|"
    r"predicted|initial|broad|class)", re.I)
# Columns that look like a per-cell identifier.
ID_HINT = re.compile(r"(barcode|cell|^index$|^id$|rowname)", re.I)


def _dec(x):
    return x.decode("utf-8", "replace") if isinstance(x, bytes) else str(x)


def read_h5ad_labels(path: Path, id_sample=5000) -> dict:
    out = {"file": str(path), "categorical": {}, "obs_ids": [], "n_obs": None}
    with h5py.File(path, "r") as f:
        obs = f.get("obs")
        if obs is None:
            return out
        idx = _dec(obs.attrs.get("_index", "_index"))
        if idx in obs:
            n = obs[idx].shape[0]
            out["n_obs"] = int(n)
            out["obs_ids"] = [_dec(v) for v in obs[idx][: min(id_sample, n)]]
        for k in obs.keys():
            if k == idx or k.startswith("__"):
                continue
            node = obs[k]
            if isinstance(node, h5py.Group) and "categories" in node:
                cats = [_dec(c) for c in node["categories"][:]]
                # counts per level, from the codes
                counts = {}
                if "codes" in node:
                    codes = node["codes"][:]
                    import collections
                    cc = collections.Counter(codes.tolist())
                    counts = {cats[i]: int(cc.get(i, 0))
                              for i in range(len(cats)) if cc.get(i, 0)}
                out["categorical"][k] = {"n_levels": len(cats),
                                         "levels": cats,
                                         "counts": counts,
                                         "label_like": bool(LABEL_HINT.search(k))}
    return out


def read_tsv_labels(path: Path, max_rows=None) -> dict:
    """Parse a TSV, handling R's write.table(row.names=TRUE) off-by-one header."""
    out = {"file": str(path)}
    with open(path, newline="") as fh:
        head = fh.readline()
        delim = "\t" if head.count("\t") >= head.count(",") else ","
        fh.seek(0)
        rdr = csv.reader(fh, delimiter=delim)
        header = next(rdr, [])
        rows = []
        widths = set()
        for i, row in enumerate(rdr):
            if row:
                widths.add(len(row))
                if max_rows is None or i < max_rows:
                    rows.append(row)
        out["n_rows"] = i + 1 if rows or widths else 0

    w = min(widths) if widths else len(header)
    rownames_shift = (w == len(header) + 1)
    out["r_rownames_quirk"] = rownames_shift
    cols = (["__rowname__"] + header) if rownames_shift else header
    out["columns"] = cols
    out["n_columns"] = len(cols)

    # Per-column value sets (only for plausible label / id columns; cheap).
    vocab, ids = {}, []
    import collections
    for ci, cname in enumerate(cols):
        if not (LABEL_HINT.search(cname) or ID_HINT.search(cname)
                or cname == "__rowname__"):
            continue
        vals = [r[ci] for r in rows if ci < len(r)]
        cnt = collections.Counter(vals)
        uniq = sorted(cnt)
        entry = {"n_unique": len(uniq),
                 "label_like": bool(LABEL_HINT.search(cname))}
        if len(uniq) <= 200:
            entry["levels"] = uniq
            entry["counts"] = {k: int(v) for k, v in cnt.most_common()}
        else:
            entry["levels_sample"] = uniq[:10]
        vocab[cname] = entry
        if cname == "__rowname__" or ID_HINT.search(cname):
            ids = vals[:5000] if len(uniq) > len(vals) * 0.5 else ids
    out["columns_of_interest"] = vocab
    out["candidate_ids"] = ids[:5000]
    return out


def jaccard(a, b):
    a, b = set(a), set(b)
    return len(a & b) / max(len(a | b), 1)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--h5ad", type=Path, nargs="+", required=True)
    ap.add_argument("--tsv", type=Path, nargs="*", default=[])
    ap.add_argument("--out", type=Path, default=Path("label_report"))
    args = ap.parse_args()

    h5 = []
    for p in args.h5ad:
        if not p.exists():
            print(f"WARNING: missing {p}", file=sys.stderr)
            continue
        print(f"reading {p} ...", flush=True)
        h5.append(read_h5ad_labels(p))
    tsvs = []
    for p in args.tsv:
        if not p.exists():
            print(f"WARNING: missing {p}", file=sys.stderr)
            continue
        print(f"reading {p} ...", flush=True)
        tsvs.append(read_tsv_labels(p))

    L = ["# Label vocabulary report", ""]

    # ---- per-file label-like columns ------------------------------------
    for r in h5:
        name = Path(r["file"]).name
        L += [f"## `{name}`  ({r['n_obs']:,} cells)" if r["n_obs"] else f"## `{name}`", ""]
        lab = {k: v for k, v in r["categorical"].items() if v["label_like"]}
        if not lab:
            L.append("**No label-like categorical column in `obs`.** "
                     "The annotation must come from a sidecar TSV.")
            L.append("")
        for k, v in sorted(lab.items(), key=lambda kv: kv[1]["n_levels"]):
            L.append(f"- **`{k}`** — {v['n_levels']} levels")
            if v["counts"]:
                top = list(v["counts"].items())
                shown = ", ".join(f"{a} ({b:,})" for a, b in top[:40])
                L.append(f"  - {shown}" + (" ..." if len(top) > 40 else ""))
        L.append("")

    # ---- TSV columns -----------------------------------------------------
    for t in tsvs:
        L += [f"## `{Path(t['file']).name}`  ({t['n_rows']:,} rows x "
              f"{t['n_columns']} cols)", ""]
        if t["r_rownames_quirk"]:
            L.append("> Header has one fewer field than the data rows — this is "
                     "R `write.table(row.names=TRUE)`. The first data field is the "
                     "**row name** (shown below as `__rowname__`); read with "
                     "`index_col=0`.")
            L.append("")
        L.append(f"- all columns: `{t['columns']}`")
        L.append("")
        for k, v in t["columns_of_interest"].items():
            if v.get("levels") is not None and v["label_like"]:
                L.append(f"- **`{k}`** — {v['n_unique']} unique")
                top = list(v["counts"].items())
                L.append("  - " + ", ".join(f"{a} ({b:,})" for a, b in top[:40])
                         + (" ..." if len(top) > 40 else ""))
            else:
                L.append(f"- `{k}` — {v['n_unique']} unique "
                         f"(e.g. {v.get('levels_sample', v.get('levels', []))[:4]})")
        L.append("")

    # ---- cross-file vocabulary comparison --------------------------------
    L += ["## Vocabulary overlap", "",
          "Pairing needs ONE column present on both sides with matching level "
          "names. Jaccard is over level-name sets.", ""]
    pools = []
    for r in h5:
        for k, v in r["categorical"].items():
            if v["label_like"] and 2 <= v["n_levels"] <= 200:
                pools.append((Path(r["file"]).name, k, v["levels"]))
    for t in tsvs:
        for k, v in t["columns_of_interest"].items():
            if v["label_like"] and v.get("levels") and 2 <= v["n_unique"] <= 200:
                pools.append((Path(t["file"]).name, k, v["levels"]))

    best = []
    for i in range(len(pools)):
        for j in range(i + 1, len(pools)):
            f1, k1, l1 = pools[i]
            f2, k2, l2 = pools[j]
            if f1 == f2:
                continue
            js = jaccard(l1, l2)
            if js > 0:
                best.append((js, f1, k1, len(l1), f2, k2, len(l2),
                             sorted(set(l1) - set(l2))[:6],
                             sorted(set(l2) - set(l1))[:6]))
    best.sort(reverse=True)
    if not best:
        L.append("**No pair of label columns shares ANY level name.** "
                 "A mapping table is required.")
    for js, f1, k1, n1, f2, k2, n2, o1, o2 in best[:12]:
        flag = ("EXACT MATCH" if js == 1.0 else
                "near match" if js > 0.8 else "partial")
        L.append(f"- **{js:.2f}** {flag} — `{f1}:{k1}` ({n1}) vs `{f2}:{k2}` ({n2})")
        if o1:
            L.append(f"  - only in {k1}: {o1}")
        if o2:
            L.append(f"  - only in {k2}: {o2}")
    L.append("")

    # ---- id alignment: TSV rows vs each h5ad's obs_names ------------------
    L += ["## Do TSV rows align to cells?", ""]
    for t in tsvs:
        tid = set(t.get("candidate_ids") or [])
        if not tid:
            L.append(f"- `{Path(t['file']).name}`: no id-like column found; "
                     "rows may be in the same ORDER as the object "
                     "(check n_rows against n_obs).")
        for r in h5:
            hid = set(r.get("obs_ids") or [])
            if not hid:
                continue
            nm = Path(r['file']).name
            if tid:
                hit = len(hid & tid)
                L.append(f"- `{Path(t['file']).name}` vs `{nm}`: "
                         f"{hit}/{len(hid)} sampled obs_names found in the TSV"
                         + ("  **-> ids align**" if hit / max(len(hid), 1) > 0.99
                            else "  **-> ids DO NOT align**"))
            if r["n_obs"] == t["n_rows"]:
                L.append(f"  - row count equals `{nm}` n_obs "
                         f"({t['n_rows']:,}) — positional merge is plausible, "
                         "but verify the order before relying on it.")
    L.append("")

    md = "\n".join(L)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with open(f"{args.out}.md", "w") as fh:
        fh.write(md)
    with open(f"{args.out}.json", "w") as fh:
        json.dump({"h5ad": h5, "tsv": tsvs}, fh, indent=2, default=str)
    print()
    print(md)
    print(f"\nwrote {args.out}.md and {args.out}.json")


if __name__ == "__main__":
    main()
