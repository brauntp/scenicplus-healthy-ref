#!/usr/bin/env python3
"""
Inspect the GLUE-integrated healthy reference (RNA + ATAC AnnData/MuData)
and report everything SCENIC+ needs to know before the pipeline is built.

READ-ONLY. Opens .h5ad in backed mode; writes nothing except the report.

Usage
-----
    python inspect_anndata.py --rna RNA.h5ad --atac ATAC.h5ad --out report_py
    python inspect_anndata.py --mudata integrated.h5mu --out report_py

Emits <out>.json (machine readable) and <out>.md (human readable).
Both contain only shapes, names and summary statistics -- no count data.
"""
import argparse
import json
import sys
from pathlib import Path

try:
    import numpy as np
except ImportError:                                              # pragma: no cover
    sys.exit(
        "ERROR: this python has no 'numpy' (and so almost certainly no anndata).\n"
        f"       interpreter: {sys.executable}\n"
        "\n"
        "       Use the environment the GLUE integration ran in. To locate it:\n"
        "           bash 00_inspect/find_inspect_env.sh\n"
        "       then call that interpreter directly:\n"
        "           /path/to/env/bin/python 00_inspect/inspect_anndata.py ...")


# ---------------------------------------------------------------- helpers
def _jsonable(x):
    if isinstance(x, np.integer):
        return int(x)
    if isinstance(x, np.floating):
        return None if np.isnan(x) else float(x)
    if isinstance(x, np.bool_):
        return bool(x)
    if isinstance(x, (list, tuple)):
        return [_jsonable(i) for i in x]
    if isinstance(x, dict):
        return {str(k): _jsonable(v) for k, v in x.items()}
    return x


def describe_obs(obs, max_card=60):
    """Per-column dtype + cardinality. Flags columns usable as key_to_group_by."""
    out = {}
    for c in obs.columns:
        col = obs[c]
        rec = {"dtype": str(col.dtype)}
        try:
            n = int(col.nunique(dropna=True))
            rec["n_unique"] = n
            rec["n_missing"] = int(col.isna().sum())
            # A grouping key must be categorical-ish, low cardinality, complete.
            rec["candidate_group_key"] = bool(
                2 <= n <= max_card and rec["n_missing"] == 0
                and (str(col.dtype) in ("object", "category", "bool")
                     or str(col.dtype).startswith("str"))
            )
            if n <= max_card:
                vc = col.value_counts(dropna=True)
                rec["levels"] = {str(k): int(v) for k, v in vc.items()}
        except Exception as e:                                   # noqa: BLE001
            rec["error"] = repr(e)
        out[c] = rec
    return out


def matrix_report(X, name):
    """Raw counts or already normalized? Decided from the data, not the layer name."""
    if X is None:
        return {"which": name, "present": False}
    rec = {"which": name, "present": True, "type": type(X).__name__,
           "shape": list(getattr(X, "shape", [])),
           "dtype": str(getattr(X, "dtype", "?"))}
    try:
        n = X.shape[0]
        idx = np.linspace(0, n - 1, num=min(n, 200), dtype=int)
        sub = X[idx]
        sub = sub.toarray() if hasattr(sub, "toarray") else np.asarray(sub)
        sub = sub.astype(np.float64, copy=False)
        finite = sub[np.isfinite(sub)]
        if finite.size:
            rec["min"] = float(finite.min())
            rec["max"] = float(finite.max())
            rec["frac_zero"] = float((finite == 0).mean())
            rec["all_integer"] = bool(np.all(finite - np.floor(finite) == 0))
            # Integer-valued means counts REGARDLESS of magnitude: a sparse peak
            # matrix or a shallow RNA library can be all-integer with max < 10.
            # Only non-integer data can be log1p/normalized.
            if rec["all_integer"]:
                rec["looks_like"] = "raw_counts"
            elif rec["max"] < 30:
                rec["looks_like"] = "log1p_or_normalized"
            else:
                rec["looks_like"] = "normalized_nonlog"
            rec["n_rows_sampled"] = int(len(idx))
    except Exception as e:                                       # noqa: BLE001
        rec["error"] = repr(e)
    return rec


PEAK_SEPS = [(":", "-"), ("-", "-"), ("_", "_")]


def parse_regions(names):
    """Parse chr:start-end / chr-start-end / chr_start_end."""
    names = list(names)
    for s1, s2 in PEAK_SEPS:
        chroms, starts, ends = [], [], []
        for nm in names:
            try:
                c, rest = str(nm).split(s1, 1)
                a, b = rest.split(s2, 1)
                starts.append(int(a))
                ends.append(int(b))
                chroms.append(c)
            except Exception:                                    # noqa: BLE001
                continue
        if starts and len(starts) >= 0.9 * len(names):
            return {"format": f"chrom{s1}start{s2}end", "n_parsed": len(starts),
                    "chroms": chroms, "starts": np.asarray(starts),
                    "ends": np.asarray(ends)}
    return None


def region_report(var_names):
    """Peak width distribution, chr prefix, and genome-build evidence."""
    rec = {"n_regions": len(var_names),
           "examples": [str(v) for v in list(var_names)[:5]]}
    p = parse_regions(var_names)
    if p is None:
        rec["parsed"] = False
        rec["note"] = "var_names are not coordinate-like; peaks may be indexed by ID"
        return rec
    rec.update(parsed=True, format=p["format"], n_parsed=p["n_parsed"])
    w = (p["ends"] - p["starts"]).astype(float)
    uw = np.unique(w)
    rec["width"] = {"min": float(w.min()), "max": float(w.max()),
                    "median": float(np.median(w)), "mean": float(w.mean()),
                    "n_distinct": int(len(uw)),
                    "fixed_width": bool(len(uw) == 1),
                    "percentiles": {str(q): float(np.percentile(w, q))
                                    for q in (1, 25, 50, 75, 99)}}
    chroms = np.asarray([str(c) for c in p["chroms"]])
    uniq = sorted(set(chroms.tolist()))
    rec["chrom"] = {
        "n_unique": len(uniq),
        "has_chr_prefix": bool(all(c.startswith("chr") for c in uniq)),
        "examples": uniq[:30],
        "nonstandard": [c for c in uniq
                        if not c.replace("chr", "").isdigit()
                        and c.replace("chr", "") not in ("X", "Y", "M", "MT")][:20]}
    # hg38 chr1 = 248,956,422 ; hg19 chr1 = 249,250,621
    m = chroms == ("chr1" if rec["chrom"]["has_chr_prefix"] else "1")
    if m.any():
        rec["chrom"]["max_end_chr1"] = int(p["ends"][m].max())
        rec["chrom"]["build_note"] = ("compare to hg38 248956422 / hg19 249250621 "
                                      "-- peaks alone give a lower bound only")
    return rec


def var_report(var_names, var):
    rec = {"n_vars": len(var_names),
           "examples": [str(v) for v in list(var_names)[:5]],
           "columns": list(map(str, var.columns))[:60]}
    head = [str(v) for v in list(var_names)[:2000]]
    ens = sum(1 for v in head if v.startswith("ENSG"))
    rec["n_ensembl_prefixed_in_first_2000"] = int(ens)
    rec["looks_like"] = "ensembl_ids" if ens > 100 else "gene_symbols"
    return rec


def inspect_modality(ad, label, kind):
    rec = {"label": label, "kind": kind,
           "n_obs": int(ad.n_obs), "n_vars": int(ad.n_vars),
           "obs": describe_obs(ad.obs),
           "obsm_keys": list(map(str, ad.obsm.keys())),
           "obsm_shapes": {str(k): list(ad.obsm[k].shape) for k in ad.obsm.keys()},
           "layers": list(map(str, ad.layers.keys())),
           "obs_names_examples": [str(x) for x in ad.obs_names[:5]]}

    # The GLUE latent space -- what we pair cells on.
    glue = [k for k in ad.obsm.keys() if "glue" in str(k).lower()]
    rec["glue_obsm_keys"] = list(map(str, glue))
    rec["has_glue_latent"] = bool(glue)
    for k in glue:
        try:
            M = np.asarray(ad.obsm[k])
            rec.setdefault("glue_latent_stats", {})[str(k)] = {
                "shape": list(M.shape), "n_dim": int(M.shape[1]),
                "any_nan": bool(np.isnan(M).any()),
                "l2_norm_median": float(np.median(np.linalg.norm(M, axis=1)))}
        except Exception as e:                                   # noqa: BLE001
            rec.setdefault("glue_latent_stats", {})[str(k)] = {"error": repr(e)}

    # SCENIC+ reads .raw by default (use_raw_for_GEX_anndata=True).
    rec["X"] = matrix_report(ad.X, "X")
    try:
        rec["raw_present"] = ad.raw is not None
        if ad.raw is not None:
            rec["raw"] = matrix_report(ad.raw.X, "raw.X")
            rec["raw"]["n_vars"] = int(ad.raw.shape[1])
    except Exception as e:                                       # noqa: BLE001
        rec["raw_present"] = f"error: {e!r}"
    for lk in list(ad.layers.keys())[:6]:
        try:
            rec.setdefault("layer_reports", {})[str(lk)] = matrix_report(ad.layers[lk], lk)
        except Exception as e:                                   # noqa: BLE001
            rec.setdefault("layer_reports", {})[str(lk)] = {"error": repr(e)}

    rec["features"] = (region_report(ad.var_names) if kind == "atac"
                       else var_report(ad.var_names, ad.var))
    return rec


def shared_keys(rna_rec, atac_rec):
    """key_to_group_by must exist in BOTH modalities with overlapping level names."""
    out = []
    r_obs, a_obs = rna_rec.get("obs", {}), atac_rec.get("obs", {})
    for k in set(r_obs) & set(a_obs):
        r, a = r_obs[k], a_obs[k]
        if not (r.get("candidate_group_key") and a.get("candidate_group_key")):
            continue
        rl, al = set(r.get("levels", {})), set(a.get("levels", {}))
        if not rl or not al:
            continue
        shared = rl & al
        out.append({
            "key": k,
            "n_levels_rna": len(rl), "n_levels_atac": len(al),
            "n_shared_levels": len(shared),
            "frac_shared": round(len(shared) / max(len(rl | al), 1), 3),
            "only_in_rna": sorted(map(str, rl - al))[:12],
            "only_in_atac": sorted(map(str, al - rl))[:12],
            "min_cells_per_shared_level": (
                min(min(r["levels"][s], a["levels"][s]) for s in shared) if shared else 0),
            "shared_levels_under_50_cells": sorted(
                str(s) for s in shared
                if min(r["levels"][s], a["levels"][s]) < 50)[:20]})
    out.sort(key=lambda d: (-d["frac_shared"], -d["n_shared_levels"]))
    return out


# ---------------------------------------------------------------- markdown
def to_md(rep):
    v = rep["versions"]
    L = ["# SCENIC+ input inspection (Python side)", "",
         f"- anndata `{v.get('anndata')}`, mudata `{v.get('mudata')}`, "
         f"numpy `{v.get('numpy')}`, python `{v.get('python')}`", ""]
    for m in rep.get("modalities", []):
        if "error" in m:
            L += [f"## {m['label']}", "", f"> FAILED: {m['error']}", ""]
            continue
        L += [f"## {m['label']}  ({m['kind']})", "",
              f"- shape: **{m['n_obs']:,} cells x {m['n_vars']:,} features**",
              f"- obs_names look like: `{', '.join(m['obs_names_examples'][:3])}`",
              f"- GLUE latent present: **{m['has_glue_latent']}** "
              f"({', '.join(m['glue_obsm_keys']) or 'none'})",
              f"- obsm: {', '.join(f'{k}{tuple(s)}' for k, s in m['obsm_shapes'].items()) or 'none'}",
              f"- layers: {', '.join(m['layers']) or 'none'}",
              f"- `.raw` present: **{m.get('raw_present')}**"]
        x = m.get("X", {})
        if x.get("present"):
            L.append(f"- `X`: dtype {x.get('dtype')}, range [{x.get('min')}, {x.get('max')}], "
                     f"{100 * (x.get('frac_zero') or 0):.1f}% zero -> **{x.get('looks_like')}**")
        if "raw" in m:
            rx = m["raw"]
            L.append(f"- `.raw.X`: range [{rx.get('min')}, {rx.get('max')}] -> "
                     f"**{rx.get('looks_like')}**, n_vars {rx.get('n_vars')}")
        for lk, lr in (m.get("layer_reports") or {}).items():
            L.append(f"- layer `{lk}`: range [{lr.get('min')}, {lr.get('max')}] -> "
                     f"**{lr.get('looks_like')}**")
        L.append("")
        f = m.get("features", {})
        if m["kind"] == "atac":
            L += ["### Peak set", ""]
            if f.get("parsed"):
                w = f["width"]
                L += [f"- {f['n_regions']:,} regions, format `{f['format']}`, "
                      f"e.g. `{(f['examples'] or [''])[0]}`",
                      f"- width: median **{w['median']:.0f} bp**, "
                      f"min {w['min']:.0f}, max {w['max']:.0f}, "
                      f"{w['n_distinct']} distinct -> fixed width: **{w['fixed_width']}**",
                      f"- width percentiles: " + ", ".join(
                          f"p{q}={w['percentiles'][q]:.0f}" for q in ("1", "25", "50", "75", "99")),
                      f"- chroms: {f['chrom']['n_unique']} unique, `chr` prefix: "
                      f"**{f['chrom']['has_chr_prefix']}**",
                      f"- non-standard contigs: {', '.join(f['chrom']['nonstandard']) or 'none'}",
                      f"- max end on chr1: {f['chrom'].get('max_end_chr1')} "
                      f"({f['chrom'].get('build_note', '')})", ""]
            else:
                L += [f"- {f.get('n_regions')} features, NOT coordinate-parseable: "
                      f"`{f.get('examples')}`", ""]
        else:
            L += ["### Genes", "",
                  f"- {f.get('n_vars', 0):,} genes, identifiers look like "
                  f"**{f.get('looks_like')}** (e.g. `{', '.join(f.get('examples', [])[:3])}`)",
                  f"- var columns: {', '.join(f.get('columns', [])[:15]) or 'none'}", ""]
        cands = [k for k, d in m["obs"].items() if d.get("candidate_group_key")]
        L += [f"- candidate grouping columns in obs: "
              f"{', '.join(f'`{c}`' for c in cands) or 'NONE'}", ""]

    sk = rep.get("shared_group_keys")
    if sk is not None:
        L += ["## Shared grouping keys (candidates for `key_to_group_by`)", ""]
        if not sk:
            L.append("**None found.** RNA and ATAC obs share no low-cardinality, complete "
                     "column with overlapping level names. A harmonized cell-type label "
                     "must be added to both before SCENIC+ can group them.")
        for d in sk:
            L.append(f"- **`{d['key']}`** -- {d['n_shared_levels']} shared levels "
                     f"(RNA {d['n_levels_rna']} / ATAC {d['n_levels_atac']}, "
                     f"{d['frac_shared']:.0%} overlap); smallest shared level has "
                     f"{d['min_cells_per_shared_level']} cells")
            if d["only_in_rna"]:
                L.append(f"    - RNA-only levels: {', '.join(d['only_in_rna'])}")
            if d["only_in_atac"]:
                L.append(f"    - ATAC-only levels: {', '.join(d['only_in_atac'])}")
            if d["shared_levels_under_50_cells"]:
                L.append("    - under 50 cells (too small to metacell): "
                         f"{', '.join(d['shared_levels_under_50_cells'])}")
        L.append("")
    if rep.get("warnings"):
        L += ["## Warnings", ""] + [f"> {w}" for w in rep["warnings"]]
    return "\n".join(L)


# ---------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rna", type=Path, help="RNA .h5ad")
    ap.add_argument("--atac", type=Path, help="ATAC (peak) .h5ad")
    ap.add_argument("--mudata", type=Path, help=".h5mu holding both modalities")
    ap.add_argument("--out", type=Path, default=Path("report_py"))
    args = ap.parse_args()
    if not (args.mudata or args.rna or args.atac):
        ap.error("give --mudata, or --rna and/or --atac")

    try:
        import anndata
    except ImportError:
        sys.exit(
            "ERROR: this python has no 'anndata'.\n"
            f"       interpreter: {sys.executable}\n"
            "\n"
            "       You already have a python with anndata -- it is the env the GLUE\n"
            "       integration ran in. Find it with:\n"
            "           bash 00_inspect/find_inspect_env.sh\n"
            "       then call that interpreter directly, e.g.\n"
            "           /path/to/env/bin/python 00_inspect/inspect_anndata.py ...\n"
            "\n"
            "       (mudata is optional here -- only --mudata input needs it.)")

    rep = {"versions": {"python": sys.version.split()[0],
                        "anndata": anndata.__version__,
                        "numpy": np.__version__},
           "inputs": {}, "modalities": [], "warnings": []}
    try:
        import mudata
        rep["versions"]["mudata"] = mudata.__version__
    except Exception:                                            # noqa: BLE001
        rep["versions"]["mudata"] = None

    mods = {}
    if args.mudata:
        import mudata
        rep["inputs"]["mudata"] = str(args.mudata)
        md = mudata.read(str(args.mudata))
        rep["mudata_modalities"] = list(md.mod.keys())
        for k, ad in md.mod.items():
            kind = "atac" if any(t in k.lower() for t in ("atac", "acc", "peak")) else "rna"
            mods[k] = (ad, kind)
    if args.rna:
        rep["inputs"]["rna"] = str(args.rna)
        mods["scRNA"] = (anndata.read_h5ad(args.rna, backed="r"), "rna")
    if args.atac:
        rep["inputs"]["atac"] = str(args.atac)
        mods["scATAC"] = (anndata.read_h5ad(args.atac, backed="r"), "atac")

    for label, (ad, kind) in mods.items():
        try:
            rep["modalities"].append(inspect_modality(ad, label, kind))
        except Exception as e:                                   # noqa: BLE001
            rep["modalities"].append({"label": label, "kind": kind, "error": repr(e)})
            rep["warnings"].append(f"{label}: inspection failed: {e!r}")

    rna = next((m for m in rep["modalities"] if m.get("kind") == "rna"), None)
    atac = next((m for m in rep["modalities"] if m.get("kind") == "atac"), None)
    if rna and atac:
        rep["shared_group_keys"] = shared_keys(rna, atac)
        for lbl, m in (("RNA", rna), ("ATAC", atac)):
            if not m.get("has_glue_latent"):
                rep["warnings"].append(
                    f"No GLUE latent found in {lbl} obsm -- kNN pairing needs a shared embedding.")
        dims = {s["n_dim"] for m in (rna, atac)
                for s in (m.get("glue_latent_stats") or {}).values() if "n_dim" in s}
        if len(dims) > 1:
            rep["warnings"].append(
                f"GLUE latent dimensionality differs between modalities: {sorted(dims)}")
        if not atac.get("features", {}).get("parsed", True):
            rep["warnings"].append(
                "ATAC var_names are not coordinate-like; a peak->coordinate map is needed.")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with open(f"{args.out}.json", "w") as fh:
        json.dump(_jsonable(rep), fh, indent=2, default=str)
    with open(f"{args.out}.md", "w") as fh:
        fh.write(to_md(rep))
    print(f"wrote {args.out}.json and {args.out}.md")
    for w in rep["warnings"]:
        print("WARNING:", w)


if __name__ == "__main__":
    main()
