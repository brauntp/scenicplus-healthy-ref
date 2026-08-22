#!/usr/bin/env python3
"""
Gate the hand-built paired MuData before a multi-hour SCENIC+ run on ARC.

Our data are unpaired scRNA + scATAC co-embedded with scGLUE, so
02_pair/glue_metacells.py builds the paired MuData itself rather than letting
SCENIC+ do it (see the header of 03_pipeline/config.template.yaml for why:
`process_non_multiome_data` pairs RNA and ATAC metacells by independent uniform
random draws within a categorical label, which destroys the within-cell-type
covariation that region-to-gene is supposed to measure). The consequence is that
nothing in SCENIC+ ever validates our file -- the pipeline just indexes into it
and assumes it is well formed. This script is that missing check.

Every requirement below is taken from the pinned scenicplus v1.0a2 source, not
from documentation:

  * modality keys "scRNA" and "scATAC", exactly, as string literals:
      cli/commands.py:549,558     prepare_menr -> cistromes + TF names
      cli/commands.py:662,663     get_search_space
      cli/commands.py:691         infer_TF_to_gene
      cli/commands.py:740,741     infer_region_to_gene
      cli/commands.py:958,959     AUCell (score_eRegulons)
      scenicplus_mudata.py:35,36  final ScenicPlusMuData assembly
    There is no fallback, no alias, and no case-insensitive lookup anywhere.

  * shared obs across modalities: enhancer_to_gene.py:216-217 takes
      EXP = df_exp_mtx[gene_names].to_numpy(); ACC = df_acc_mtx.to_numpy()
    and regresses column-of-ACC on column-of-EXP positionally, ROW BY ROW.
    MuData does not force the two AnnData objects to share row order, so
    mismatched obs_names would silently regress metacell i's expression on
    metacell j's accessibility.

  * no duplicate var_names: enhancer_to_gene.py:199-202 raises
    "Expression matrix contains duplicate gene names" / "Chromatin
    accessibility matrix contains duplicate gene names" -- but only once
    region_to_gene starts, i.e. hours in. We check up front.

  * region var_names parseable as chrom:start-end: pycisTopic's
    region_names_to_coordinates (mirrored at utils.py:215-224) does
      chrom = i.split(':', 1)[0]; coor = i.split(':', 1)[1]
      start = int(coor.split('-', 1)[0]); end = int(coor.split('-', 1)[1])
    and -- critically -- it FILTERS on `if ':' in i`. A region whose name lacks
    a colon is dropped from the search space WITHOUT WARNING; a region with a
    colon but a non-integer coordinate raises ValueError deep inside
    get_search_space. Both are checked here.

  * finite values: the GBM regressors and the Spearman correlation will happily
    consume NaN/inf and produce garbage or opaque failures.

Exit code 0 = safe to submit. Nonzero = do not submit.

Usage:
    python 03_pipeline/validate_h5mu.py <path/to/ACC_GEX.h5mu> [--max-report N]
"""

from __future__ import annotations

import argparse
import sys
from typing import List

# Literal modality keys required by scenicplus v1.0a2. Do not "fix" these.
RNA_MOD = "scRNA"
ACC_MOD = "scATAC"

# ANSI is only emitted when stdout is a tty, so SLURM logs stay clean.
_TTY = sys.stdout.isatty()


def _c(code: str, s: str) -> str:
    return f"\033[{code}m{s}\033[0m" if _TTY else s


class Report:
    """Accumulates failures instead of aborting, so one run surfaces every problem."""

    def __init__(self) -> None:
        self.failures: List[str] = []
        self.warnings: List[str] = []

    def check(self, ok: bool, label: str, detail: str = "") -> bool:
        if ok:
            print(f"  {_c('32', 'PASS')}  {label}")
        else:
            print(f"  {_c('31', 'FAIL')}  {label}")
            if detail:
                for line in detail.rstrip().splitlines():
                    print(f"         {line}")
            self.failures.append(label)
        return ok

    def warn(self, label: str, detail: str = "") -> None:
        print(f"  {_c('33', 'WARN')}  {label}")
        if detail:
            for line in detail.rstrip().splitlines():
                print(f"         {line}")
        self.warnings.append(label)

    def info(self, label: str) -> None:
        print(f"  {_c('36', 'INFO')}  {label}")


def _fmt_examples(items, n: int) -> str:
    items = list(items)
    shown = ", ".join(repr(x) for x in items[:n])
    if len(items) > n:
        shown += f", ... (+{len(items) - n} more)"
    return shown


def _nonfinite_count(X) -> int:
    """Count non-finite entries in a dense or sparse matrix without densifying."""
    import numpy as np
    import scipy.sparse as sp

    if sp.issparse(X):
        data = X.data
    else:
        data = np.asarray(X).ravel()
    if data.dtype.kind not in "fc":
        # Integer / boolean matrices cannot hold NaN or inf.
        return 0
    return int((~np.isfinite(data)).sum())


def main(argv=None) -> int:
    p = argparse.ArgumentParser(
        description="Validate the hand-built paired MuData for SCENIC+ v1.0a2.",
    )
    p.add_argument("h5mu", help="Path to the paired MuData (.h5mu) built by 02_pair/glue_metacells.py")
    p.add_argument(
        "--max-report",
        type=int,
        default=5,
        metavar="N",
        help="How many offending names to print per failure (default: 5).",
    )
    args = p.parse_args(argv)

    import os

    print("=" * 74)
    print("SCENIC+ paired-MuData validation")
    print("=" * 74)
    print(f"file: {args.h5mu}")

    if not os.path.exists(args.h5mu):
        print(f"\n{_c('31', 'FATAL')}  file does not exist: {args.h5mu}")
        print("        02_pair/glue_metacells.py must write the paired MuData to the")
        print("        path given by output_data.combined_GEX_ACC_mudata in the config,")
        print("        BEFORE snakemake is invoked.")
        return 2

    try:
        import mudata
        import numpy as np  # noqa: F401  (used by helpers)
        import pandas as pd
    except ImportError as e:
        print(f"\n{_c('31', 'FATAL')}  missing dependency: {e}")
        print("        Activate the scenicplus conda env (03_pipeline/environment.yml).")
        return 2

    try:
        # Older mudata (0.2.3, our pin) exposes read(); it warns on some files
        # but returns a MuData regardless.
        mdata = mudata.read(args.h5mu)
    except Exception as e:  # noqa: BLE001 -- any read failure is fatal and we want the message
        print(f"\n{_c('31', 'FATAL')}  could not read as MuData: {type(e).__name__}: {e}")
        return 2

    r = Report()

    # ---------------------------------------------------------------- modalities
    print("\n[1] modality keys (literal, required by scenicplus v1.0a2)")
    present = list(mdata.mod.keys())
    r.info(f"modalities found: {present}")

    have_rna = r.check(
        RNA_MOD in mdata.mod,
        f'modality "{RNA_MOD}" present',
        detail=(
            f"found {present} instead.\n"
            f'scenicplus indexes mdata["{RNA_MOD}"] as a bare string literal in\n'
            "cli/commands.py:558,691,740,958 and scenicplus_mudata.py:35.\n"
            "Rename the modality; there is no configurable alias."
        ),
    )
    have_acc = r.check(
        ACC_MOD in mdata.mod,
        f'modality "{ACC_MOD}" present',
        detail=(
            f"found {present} instead.\n"
            f'scenicplus indexes mdata["{ACC_MOD}"] in\n'
            "cli/commands.py:549,662,741,959 and scenicplus_mudata.py:36."
        ),
    )
    extra = [m for m in present if m not in (RNA_MOD, ACC_MOD)]
    if extra:
        r.warn(
            f"unexpected extra modalities: {extra}",
            "Harmless for the pipeline (they are never indexed), but they are carried\n"
            "through to the final object and inflate file size.",
        )

    if not (have_rna and have_acc):
        print("\n" + "=" * 74)
        print(f"{_c('31', 'RESULT: FAIL')} -- modality keys wrong; nothing else can be checked.")
        print("=" * 74)
        return 1

    rna = mdata[RNA_MOD]
    acc = mdata[ACC_MOD]

    # -------------------------------------------------------------------- shapes
    print("\n[2] shapes")
    r.info(f"{RNA_MOD}:  {rna.n_obs} metacells x {rna.n_vars} genes")
    r.info(f"{ACC_MOD}: {acc.n_obs} metacells x {acc.n_vars} regions")
    r.info(f"MuData:  {mdata.n_obs} obs x {mdata.n_vars} vars (union)")
    r.check(rna.n_obs > 0 and rna.n_vars > 0, f"{RNA_MOD} is non-empty")
    r.check(acc.n_obs > 0 and acc.n_vars > 0, f"{ACC_MOD} is non-empty")

    # ------------------------------------------------------------------ obs align
    print("\n[3] obs_names identical across modalities (positional pairing)")
    rna_obs = list(rna.obs_names)
    acc_obs = list(acc.obs_names)

    same_len = r.check(
        len(rna_obs) == len(acc_obs),
        "same number of metacells in both modalities",
        detail=f"{RNA_MOD}={len(rna_obs)}  {ACC_MOD}={len(acc_obs)}",
    )
    if same_len:
        identical = rna_obs == acc_obs
        if identical:
            r.check(True, "obs_names identical AND in the same order")
        else:
            as_sets_equal = set(rna_obs) == set(acc_obs)
            if as_sets_equal:
                first_bad = next(
                    i for i, (a, b) in enumerate(zip(rna_obs, acc_obs)) if a != b
                )
                r.check(
                    False,
                    "obs_names in the same ORDER",
                    detail=(
                        "the two modalities contain the same metacell names but in a\n"
                        f"different order (first divergence at row {first_bad}: "
                        f"{rna_obs[first_bad]!r} vs {acc_obs[first_bad]!r}).\n"
                        "enhancer_to_gene.py:216-217 converts both to numpy and regresses\n"
                        "row-for-row, so this would pair metacell i's expression with\n"
                        "metacell j's accessibility. Reindex ATAC to the RNA order."
                    ),
                )
            else:
                only_rna = [x for x in rna_obs if x not in set(acc_obs)]
                only_acc = [x for x in acc_obs if x not in set(rna_obs)]
                r.check(
                    False,
                    "obs_names identical across modalities",
                    detail=(
                        f"only in {RNA_MOD} ({len(only_rna)}): {_fmt_examples(only_rna, args.max_report)}\n"
                        f"only in {ACC_MOD} ({len(only_acc)}): {_fmt_examples(only_acc, args.max_report)}"
                    ),
                )

    dup_rna_obs = pd.Index(rna_obs).duplicated()
    dup_acc_obs = pd.Index(acc_obs).duplicated()
    r.check(
        not dup_rna_obs.any(),
        f"{RNA_MOD} obs_names unique",
        detail=_fmt_examples(pd.Index(rna_obs)[dup_rna_obs], args.max_report),
    )
    r.check(
        not dup_acc_obs.any(),
        f"{ACC_MOD} obs_names unique",
        detail=_fmt_examples(pd.Index(acc_obs)[dup_acc_obs], args.max_report),
    )

    # ----------------------------------------------------------------- gene names
    print("\n[4] gene names (scRNA.var_names)")
    genes = pd.Index(rna.var_names.astype(str))
    empty_genes = [g for g in genes if g.strip() == "" or g.lower() in ("nan", "none")]
    r.check(
        len(empty_genes) == 0,
        "no empty / placeholder gene names",
        detail=f"{len(empty_genes)} offending: {_fmt_examples(empty_genes, args.max_report)}",
    )
    dup_genes = genes[genes.duplicated()]
    r.check(
        len(dup_genes) == 0,
        "gene names unique",
        detail=(
            f"{len(dup_genes)} duplicated: {_fmt_examples(dup_genes.unique(), args.max_report)}\n"
            "enhancer_to_gene.py:199-200 raises on this, but only after region_to_gene\n"
            "has already started. Deduplicate now (e.g. var_names_make_unique() is NOT\n"
            "enough -- collapse or drop, since a suffixed duplicate will not match the\n"
            "Ensembl gene annotation and silently drops out of the search space)."
        ),
    )
    ws_genes = [g for g in genes if g != g.strip()]
    if ws_genes:
        r.warn(
            f"{len(ws_genes)} gene names have leading/trailing whitespace",
            f"{_fmt_examples(ws_genes, args.max_report)}\n"
            "These will not match the biomart gene annotation and will be dropped\n"
            "from the search space without an error.",
        )

    # --------------------------------------------------------------- region names
    print("\n[5] region names (scATAC.var_names) parseable as chrom:start-end")
    regions = pd.Index(acc.var_names.astype(str))

    dup_regions = regions[regions.duplicated()]
    r.check(
        len(dup_regions) == 0,
        "region names unique",
        detail=(
            f"{len(dup_regions)} duplicated: {_fmt_examples(dup_regions.unique(), args.max_report)}\n"
            "enhancer_to_gene.py:201-202 raises on this during region_to_gene."
        ),
    )

    no_colon = [x for x in regions if ":" not in x]
    r.check(
        len(no_colon) == 0,
        "every region name contains ':'",
        detail=(
            f"{len(no_colon)} offending: {_fmt_examples(no_colon, args.max_report)}\n"
            "region_names_to_coordinates (pycisTopic.utils; mirrored at\n"
            "scenicplus/utils.py:215-224) builds its dataframe with a\n"
            "`for i in region_names if ':' in i` filter -- colon-less regions are\n"
            "SILENTLY DROPPED from the search space, so this never errors, it just\n"
            "quietly shrinks the analysis."
        ),
    )

    bad_coord = []
    neg_or_inverted = []
    for x in regions:
        if ":" not in x:
            continue
        coor = x.split(":", 1)[1]
        if "-" not in coor:
            bad_coord.append(x)
            continue
        s, e = coor.split("-", 1)
        try:
            si, ei = int(s), int(e)
        except ValueError:
            bad_coord.append(x)
            continue
        if si < 0 or ei <= si:
            neg_or_inverted.append(x)

    r.check(
        len(bad_coord) == 0,
        "coordinates parse as int(start) and int(end)",
        detail=(
            f"{len(bad_coord)} offending: {_fmt_examples(bad_coord, args.max_report)}\n"
            "utils.py:219-220 calls int() on these directly -- a non-integer raises\n"
            "ValueError inside get_search_space. Note that names like\n"
            "'chr1:1000-2000_peak3' or 'chr1-1000-2000' both land here."
        ),
    )
    r.check(
        len(neg_or_inverted) == 0,
        "coordinates are non-negative with end > start",
        detail=(
            f"{len(neg_or_inverted)} offending: {_fmt_examples(neg_or_inverted, args.max_report)}\n"
            "PyRanges will accept these and then behave unpredictably in the\n"
            "search-space overlap join."
        ),
    )

    chroms = sorted({x.split(":", 1)[0] for x in regions if ":" in x})
    r.info(f"{len(chroms)} distinct chromosome names; first few: {chroms[:6]}")
    if chroms and not any(c.startswith("chr") for c in chroms):
        r.warn(
            "no chromosome name starts with 'chr' (Ensembl-style names?)",
            "download_genome_annotations is invoked WITHOUT\n"
            "--do_not_use_ucsc_chromosome_style, so the Snakefile always produces\n"
            "UCSC-style ('chr1') annotation and chromsizes. Ensembl-style region\n"
            "names ('1') will not join against it and the search space comes back\n"
            "empty or near-empty. Rename regions to UCSC style.",
        )

    # ---------------------------------------------------------------- finiteness
    print("\n[6] matrix values finite")
    n_bad_rna = _nonfinite_count(rna.X)
    n_bad_acc = _nonfinite_count(acc.X)
    r.check(
        n_bad_rna == 0,
        f"{RNA_MOD}.X has no NaN/inf",
        detail=f"{n_bad_rna} non-finite entries",
    )
    r.check(
        n_bad_acc == 0,
        f"{ACC_MOD}.X has no NaN/inf",
        detail=f"{n_bad_acc} non-finite entries",
    )

    # NaN in obs is not fatal for the pipeline but usually signals a botched
    # metacell merge, which is worth knowing before a 12-hour job.
    for name, ad in ((RNA_MOD, rna), (ACC_MOD, acc)):
        if ad.obs.shape[1]:
            nan_cols = [c for c in ad.obs.columns if ad.obs[c].isna().any()]
            if nan_cols:
                r.warn(
                    f"{name}.obs has NaN in: {nan_cols[: args.max_report]}",
                    "Not read by the pipeline, but usually a sign the metacell join\n"
                    "in 02_pair/glue_metacells.py lost rows.",
                )

    # ------------------------------------------------------------------- summary
    print("\n" + "=" * 74)
    if r.failures:
        print(f"{_c('31', 'RESULT: FAIL')} -- {len(r.failures)} check(s) failed:")
        for f in r.failures:
            print(f"  - {f}")
        if r.warnings:
            print(f"({len(r.warnings)} warning(s) also raised.)")
        print("Do NOT submit the cluster job. Fix 02_pair/glue_metacells.py output first.")
        print("=" * 74)
        return 1

    print(f"{_c('32', 'RESULT: PASS')} -- MuData is safe to hand to SCENIC+.")
    if r.warnings:
        print(f"{len(r.warnings)} warning(s) -- review above; none are hard blockers.")
    print(
        f"Place this file at output_data.combined_GEX_ACC_mudata so Snakemake treats\n"
        f"prepare_GEX_ACC as already satisfied."
    )
    print("=" * 74)
    return 0


if __name__ == "__main__":
    sys.exit(main())
