#!/usr/bin/env python3
"""
Export the peak-to-gene links as a portable, project-independent resource.

WHY THIS IS THE ASSET
---------------------
`region_to_gene_adj.tsv` is the **TF-agnostic** link table: every (peak, gene)
pair inside the TSS search space, scored by gradient boosting importance and by
correlation. It is not restricted to the 229 TFs that got eRegulons, and it does
not depend on any motif database -- so it is reusable in projects that have
nothing to do with this reference's TF list. That is the original goal of the
project: peak-to-gene regulatory linkages.

The eRegulon tables are a TF-filtered, motif-gated VIEW of this table. If you
want links, take these; if you want regulons, take those.

SCHEMA (read from scenicplus v1.0a2 `calculate_regions_to_genes_relationships`)
------------------------------------------------------------------------------
    target                  gene symbol
    region                  `chr:start-end`, the consensus peak
    importance              GBM feature importance of the region for that gene
    rho                     Spearman correlation (`--correlation-method SR`)
    importance_x_rho        importance * rho          (signed)
    importance_x_abs_rho    importance * |rho|        (magnitude)
    Distance                peak-to-TSS distance; present because `add_distance`
                            defaults True and the CLI does not override it
    + merged search-space columns (Chromosome, Start, End, Strand, Gene,
      Transcription_Start_Site, Transcript_type)

`rho` sign is the biologically interesting field: positive = accessibility and
expression move together (candidate enhancer), negative = they move oppositely
(candidate silencer, or a repressor-bound element).

OUTPUTS
-------
    peak_gene_links.parquet     everything, typed and coordinate-parsed
    peak_gene_links.bedpe       BEDPE: peak block <-> TSS block, for IGV/UCSC
    peak_gene_links_high.tsv    filtered set (see --min-abs-rho / --top-n-per-gene)
    peak_gene_links.README.md   provenance and caveats, written next to the data

MEMORY: chunked read, so a large adjacency never lands in memory whole.

Usage
-----
    python 05_report/export_peak_gene_links.py
    python 05_report/export_peak_gene_links.py --min-abs-rho 0.1 --top-n-per-gene 10
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

try:
    import numpy as np
    import pandas as pd
except ImportError as exc:
    sys.exit(f"ERROR: missing {exc.name}.\n"
             f"       interpreter: {sys.executable}\n"
             "       conda activate scplus-pairing")

REGION_RE = re.compile(r"^(?P<chrom>[^:]+):(?P<start>\d+)-(?P<end>\d+)$")

def write_table(df, out_dir, stem, index=False):
    """Parquet if pyarrow/fastparquet is present, else CSV. Returns the path.

    Do NOT make parquet a hard requirement: the pairing env is deliberately
    minimal (numpy, scipy, pandas, h5py, anndata, mudata) and does not ship an
    engine. 04_db/fetch_db_regions.py hit exactly this and already falls back --
    the first version of this script did not, and a completed extraction died on
    the write. Every consumer here reads either format.
    """
    from pathlib import Path
    p = Path(out_dir) / f"{stem}.parquet"
    try:
        df.to_parquet(p, index=index)
        return p
    except Exception as e:                                    # ImportError, ValueError
        p.unlink(missing_ok=True)
        alt = Path(out_dir) / f"{stem}.csv.gz"
        df.to_csv(alt, index=index, compression="gzip")
        print(f"    parquet unavailable ({type(e).__name__}: no pyarrow/"
              f"fastparquet) -- wrote {alt.name} instead")
        return alt



def parse_regions(regions: pd.Series) -> pd.DataFrame:
    """`chr1:1000-1500` -> columns. Fails loudly rather than silently dropping."""
    ex = regions.str.extract(REGION_RE)
    bad = ex["chrom"].isna()
    if bad.any():
        eg = regions[bad].unique()[:3].tolist()
        sys.exit(f"ERROR: {bad.sum():,} region names do not match "
                 f"'chrom:start-end', e.g. {eg}\n"
                 "       The exporter assumes the consensus-peak naming that\n"
                 "       aggregate_atac_sparse.py writes into scATAC var_names.")
    return pd.DataFrame({
        "peak_chrom": ex["chrom"].to_numpy(),
        "peak_start": ex["start"].astype(np.int64).to_numpy(),
        "peak_end": ex["end"].astype(np.int64).to_numpy(),
    }, index=regions.index)


def load_adjacency(path: Path, chunk: int) -> pd.DataFrame:
    if not path.exists():
        sys.exit(f"ERROR: adjacency not found: {path}\n"
                 "       It is the output of the region_to_gene rule; check\n"
                 "       output_data.region_to_gene_adjacencies in your config.")
    need = {"target", "region", "importance", "rho"}
    parts = []
    n = 0
    for i, ch in enumerate(pd.read_table(path, chunksize=chunk)):
        if i == 0:
            missing = need - set(ch.columns)
            if missing:
                sys.exit(f"ERROR: {path} lacks columns {sorted(missing)}\n"
                         f"       found: {sorted(ch.columns)}\n"
                         "       Schema is from scenicplus v1.0a2; another\n"
                         "       version may differ.")
            print(f"  columns: {', '.join(ch.columns)}")
        parts.append(ch)
        n += len(ch)
        if i and i % 20 == 0:
            print(f"    ... {n:,} rows")
    df = pd.concat(parts, ignore_index=True)
    print(f"  {n:,} links read")
    return df


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--adj", type=Path,
                    default=Path("03_pipeline/region_to_gene_adj.tsv"))
    ap.add_argument("--out-dir", type=Path, default=Path("05_report/plot_bundle"))
    ap.add_argument("--chunk", type=int, default=2_000_000)
    ap.add_argument("--min-abs-rho", type=float, default=0.05,
                    help="filter for the 'high' set; 0 disables")
    ap.add_argument("--top-n-per-gene", type=int, default=20,
                    help="per-gene cap in the 'high' set by importance; 0 disables")
    ap.add_argument("--no-bedpe", action="store_true")
    args = ap.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    print("=" * 74)
    print("peak-to-gene links")
    print("=" * 74)
    print("-- reading the adjacency (chunked) ---------------------------------")
    df = load_adjacency(args.adj, args.chunk)
    print()

    print("-- shape of the resource -------------------------------------------")
    n_link, n_gene, n_peak = len(df), df["target"].nunique(), df["region"].nunique()
    print(f"  links {n_link:,}   genes {n_gene:,}   peaks {n_peak:,}")
    print(f"  links per gene : median "
          f"{df.groupby('target').size().median():.0f}")
    print(f"  links per peak : median "
          f"{df.groupby('region').size().median():.0f}")
    # Count negatives DIRECTLY. `n_link - pos` buckets rho==0 and NaN into
    # "negative": on the real table that mislabelled 3,976 NaN links, inflating
    # the silencer count by ~1%.
    pos = int((df["rho"] > 0).sum())
    neg = int((df["rho"] < 0).sum())
    nan = int(df["rho"].isna().sum())
    zero = int((df["rho"] == 0).sum())
    print(f"  rho > 0 (candidate enhancer) : {pos:,} ({pos/n_link:.1%})")
    print(f"  rho < 0 (candidate silencer) : {neg:,} ({neg/n_link:.1%})")
    if nan or zero:
        print(f"  rho undefined / exactly 0    : {nan:,} NaN + {zero:,} zero "
              f"({(nan+zero)/n_link:.2%}) -- neither class")
    if "Distance" in df.columns:
        # Upstream ships this as a stringified 1-element list ("[-126154]"):
        # enhancer_to_gene.py leaves `result_df['Distance'].map(lambda x: x[0])`
        # COMMENTED OUT. A bare pd.to_numeric returns all-NaN and the summary
        # line silently printed "nan" -- observed on the real table.
        df["Distance"] = (df["Distance"].astype(str)
                            .str.strip("[]").replace("", np.nan)
                            .astype("Float64").astype("Int64"))
        d = df["Distance"].abs()
        print(f"  |Distance| to TSS : median {d.median():,.0f} bp, "
              f"max {d.max():,.0f} bp")
    print()

    print("-- writing --------------------------------------------------------")
    coords = parse_regions(df["region"])
    out = pd.concat([df, coords], axis=1)
    pq = write_table(out, args.out_dir, "peak_gene_links")
    print(f"  {pq.name:<30} {pq.stat().st_size/1024**2:>8.1f} MB  (all {n_link:,})")

    # -- filtered set --------------------------------------------------------
    hi = out
    if args.min_abs_rho > 0:
        hi = hi[hi["rho"].abs() >= args.min_abs_rho]
    if args.top_n_per_gene > 0:
        hi = (hi.sort_values("importance", ascending=False)
                .groupby("target", sort=False)
                .head(args.top_n_per_gene))
    hi = hi.sort_values(["peak_chrom", "peak_start"])
    ht = args.out_dir / "peak_gene_links_high.tsv"
    hi.to_csv(ht, sep="\t", index=False)
    print(f"  {ht.name:<30} {ht.stat().st_size/1024**2:>8.1f} MB  "
          f"({len(hi):,} links, {hi['target'].nunique():,} genes)")

    # -- BEDPE ---------------------------------------------------------------
    bedpe_note = "- `peak_gene_links.bedpe` — not produced (`--no-bedpe`)."
    if not args.no_bedpe:
        tss_col = ("Transcription_Start_Site"
                   if "Transcription_Start_Site" in hi.columns else None)
        if tss_col is None:
            print("  BEDPE skipped: the adjacency has no "
                  "Transcription_Start_Site column.")
            print("    `Distance` cannot substitute: it is strand-RELATIVE")
            print("    (location x min_distance, where location is multiplied by")
            print("    the gene's strand), so its sign gives upstream/downstream")
            print("    in gene orientation, not genomic direction. Without")
            print("    `Strand` the TSS is genuinely unrecoverable -- join")
            print("    03_pipeline/genome_annotation.tsv on the gene to get it.")
            bedpe_note = (
                "- `peak_gene_links.bedpe` — **NOT produced.** The adjacency\n"
                "  lacked `Transcription_Start_Site`, and `Distance` cannot\n"
                "  substitute: its sign is strand-relative, so TSS position is\n"
                "  unrecoverable without `Strand`. Join\n"
                "  `03_pipeline/genome_annotation.tsv` on the gene to build one.")
        else:
            tss = pd.to_numeric(hi[tss_col], errors="coerce")
            keep = tss.notna()
            bp = pd.DataFrame({
                "chrom1": hi.loc[keep, "peak_chrom"],
                "start1": hi.loc[keep, "peak_start"],
                "end1": hi.loc[keep, "peak_end"],
                "chrom2": hi.loc[keep, "peak_chrom"],
                "start2": tss[keep].astype(np.int64),
                "end2": tss[keep].astype(np.int64) + 1,
                "name": hi.loc[keep, "target"],
                # BEDPE score is an integer 0-1000; scale |rho| into it and say so.
                "score": (hi.loc[keep, "rho"].abs().clip(0, 1) * 1000)
                           .round().astype(int),
                "strand1": ".",
                "strand2": hi.loc[keep, "Strand"] if "Strand" in hi.columns else ".",
            })
            bpe = args.out_dir / "peak_gene_links.bedpe"
            bp.to_csv(bpe, sep="\t", index=False, header=False)
            print(f"  {bpe.name:<30} {bpe.stat().st_size/1024**2:>8.1f} MB  "
                  f"({len(bp):,} pairs)")
            bedpe_note = (
                "- `peak_gene_links.bedpe` — peak block to TSS block, for IGV /\n"
                "  UCSC. **Its score column is `|rho| * 1000` rounded**, because\n"
                "  BEDPE requires an integer 0–1000 — that discards the sign, so\n"
                "  read the sign from the table above.")

    # -- provenance, written beside the data ---------------------------------
    readme = args.out_dir / "peak_gene_links.README.md"
    readme.write_text(f"""# Peak-to-gene links

Produced by `05_report/export_peak_gene_links.py` from
`{args.adj}` (SCENIC+ v1.0a2 `region_to_gene` rule).

## What these are

Every (peak, gene) pair inside the TSS search space, scored by GBM importance and
Spearman correlation across the paired metacells. **TF-agnostic** — no
motif database involved, not restricted to the TFs that got eRegulons. Reusable
outside this project.

| | |
|---|---|
| links | {n_link:,} |
| genes | {n_gene:,} |
| peaks | {n_peak:,} |
| rho > 0 | {pos:,} ({pos/n_link:.1%}) |
| rho < 0 | {neg:,} ({neg/n_link:.1%}) |
| rho NaN or exactly 0 | {nan + zero:,} ({(nan+zero)/n_link:.2%}) — neither class |

`rho > 0` = accessibility and expression covary (candidate enhancer);
`rho < 0` = they move oppositely (candidate silencer or repressor-bound element).

## Files

- `{pq.name}` — all links, coordinates parsed out of the region name
- `peak_gene_links_high.tsv` — `|rho| >= {args.min_abs_rho}`, top
  {args.top_n_per_gene} per gene by importance
{bedpe_note}

## Caveats that travel with the data

1. **These are correlations across metacells, not measured contacts.** No 3C
   evidence. A link means accessibility and expression covary within the search
   space, which is consistent with regulation but does not establish it.
2. **The metacells are GLUE-paired, not true multiome.** RNA and ATAC come from
   different cells matched through a shared embedding, so every link inherits the
   pairing's error. Pairing QC flagged both erythroid groups as the worst
   cross-modal alignments — treat erythroid-specific links with more caution.
3. **The search space bounds what can be found**: ±150 kb from TSS, capped by
   chromosome boundaries. Longer-range interactions are absent by construction,
   not by evidence.
4. **Aggregation limits resolution.** Metacells average ~50 cells, so a link
   present in a small subpopulation inside a cell type will be diluted.
5. **Peaks are the reference's own consensus set** (393,832 peaks). To use these
   links against another dataset the peaks must be re-mapped, and partial overlap
   will lose links.
""")
    print(f"  {readme.name:<30} {readme.stat().st_size/1024:>8.1f} KB")
    print()
    print(f"  bundle dir: {args.out_dir}")
    print("=" * 74)


if __name__ == "__main__":
    main()
