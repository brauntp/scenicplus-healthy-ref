#!/usr/bin/env python3
"""
Build genome_annotation.tsv and chromsizes.tsv, so the DAG never calls BioMart.

WHY THIS EXISTS
---------------
`download_genome_annotations` failed with:

    xml.etree.ElementTree.ParseError: mismatched tag: line 62, column 2

That is not our bug and not a pin problem: pybiomart asked
`http://www.ensembl.org` for its dataset configuration and got malformed XML
back. Ensembl's BioMart returns HTML error pages, truncated responses and
maintenance notices under load, and `pybiomart` feeds the body straight to an
XML parser. Worse, scenicplus's `download_gene_annotation_and_chromsizes` then
makes a SECOND network call chain -- NCBI esearch plus an assembly report -- to
derive chromosome sizes, with its own retry ceiling.

So the rule is two flaky third-party services standing between a 24-hour job and
its first real computation, for two small static files. This script builds them
from UCSC's API -- one JSON call per chromosome, whole-chromosome responses in
well under a second -- and writes them in exactly the format the Snakefile's
downstream rules read. Once they exist, snakemake treats
`download_genome_annotations` as satisfied and never runs it.

FORMATS (read from the pinned v1.0a2 source, not guessed)
---------------------------------------------------------
`data_wrangling/gene_search_space.py::download_gene_annotation_and_chromsizes`
returns, and `get_search_space` consumes:

    genome_annotation : Chromosome, Start, End, Strand, Gene,
                        Transcription_Start_Site, Transcript_type
                        - Strand as '+'/'-' (NOT 1/-1)
                        - filtered to Transcript_type == "protein_coding"
                        - Chromosome in UCSC style ("chr1") when
                          use_ucsc_chromosome_style=True, which is the default
    chromsizes        : Chromosome, Start, End    (Start always 0)

Usage:
    python 03_pipeline/build_genome_annotation.py \
        --out-dir 03_pipeline \
        --genome hg38

    # then verify the DAG no longer wants the rule:
    bash 03_pipeline/run_pipeline.sh --config 03_pipeline/config.yaml \
        --cores 4 -n
"""
from __future__ import annotations

import argparse
import json
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

# UCSC's API, not Ensembl REST. Ensembl's /overlap/region endpoint caps requests
# at 5 Mb ("46709983 is greater than the maximum allowed length of 5000000"), so
# a whole genome would be ~620 paginated calls -- more moving parts than the
# BioMart path this replaces. UCSC returns a whole chromosome's knownGene track
# in one call, in well under a second.
UCSC = "https://api.genome.ucsc.edu"

# The primary assembly only. Scaffolds and patches are excluded exactly as the
# upstream code excludes anything outside the assembled molecules, and a TSS on
# an unplaced scaffold cannot anchor a search space against our peaks anyway.
MAIN_CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"]


def get(url: str, tries: int = 4, pause: float = 1.5) -> dict:
    """GET a JSON endpoint, retrying transient failures.

    Unlike the BioMart XML path this replaces, a failure here surfaces as an
    HTTP status or a JSON decode error on a small body -- not as a parse error
    sixty lines into a malformed document -- so a retry loop is meaningful.
    """
    last = None
    for attempt in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=300) as fh:
                return json.loads(fh.read().decode())
        except (urllib.error.HTTPError, urllib.error.URLError,
                TimeoutError, json.JSONDecodeError) as exc:
            last = exc
            if attempt < tries - 1:
                time.sleep(pause * (attempt + 1))
    sys.exit(f"ERROR: {url}\n       failed after {tries} attempts: {last!r}")


def fetch_chromsizes(genome: str) -> list[tuple[str, int, int]]:
    """Chromosome, Start(0), End for the primary assembly, UCSC-style names."""
    d = get(f"{UCSC}/list/chromosomes?genome={genome}")
    lengths = d.get("chromosomes", {})
    if not lengths:
        sys.exit("ERROR: no 'chromosomes' key in the response; the API's "
                 "response shape may have changed.")
    rows = []
    for c in MAIN_CHROMS:
        if c in lengths:
            rows.append((c, 0, int(lengths[c])))
    missing = [c for c in MAIN_CHROMS if c not in lengths]
    if missing:
        sys.exit(f"ERROR: primary-assembly chromosomes absent from {genome}: "
                 f"{missing}")
    return rows


def fetch_annotation(genome: str, track: str,
                     chroms: list[str]) -> list[tuple]:
    """One row per protein-coding TRANSCRIPT, carrying that transcript's TSS.

    Not per gene: the search space is anchored on transcription start sites, and
    a gene with several protein-coding transcripts has several distinct TSSs.
    Collapsing to one per gene would silently narrow every search space.

    Field semantics were verified against the API rather than assumed:
    `geneName` is the HGNC symbol (RUNX1, APP, SOD1 ... all present on chr21),
    while `geneName2` holds UniProt accessions and is NOT what the Snakefile's
    downstream rules match against our RNA var_names.
    """
    rows = []
    for c in chroms:
        d = get(f"{UCSC}/getData/track?genome={genome};track={track};chrom={c}")
        recs = d.get(track) or d.get(c) or []
        kept = 0
        for t in recs:
            if t.get("transcriptType") != "protein_coding":
                continue
            sym = t.get("geneName")
            if not sym or sym == "none":
                continue
            strand = t.get("strand")
            if strand not in ("+", "-"):
                continue
            start, end = int(t["chromStart"]), int(t["chromEnd"])
            # chromStart is 0-based (BED); the TSS is the 5' end, which is
            # chromStart on + and chromEnd on -.
            tss = start if strand == "+" else end
            rows.append((c, start, end, strand, sym, tss, "protein_coding"))
            kept += 1
        print(f"  {c:<6} protein-coding transcripts: {kept:>6,}"
              f"   (of {len(recs):>6,} in track)", flush=True)
    if not rows:
        sys.exit(f"ERROR: no protein-coding transcripts found in track "
                 f"'{track}'. Check the track name against "
                 f"{UCSC}/list/tracks?genome={genome}")
    return rows


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out-dir", type=Path, required=True,
                    help="Where to write the two TSVs. Must be snakemake's "
                         "working directory -- the Snakefile names them "
                         "relative to --directory, i.e. 03_pipeline/.")
    ap.add_argument("--genome", default="hg38",
                    help="UCSC assembly name. Must match the assembly the "
                         "peaks and the cisTarget database use.")
    ap.add_argument("--track", default="knownGene",
                    help="UCSC gene track. knownGene is GENCODE-backed; a "
                         "versioned alias (knownGeneV47 etc.) pins the GENCODE "
                         "release if reproducibility across UCSC updates "
                         "matters.")
    ap.add_argument("--annotation-name", default="genome_annotation.tsv")
    ap.add_argument("--chromsizes-name", default="chromsizes.tsv")
    ap.add_argument("--check-against", type=Path, default=None,
                    help="Paired .h5mu. Reports how many of its scRNA "
                         "var_names match the Gene column. get_search_space "
                         "does `gene_annotation.query(\"Gene in "
                         "@scplus_genes\")` -- a plain intersection that "
                         "returns EMPTY without erroring if the two use "
                         "different symbol namespaces, so a low overlap here "
                         "is the difference between a result and silence.")
    ap.add_argument("--force", action="store_true",
                    help="Overwrite existing files. Without this an existing "
                         "file is left alone, since it may be a good one from "
                         "a successful BioMart call.")
    args = ap.parse_args()

    ann_p = args.out_dir / args.annotation_name
    chr_p = args.out_dir / args.chromsizes_name
    args.out_dir.mkdir(parents=True, exist_ok=True)

    for p in (ann_p, chr_p):
        if p.exists() and not args.force:
            sys.exit(f"ERROR: {p} exists. Pass --force to overwrite, or delete "
                     f"it first.\n       (An existing file may be a valid one "
                     f"from a successful run.)")

    print("=" * 70)
    print("building genome annotation + chromsizes from the Ensembl REST API")
    print("=" * 70)
    print(f"  genome  : {args.genome}   track: {args.track}")
    print(f"  out dir : {args.out_dir}")
    print()

    print("-- chromosome sizes ------------------------------------------------")
    cs = fetch_chromsizes(args.genome)
    with open(chr_p, "w") as fh:
        fh.write("Chromosome\tStart\tEnd\n")
        for c, s, e in cs:
            fh.write(f"{c}\t{s}\t{e}\n")
    print(f"  {len(cs)} chromosomes -> {chr_p}")
    print()

    print("-- gene annotation (protein-coding transcripts) --------------------")
    ann = fetch_annotation(args.genome, args.track,
                           [c for c, _, _ in cs])
    with open(ann_p, "w") as fh:
        fh.write("Chromosome\tStart\tEnd\tStrand\tGene\t"
                 "Transcription_Start_Site\tTranscript_type\n")
        for r in ann:
            fh.write("\t".join(str(x) for x in r) + "\n")
    genes = {r[4] for r in ann}
    print(f"  {len(ann):,} transcripts over {len(genes):,} genes -> {ann_p}")
    print()

    if args.check_against:
        print("-- overlap with the paired object's RNA var_names ------------------")
        try:
            import h5py
            with h5py.File(args.check_against, "r") as fh:
                idx = fh["mod/scRNA/var"].attrs.get("_index", "_index")
                if isinstance(idx, bytes):
                    idx = idx.decode()
                raw = fh[f"mod/scRNA/var/{idx}"][:]
            rna = {x.decode() if isinstance(x, bytes) else str(x) for x in raw}
        except Exception as exc:                          # pragma: no cover
            print(f"  could not read scRNA var_names: {type(exc).__name__}: {exc}")
            print("  (skipping the overlap check; it is advisory)")
        else:
            hit = rna & genes
            frac = len(hit) / max(len(rna), 1)
            print(f"  RNA var_names        : {len(rna):,}")
            print(f"  annotation genes     : {len(genes):,}")
            print(f"  matched              : {len(hit):,}  ({frac:.1%} of RNA)")
            if frac < 0.5:
                print()
                print("  WARNING: fewer than half the RNA genes are in the")
                print("  annotation. get_search_space intersects on the Gene")
                print("  column WITHOUT erroring on an empty result, so this")
                print("  becomes a silently tiny search space and few or no")
                print("  region-to-gene links. Check whether the RNA object")
                print("  uses Ensembl IDs (ENSG...) rather than HGNC symbols.")
            else:
                miss = sorted(rna - genes)[:8]
                print(f"  unmatched sample     : {miss}")
                print("  (non-coding genes and deprecated symbols are expected")
                print("   here -- the annotation is protein-coding only.)")
        print()

    print("=" * 70)
    print("Now confirm the DAG no longer plans download_genome_annotations:")
    print("    bash 03_pipeline/run_pipeline.sh "
          "--config 03_pipeline/config.yaml --cores 4 -n")
    print("It should drop from 13 jobs to 12.")
    print("=" * 70)


if __name__ == "__main__":
    main()
