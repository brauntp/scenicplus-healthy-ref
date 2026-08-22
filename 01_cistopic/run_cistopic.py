#!/usr/bin/env python3
"""
01_cistopic/run_cistopic.py -- canonical pycisTopic path for an unpaired
ArchR-derived hematopoietic reference.

Consumes the output of 01_cistopic/export_from_archr.R and runs, in order:

  1. create_cistopic_object()     count matrix + region/cell names + blacklist
  2. add_cell_data()              ArchR cellColData
  3. run_cgs_models_mallet()      LDA via MALLET, one model per --n-topics value
  4. evaluate_models()            model selection (or forced with --select-model)
  5. add_LDA_model()              attach the winner to the object
  6. compute_topic_metrics()      per-topic QC table
  7. impute_accessibility()       topic-imputed accessibility -> .h5ad for GLUE
  8. binarize_topics()            per-topic region sets -> BED
  9. find_highly_variable_features() + find_diff_features()
                                  per-cell-type DAR region sets -> BED

Every pycisTopic name and keyword below was checked against the cloned source at
./src/pycisTopic (commit 787ce42, pyproject version "2.0a"). Deviations from the
published tutorial are called out inline with "TUTORIAL MISMATCH".


=============================== THE TRAPS ===============================

TRAP 1 -- MALLET IS NOT CONFIGURED BY ARGUMENTS ALONE.
  There are two independent requirements and only one of them is a function
  keyword.
    (a) mallet_path=  is a real keyword on run_cgs_models_mallet(). It must
        point at the mallet *launcher script*, e.g.
        /path/Mallet-202108/bin/mallet -- not the directory, not the jar.
    (b) Heap size is NOT a keyword anywhere in pycisTopic. The bin/mallet shell
        script reads the environment variable MALLET_MEMORY and passes it to
        java as -Xmx. pycisTopic's own CLI does exactly
        `os.environ["MALLET_MEMORY"] = f"{gb}G"` (cli/subcommand/topic_modeling.py
        line 121) before calling run_cgs_models_mallet, and nothing inside
        lda_models.py ever touches it. So this script sets os.environ itself,
        before the import that matters, from --mallet-memory. Forget it and
        MALLET runs on the JVM default heap (~256 MB-1 GB) and dies with a Java
        OutOfMemoryError surfaced as an opaque
        `RuntimeError: command '[...]' return with error (code 1)`.
    Also: MALLET is Java. `module load java`/a JDK on PATH is a third
    requirement that neither pycisTopic nor this script can satisfy for you; we
    check for it and refuse early.

TRAP 2 -- MALLET TEMP FILES ARE HUGE AND DEFAULT TO /tmp.
  LDAMallet.__init__ does `self.tmp_dir = tmp_dir if tmp_dir else
  tempfile.gettempdir()`. The corpus text file it writes is one line per cell
  listing every accessible region id; for 100k cells x 300k peaks that is tens
  of GB. On a shared SLURM node /tmp is small and shared, and filling it takes
  down other people's jobs. Always pass --tmp-path to node-local scratch
  ($SLURM_TMPDIR / $TMPDIR). This script requires it to exist and be writable.
  Note the corpus lands at <tmp_path>/corpus.mallet with a FIXED name (no
  per-model randomisation -- only state/doctopics/inferencer/topickeys get a
  random label), so two concurrent runs sharing one --tmp-path will corrupt
  each other. Use a per-job directory.

TRAP 3 -- MODEL SCAN COST IS NOT WHAT run_cgs_models_mallet's n_cpu SUGGESTS.
  Despite the docstring for the non-MALLET sibling ("parallelization is done
  per model"), run_cgs_models_mallet builds its model list with a plain list
  comprehension -- the models run STRICTLY SEQUENTIALLY in this process, and
  n_cpu is forwarded to MALLET's own --num-threads, i.e. it parallelises
  *within* each model. So scanning 11 topic counts costs 11 sequential MALLET
  runs. Budget wall time accordingly, use --save-path so completed models are
  pickled as they finish (Topic<N>.pkl), and consider splitting the scan across
  SLURM array tasks with one --n-topics value each, then a final pass with
  --models-from to select.

TRAP 4 -- create_cistopic_object() RENAMES YOUR CELLS BY DEFAULT.
  tag_cells=True (the default) rewrites every cell name to
  f"{cell}{split_pattern}{project}". Our barcodes are already ArchR's
  "Sample#Barcode" and our metadata is keyed on exactly that, so this script
  passes tag_cells=False and split_pattern="___" (the pycisTopic default, kept
  so add_cell_data's prepare_tag_cells() fallback never fires -- with
  split_pattern="-" that helper applies barcode-shaped regexes that would
  mangle "Sample#Barcode").

TRAP 5 -- create_cistopic_object() SILENTLY DROPS REGIONS AND ITS OWN
  path_to_blacklist BRANCH IS ORDER-SENSITIVE. It (i) removes blacklist
  overlaps via pyranges, (ii) drops regions with zero accessible cells, and
  (iii) rebuilds region_names from pyranges output -- which is re-sorted
  relative to your input. That is fine, but it means cistopic_obj.region_names
  is NOT your regions.txt, and any file you wrote earlier keyed on input order
  is no longer aligned. Always read names back off the object. This script does,
  and logs how many regions were lost at each step.

TRAP 6 -- impute_accessibility() RETURNS AN OBJECT, NOT AN AnnData, AND ITS
  MATRIX IS REGIONS x CELLS. It returns CistopicImputedFeatures with
  .mtx (features x cells), .feature_names, .cell_names. pycisTopic has no h5ad
  writer at all (grep: anndata appears only in label_transfer.py). AnnData is
  cells x vars, so we transpose exactly once, on write, and set
  var_names = feature_names so region names travel to the GLUE pairing step.
  Also: with the default scale_factor=10**6 the returned matrix is DENSE int32,
  regions x cells. For 250k regions x 100k cells that is 100 GB. This script
  refuses to proceed past a --max-impute-gb guard and offers
  --impute-on-variable-regions to impute only the HVF subset, plus
  --impute-scale-factor 1 to get float32 instead.

TRAP 7 -- find_diff_features() USES RAY, AND RAY ON SLURM NEEDS A TEMP DIR.
  find_diff_features(..., n_cpu>1) calls ray.init(num_cpus=n_cpu, **kwargs) and
  forwards **kwargs straight to ray.init. The tutorial passes
  _temp_dir='/scratch/...' -- that is a ray.init argument, not a pycisTopic one,
  which is why it does not appear in the signature. Ray's default temp path is
  /tmp/ray and will exhaust a shared node. This script forwards
  _temp_dir=<--ray-temp-dir or --tmp-path>. With --dar-n-cpu 1 ray is never
  initialised at all, which on a contended cluster is often the safer choice.

TRAP 8 -- binarize_topics() MUTATES THE MODEL AND ITS 'ntop' THRESHOLD IS AN
  INDEX. It writes selected_model.topic_ass["Regions_in_binarized_topic"] (or
  "Cells_in_..."), so call it before pickling if you want those counts stored.
  With method="ntop" the threshold is `data.iloc[ntop]` on the sorted
  normalised distribution, so you get *approximately* ntop regions, and asking
  for ntop >= n_regions raises IndexError. Keys of the returned dict are
  "Topic1".."TopicN", 1-based.

TRAP 9 -- matplotlib MUST BE HEADLESS AND PLOTTING MUST BE OFF.
  evaluate_models(plot=True) and binarize_topics(plot=True) call plt.show().
  Under SLURM with no display that is at best a warning and at worst a hang.
  This script forces the Agg backend before importing pycisTopic and passes
  plot=False everywhere, using the functions' own save= arguments to get PDFs.

TRAP 10 -- REGION-SET BEDs FOR SCENIC+ MUST BE 3-COLUMN AND SORTED.
  region_names_to_coordinates() returns a DataFrame indexed by region name with
  Chromosome/Start/End columns. pycistarget wants headerless, index-less,
  coordinate-sorted BED. We sort by (Chromosome, Start, End) and write
  header=False, index=False. An empty region set is written as a 0-byte file and
  logged loudly rather than skipped -- a missing file downstream looks like a
  pipeline bug, an empty one looks like what it is.


=========================== TYPICAL INVOCATION ==========================

  export MALLET_HOME=/path/Mallet-202108        # optional, for your own sanity
  module load java

  python 01_cistopic/run_cistopic.py \
    --matrix         archr_export/peak_matrix.mtx \
    --barcodes       archr_export/barcodes.tsv \
    --regions        archr_export/regions.txt \
    --cell-metadata  archr_export/cell_metadata.tsv \
    --blacklist      resources/hg38-blacklist.v2.bed \
    --out-dir        01_cistopic/out \
    --group-col      CellType \
    --n-topics       2 5 10 15 20 25 30 35 40 45 50 \
    --n-iter         500 \
    --mallet-path    /path/Mallet-202108/bin/mallet \
    --mallet-memory  200G \
    --tmp-path       "$SLURM_TMPDIR/mallet" \
    --n-cpu          "$SLURM_CPUS_PER_TASK"

  hg38-blacklist.v2.bed ships inside the pycisTopic repo at blacklist/.

Exit 0 on success. Any missing input, unwritable path, absent MALLET/java, or
absent --group-col is a hard exit 1 with a "FATAL:" line on stderr. There are no
silent fallbacks.
"""

from __future__ import annotations

import argparse
import logging
import os
import pickle
import shutil
import subprocess
import sys
import time
from pathlib import Path

# --- headless before anything imports pyplot (TRAP 9) -----------------------
import matplotlib

matplotlib.use("Agg")

LOG_FMT = "[%(asctime)s] %(levelname)-7s %(message)s"
logging.basicConfig(
    level=logging.INFO, format=LOG_FMT, datefmt="%Y-%m-%d %H:%M:%S",
    stream=sys.stderr,
)
log = logging.getLogger("run_cistopic")


def die(msg: str, code: int = 1) -> None:
    """Hard stop with a greppable prefix. Never returns."""
    print(f"FATAL: {msg}", file=sys.stderr, flush=True)
    sys.exit(code)


def need_file(path: str | None, what: str) -> Path:
    if path is None:
        die(f"{what} not given")
    p = Path(path).expanduser()
    if not p.exists():
        die(f"{what} does not exist: {p}")
    if not p.is_file():
        die(f"{what} is not a file: {p}")
    if p.stat().st_size == 0:
        die(f"{what} is empty: {p}")
    return p


def need_dir(path: str, what: str, create: bool = False) -> Path:
    p = Path(path).expanduser()
    if create:
        p.mkdir(parents=True, exist_ok=True)
    if not p.is_dir():
        die(f"{what} is not an existing directory: {p}")
    if not os.access(p, os.W_OK):
        die(f"{what} is not writable: {p}")
    return p


# ===========================================================================
# arguments
# ===========================================================================
def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        prog="run_cistopic.py",
        description="LDA topic modelling on an ArchR-exported peak matrix (pycisTopic 2.0a).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    g = ap.add_argument_group("inputs (from export_from_archr.R)")
    g.add_argument("--matrix", required=True,
                   help="Matrix Market .mtx, regions(rows) x cells(cols).")
    g.add_argument("--barcodes", required=True,
                   help="one cell name per line, in .mtx column order.")
    g.add_argument("--regions", required=True,
                   help="one region name 'chr:start-end' per line, in .mtx row order. "
                        "peaks.bed also accepted (column 4 is used).")
    g.add_argument("--cell-metadata", required=True,
                   help="TSV of cellColData; first column is the cell id.")
    g.add_argument("--blacklist", required=True,
                   help="ENCODE blacklist BED, e.g. hg38-blacklist.v2.bed. Passed to "
                        "create_cistopic_object(path_to_blacklist=...). Required "
                        "explicitly so it is never skipped by accident.")

    g = ap.add_argument_group("output")
    g.add_argument("--out-dir", required=True)
    g.add_argument("--region-set-folder", default=None,
                   help="where per-topic and per-cell-type BEDs go. "
                        "Default: <out-dir>/region_sets")
    g.add_argument("--project", default="cisTopic_hema",
                   help="cisTopic project name (stored on the object).")

    g = ap.add_argument_group("object construction")
    g.add_argument("--group-col", default=None,
                   help="cell_metadata column with the cell-type label. Required for "
                        "per-cell-type DARs; omit to skip that step.")
    g.add_argument("--min-frag", type=int, default=1,
                   help="create_cistopic_object(min_frag=). Cells were already QC'd in "
                        "ArchR; leave at 1 unless you mean it.")
    g.add_argument("--min-cell", type=int, default=1,
                   help="create_cistopic_object(min_cell=).")
    g.add_argument("--is-acc", type=int, default=1,
                   help="create_cistopic_object(is_acc=): fragments needed for a region "
                        "to count as accessible when binarising.")
    g.add_argument("--split-pattern", default="___",
                   help="pycisTopic barcode/sample separator. Keep the default: our "
                        "cell names are ArchR 'Sample#Barcode' and '-' triggers "
                        "prepare_tag_cells' barcode regexes (TRAP 4).")

    g = ap.add_argument_group("LDA / MALLET")
    g.add_argument("--n-topics", type=int, nargs="+",
                   default=[2, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50],
                   help="topic counts to scan. Each is a separate sequential MALLET run "
                        "(TRAP 3).")
    g.add_argument("--n-iter", type=int, default=500,
                   help="Gibbs iterations (pycisTopic default is 150; the tutorial uses 500).")
    g.add_argument("--random-state", type=int, default=555)
    g.add_argument("--alpha", type=float, default=50.0)
    g.add_argument("--alpha-by-topic", dest="alpha_by_topic",
                   action="store_true", default=True)
    g.add_argument("--no-alpha-by-topic", dest="alpha_by_topic", action="store_false")
    g.add_argument("--eta", type=float, default=0.1)
    g.add_argument("--eta-by-topic", dest="eta_by_topic",
                   action="store_true", default=False)
    g.add_argument("--top-topics-coh", type=int, default=5,
                   help="topics averaged for the coherence metric.")
    g.add_argument("--n-cpu", type=int, default=1,
                   help="forwarded to run_cgs_models_mallet(n_cpu=) -> MALLET "
                        "--num-threads, i.e. threads WITHIN one model (TRAP 3).")
    g.add_argument("--mallet-path", default="mallet",
                   help="path to the mallet launcher script, e.g. "
                        "/path/Mallet-202108/bin/mallet.")
    g.add_argument("--mallet-memory", default="200G",
                   help="JVM heap for MALLET. Exported as MALLET_MEMORY; NOT a "
                        "pycisTopic keyword (TRAP 1). Accepts '200G' or '200'.")
    g.add_argument("--tmp-path", required=True,
                   help="node-local scratch for the MALLET corpus. Tens of GB "
                        "(TRAP 2). Must be per-job.")
    g.add_argument("--save-path", default=None,
                   help="directory where each finished model is pickled as "
                        "Topic<N>.pkl. Default: <out-dir>/models. Strongly "
                        "recommended for long scans.")
    g.add_argument("--reuse-corpus", action="store_true",
                   help="reuse <tmp-path>/corpus.mallet if present. Only safe if the "
                        "cisTopic object is byte-identical to the previous run.")
    g.add_argument("--select-model", type=int, default=None,
                   help="force the model with this many topics instead of letting "
                        "evaluate_models pick.")
    g.add_argument("--models-from", default=None,
                   help="skip LDA: load a pickled list of models (or a directory of "
                        "Topic<N>.pkl) and go straight to selection.")

    g = ap.add_argument_group("imputation / h5ad")
    g.add_argument("--impute-scale-factor", type=int, default=10 ** 6,
                   help="impute_accessibility(scale_factor=). An int != 1 gives an "
                        "int32 matrix; pass 1 for float32 (TRAP 6).")
    g.add_argument("--impute-chunk-size", type=int, default=20000,
                   help="impute_accessibility(chunk_size=), in regions.")
    g.add_argument("--impute-on-variable-regions", action="store_true",
                   help="impute only highly variable regions. Cuts the dense matrix "
                        "down hard; costs a first pass over all regions.")
    g.add_argument("--max-impute-gb", type=float, default=64.0,
                   help="refuse to build a dense imputed matrix larger than this.")
    g.add_argument("--h5ad-out", default=None,
                   help="Default: <out-dir>/topic_imputed_accessibility.h5ad")

    g = ap.add_argument_group("binarization / DARs")
    g.add_argument("--binarize-methods", nargs="+", default=["otsu", "ntop"],
                   choices=["otsu", "yen", "li", "aucell", "ntop"],
                   help="region-topic binarization methods; one BED folder each.")
    g.add_argument("--ntop", type=int, default=3000,
                   help="regions per topic for method 'ntop'.")
    g.add_argument("--binarize-cells", action="store_true",
                   help="also binarize the cell-topic distribution (method 'li'), "
                        "written as TSV not BED.")
    g.add_argument("--skip-dars", action="store_true")
    g.add_argument("--dar-n-cpu", type=int, default=1,
                   help="find_diff_features(n_cpu=). >1 initialises ray (TRAP 7).")
    g.add_argument("--dar-adjpval", type=float, default=0.05)
    g.add_argument("--dar-log2fc", type=float, default=None,
                   help="Default: log2(1.5), matching pycisTopic's own default.")
    g.add_argument("--hvf-min-disp", type=float, default=0.05)
    g.add_argument("--hvf-min-mean", type=float, default=0.0125)
    g.add_argument("--hvf-max-mean", type=float, default=3.0)
    g.add_argument("--hvf-n-top", type=int, default=None,
                   help="find_highly_variable_features(n_top_features=). If set, the "
                        "dispersion/mean thresholds are ignored by pycisTopic.")
    g.add_argument("--ray-temp-dir", default=None,
                   help="ray _temp_dir for find_diff_features. Default: <tmp-path>/ray")

    ap.add_argument("--skip-mallet-check", action="store_true",
                    help="do not probe the mallet binary/java before running. Only for "
                         "environments where the probe itself is the problem.")
    return ap


# ===========================================================================
# MALLET environment (TRAP 1)
# ===========================================================================
def configure_mallet(mallet_path: str, mallet_memory: str, tmp_path: Path,
                     skip_check: bool) -> str:
    """Set MALLET_MEMORY and verify the binary + a JVM exist.

    MALLET_MEMORY is the ONLY way to give MALLET heap: bin/mallet reads it and
    passes it to java as -Xmx. pycisTopic never sets or reads it (only its CLI
    does, in cli/subcommand/topic_modeling.py). Nothing in
    run_cgs_models_mallet's signature controls memory.
    """
    mem = mallet_memory.strip()
    if not mem:
        die("--mallet-memory is empty")
    if mem[-1].isdigit():
        mem = f"{mem}G"
        log.info("--mallet-memory had no unit; interpreting as %s", mem)
    if mem[-1].upper() not in ("G", "M"):
        die(f"--mallet-memory must end in G or M (or be a bare number of GB), got {mallet_memory!r}")
    os.environ["MALLET_MEMORY"] = mem
    log.info("MALLET_MEMORY=%s exported (bin/mallet -> java -Xmx%s)", mem, mem)

    resolved = mallet_path
    if os.sep in mallet_path or mallet_path.startswith("."):
        p = Path(mallet_path).expanduser()
        if not p.exists():
            die(f"--mallet-path does not exist: {p}\n"
                f"  Download a release from https://github.com/mimno/Mallet/releases and "
                f"point at <dir>/bin/mallet.")
        if p.is_dir():
            die(f"--mallet-path is a directory: {p}\n"
                f"  Point at the launcher script itself, probably {p / 'bin' / 'mallet'}")
        if not os.access(p, os.X_OK):
            die(f"--mallet-path is not executable: {p}  (chmod +x it)")
        resolved = str(p.resolve())
    else:
        found = shutil.which(mallet_path)
        if found is None:
            die(f"--mallet-path {mallet_path!r} is not on PATH and is not a path. "
                f"Give the full path to bin/mallet.")
        resolved = found
    log.info("mallet binary: %s", resolved)

    if not skip_check:
        if shutil.which("java") is None:
            die("no 'java' on PATH. MALLET is a Java program; pycisTopic cannot "
                "supply a JVM. On ARC: module load java (or openjdk), then re-run. "
                "Use --skip-mallet-check only if you know java is resolved "
                "some other way inside bin/mallet.")
        # bin/mallet with no subcommand prints usage and exits non-zero; we only
        # care that it executes at all and that the JVM starts.
        try:
            proc = subprocess.run([resolved], capture_output=True, text=True, timeout=180)
        except (OSError, subprocess.TimeoutExpired) as exc:
            die(f"could not execute the mallet binary ({exc}). Check the shebang and "
                f"that MALLET's lib/*.jar are present next to bin/.")
        blob = (proc.stdout or "") + (proc.stderr or "")
        # Only patterns the JVM emits when it cannot reserve the heap. Do not
        # broaden to a bare "heap" match: MALLET's own usage text is printed here
        # and a substring hit would fail a working setup.
        heap_fail = (
            "OutOfMemoryError" in blob
            or "Unable to allocate" in blob
            or "Could not reserve enough space" in blob
            or "Invalid maximum heap size" in blob
            or "Too small maximum heap" in blob
        )
        if heap_fail:
            die(f"the JVM refused MALLET_MEMORY={mem} at startup. Lower "
                f"--mallet-memory below the SLURM allocation for this job.\n"
                f"  mallet said: {blob.strip()[:500]}")
        if "Error occurred during initialization" in blob or "ClassNotFound" in blob:
            die(f"mallet started but the JVM failed:\n  {blob.strip()[:800]}")
        log.info("mallet + java probe ok")

    need_dir(str(tmp_path), "--tmp-path", create=True)
    corpus = tmp_path / "corpus.mallet"
    if corpus.exists():
        log.warning("%s already exists. run_cgs_models_mallet uses a FIXED corpus "
                    "filename, so a concurrent job sharing this --tmp-path will "
                    "corrupt it (TRAP 2). %s", corpus,
                    "It will be reused (--reuse-corpus)." if
                    "--reuse-corpus" in sys.argv else "It will be overwritten.")
    free_gb = shutil.disk_usage(tmp_path).free / 1024 ** 3
    log.info("--tmp-path %s has %.1f GB free", tmp_path, free_gb)
    if free_gb < 20:
        log.warning("under 20 GB free on --tmp-path; the MALLET corpus for a "
                    "100k-cell object will not fit (TRAP 2)")
    return resolved


# ===========================================================================
# input loading
# ===========================================================================
def read_region_names(path: Path) -> list[str]:
    """regions.txt (one name per line) or peaks.bed (col 4, else built from 1-3)."""
    names: list[str] = []
    with path.open() as fh:
        for ln, raw in enumerate(fh, 1):
            line = raw.rstrip("\n")
            if not line or line.startswith(("#", "track", "browser")):
                continue
            parts = line.split("\t")
            if len(parts) == 1:
                names.append(parts[0].strip())
            elif len(parts) >= 4 and ":" in parts[3]:
                names.append(parts[3].strip())
            elif len(parts) >= 3:
                names.append(f"{parts[0]}:{parts[1]}-{parts[2]}")
            else:
                die(f"{path}:{ln}: cannot parse a region name from {line!r}")
    if not names:
        die(f"no region names read from {path}")
    bad = [n for n in names[:1000] if ":" not in n or "-" not in n.split(":", 1)[1]]
    if bad:
        die(f"region names must be 'chr:start-end' for pycisTopic "
            f"(region_names_to_coordinates splits on the first ':' then '-'). "
            f"Offenders: {bad[:3]}")
    return names


def write_region_set_bed(region_names, out_path: Path, label: str) -> int:
    """3-column sorted BED from region names (TRAP 10)."""
    from pycisTopic.utils import region_names_to_coordinates

    names = list(region_names)
    if not names:
        out_path.write_text("")
        log.warning("region set %r is EMPTY; wrote a 0-byte %s so the file exists",
                    label, out_path.name)
        return 0
    df = region_names_to_coordinates(names)
    df = df.sort_values(["Chromosome", "Start", "End"])
    df.to_csv(out_path, sep="\t", header=False, index=False)
    return len(df)


# ===========================================================================
# main
# ===========================================================================
def main() -> None:
    t0 = time.time()
    args = build_parser().parse_args()

    log.info("run_cistopic.py starting; python %s", sys.version.split()[0])

    f_mtx = need_file(args.matrix, "--matrix")
    f_bc = need_file(args.barcodes, "--barcodes")
    f_reg = need_file(args.regions, "--regions")
    f_meta = need_file(args.cell_metadata, "--cell-metadata")
    f_bl = need_file(args.blacklist, "--blacklist")

    out_dir = need_dir(args.out_dir, "--out-dir", create=True)
    rs_dir = need_dir(args.region_set_folder or str(out_dir / "region_sets"),
                      "--region-set-folder", create=True)
    save_path = need_dir(args.save_path or str(out_dir / "models"),
                         "--save-path", create=True)
    tmp_path = Path(args.tmp_path).expanduser()
    h5ad_out = Path(args.h5ad_out) if args.h5ad_out else out_dir / "topic_imputed_accessibility.h5ad"
    f_obj = out_dir / "cistopic_obj.pkl"

    if args.min_frag < 1 or args.min_cell < 1 or args.is_acc < 1:
        die("--min-frag/--min-cell/--is-acc must all be >= 1")
    if any(n < 2 for n in args.n_topics):
        die(f"--n-topics values must be >= 2, got {args.n_topics}")
    if len(set(args.n_topics)) != len(args.n_topics):
        die(f"--n-topics has duplicates: {args.n_topics}")
    if args.select_model is not None and args.select_model not in args.n_topics:
        die(f"--select-model {args.select_model} is not among --n-topics {args.n_topics}. "
            f"evaluate_models does `models[all_topics.index(select_model)]` and "
            f"would raise ValueError.")
    # Two ways an automatic evaluate_models selection blows up, both worth
    # catching here rather than after hours of MALLET:
    #   1. every scanned topic count below min_topics_coh (default 5) leaves the
    #      Minmo_2011 index list empty, and the rescale divides by an empty range;
    #   2. a single model gives max()==min() in every rescale -> nan metrics.
    if args.select_model is None:
        if max(args.n_topics) < 5:
            die(f"--n-topics {args.n_topics} are all below evaluate_models' "
                f"min_topics_coh default of 5, so the coherence metric has no "
                f"eligible models and selection fails. Scan at least one value "
                f">= 5, or pass --select-model to skip metric-based selection.")
        if len(args.n_topics) == 1:
            die(f"--n-topics has a single value ({args.n_topics[0]}) and no "
                f"--select-model. evaluate_models rescales each metric by "
                f"(max - min) across models, which is 0 here. Either scan "
                f"several topic counts or pass --select-model {args.n_topics[0]}.")

    mallet_bin = configure_mallet(args.mallet_path, args.mallet_memory,
                                  tmp_path, args.skip_mallet_check)

    # -- imports AFTER matplotlib.use and MALLET_MEMORY are settled ---------
    import numpy as np
    import pandas as pd
    import scipy.io
    import scipy.sparse as sp

    from pycisTopic.cistopic_class import create_cistopic_object
    from pycisTopic.lda_models import evaluate_models, run_cgs_models_mallet
    from pycisTopic.topic_binarization import binarize_topics
    from pycisTopic.topic_qc import compute_topic_metrics
    from pycisTopic.diff_features import (
        impute_accessibility,
        normalize_scores,
        find_highly_variable_features,
        find_diff_features,
    )

    try:
        import pycisTopic
        log.info("pycisTopic %s", getattr(pycisTopic, "__version__", "unknown"))
    except Exception:  # noqa: BLE001 - version reporting must never be fatal
        pass

    # =====================================================================
    # 1. matrix + names
    # =====================================================================
    log.info("reading %s", f_mtx)
    mtx = scipy.io.mmread(str(f_mtx))
    # create_cistopic_object indexes rows (fragment_matrix[index,]) and columns,
    # and hands the matrix to sklearn.preprocessing.binarize. COO supports
    # neither cleanly, so convert to CSR up front.
    frag_matrix = sp.csr_matrix(mtx, dtype=np.int32)
    del mtx

    cell_names = [l.strip() for l in f_bc.read_text().splitlines() if l.strip()]
    region_names = read_region_names(f_reg)
    log.info("matrix %d regions x %d cells, %s nonzeros",
             frag_matrix.shape[0], frag_matrix.shape[1],
             f"{frag_matrix.nnz:,}")

    if frag_matrix.shape[0] != len(region_names):
        die(f"--matrix has {frag_matrix.shape[0]} rows but --regions has "
            f"{len(region_names)} names. The exporter writes regions as ROWS; "
            f"if this file is cells x regions, transpose it -- do not swap the "
            f"name lists.")
    if frag_matrix.shape[1] != len(cell_names):
        die(f"--matrix has {frag_matrix.shape[1]} columns but --barcodes has "
            f"{len(cell_names)} names.")
    if len(set(cell_names)) != len(cell_names):
        die("duplicate cell names in --barcodes")
    if len(set(region_names)) != len(region_names):
        die("duplicate region names in --regions")

    # =====================================================================
    # 2. metadata
    # =====================================================================
    log.info("reading %s", f_meta)
    cell_data = pd.read_table(f_meta, index_col=0, low_memory=False)
    log.info("metadata: %d rows x %d columns", *cell_data.shape)

    if args.group_col is not None and args.group_col not in cell_data.columns:
        print(f"FATAL: --group-col {args.group_col!r} not in {f_meta}", file=sys.stderr)
        print(f"Available columns ({len(cell_data.columns)}):", file=sys.stderr)
        for c in cell_data.columns:
            print(f"  - {c}", file=sys.stderr)
        sys.exit(1)

    overlap = len(set(cell_names) & set(cell_data.index))
    log.info("%d/%d matrix cells found in the metadata index", overlap, len(cell_names))
    if overlap == 0:
        die(f"no overlap between --barcodes and the --cell-metadata index.\n"
            f"  barcode example : {cell_names[0]!r}\n"
            f"  metadata example: {str(cell_data.index[0])!r}\n"
            f"  add_cell_data would fill every column with NaN.")
    if overlap < len(cell_names):
        log.warning("%d matrix cells have no metadata row; add_cell_data fills "
                    "them with NaN and find_diff_features drops them",
                    len(cell_names) - overlap)

    # =====================================================================
    # 3. CistopicObject
    # =====================================================================
    # VERIFIED create_cistopic_object(fragment_matrix, cell_names=None,
    #   region_names=None, path_to_blacklist=None, min_frag=1, min_cell=1,
    #   is_acc=1, path_to_fragments=None, project='cisTopic', tag_cells=True,
    #   split_pattern='___')
    # tag_cells=False keeps ArchR's "Sample#Barcode" intact (TRAP 4).
    log.info("create_cistopic_object(); blacklist=%s", f_bl)
    n_regions_in = len(region_names)
    cistopic_obj = create_cistopic_object(
        fragment_matrix=frag_matrix,
        cell_names=cell_names,
        region_names=region_names,
        path_to_blacklist=str(f_bl),
        min_frag=args.min_frag,
        min_cell=args.min_cell,
        is_acc=args.is_acc,
        project=args.project,
        tag_cells=False,
        split_pattern=args.split_pattern,
    )
    log.info("%s", cistopic_obj)
    log.info("regions: %d in -> %d kept (%d dropped by blacklist overlap + "
             "zero-accessibility filtering; region order is now pyranges' order, "
             "NOT regions.txt order -- TRAP 5)",
             n_regions_in, len(cistopic_obj.region_names),
             n_regions_in - len(cistopic_obj.region_names))
    if len(cistopic_obj.region_names) == 0:
        die("all regions were filtered out. Check that --blacklist and the peak "
            "set use the same chromosome naming ('chr1' vs '1').")
    if len(cistopic_obj.cell_names) == 0:
        die("all cells were filtered out (--min-frag too high?)")
    del frag_matrix

    # VERIFIED CistopicObject.add_cell_data(cell_data, split_pattern='___')
    # In-place; returns None. Overlapping columns are overwritten with a printed
    # notice. Index must be cell names.
    log.info("add_cell_data()")
    cistopic_obj.add_cell_data(cell_data, split_pattern=args.split_pattern)
    if args.group_col is not None:
        if args.group_col not in cistopic_obj.cell_data.columns:
            die(f"{args.group_col!r} did not survive add_cell_data; "
                f"cell_data columns are {list(cistopic_obj.cell_data.columns)}")
        n_lab = int(cistopic_obj.cell_data[args.group_col].notna().sum())
        log.info("group column %r attached: %d/%d cells labelled, %d levels",
                 args.group_col, n_lab, len(cistopic_obj.cell_names),
                 cistopic_obj.cell_data[args.group_col].nunique(dropna=True))
        if n_lab == 0:
            die(f"{args.group_col!r} is all-NaN on the object after add_cell_data; "
                f"the metadata index does not match the object's cell names")

    with open(out_dir / "cistopic_obj_premodel.pkl", "wb") as fh:
        pickle.dump(cistopic_obj, fh)
    log.info("wrote cistopic_obj_premodel.pkl (restart point: pass it to "
             "--models-from workflows or reuse it for a wider --n-topics scan)")

    # =====================================================================
    # 4. LDA
    # =====================================================================
    if args.models_from:
        src = Path(args.models_from).expanduser()
        if src.is_dir():
            pkls = sorted(src.glob("Topic*.pkl"))
            if not pkls:
                die(f"--models-from {src} contains no Topic*.pkl")
            models = []
            for p in pkls:
                with p.open("rb") as fh:
                    models.append(pickle.load(fh))
            log.info("loaded %d models from %s", len(models), src)
        else:
            with need_file(str(src), "--models-from").open("rb") as fh:
                models = pickle.load(fh)
            if not isinstance(models, list):
                models = [models]
            log.info("loaded %d models from %s", len(models), src)
    else:
        # VERIFIED run_cgs_models_mallet(cistopic_obj, n_topics, n_cpu=1,
        #   n_iter=150, random_state=555, alpha=50.0, alpha_by_topic=True,
        #   eta=0.1, eta_by_topic=False, top_topics_coh=5, tmp_path=None,
        #   save_path=None, reuse_corpus=False, mallet_path='mallet')
        # -> list[CistopicLDAModel]. NOTE: there is no memory keyword; heap comes
        # from MALLET_MEMORY (TRAP 1). Models run sequentially (TRAP 3).
        log.info("run_cgs_models_mallet(): %d models %s, %d iterations, "
                 "n_cpu=%d (MALLET --num-threads)",
                 len(args.n_topics), args.n_topics, args.n_iter, args.n_cpu)
        log.info("models are trained SEQUENTIALLY; each finished model is "
                 "pickled to %s/Topic<N>.pkl", save_path)
        t_lda = time.time()
        try:
            models = run_cgs_models_mallet(
                cistopic_obj,
                n_topics=list(args.n_topics),
                n_cpu=args.n_cpu,
                n_iter=args.n_iter,
                random_state=args.random_state,
                alpha=args.alpha,
                alpha_by_topic=args.alpha_by_topic,
                eta=args.eta,
                eta_by_topic=args.eta_by_topic,
                top_topics_coh=args.top_topics_coh,
                tmp_path=str(tmp_path),
                save_path=str(save_path),
                reuse_corpus=args.reuse_corpus,
                mallet_path=mallet_bin,
            )
        except RuntimeError as exc:
            die(f"MALLET failed. pycisTopic wraps the subprocess error, so the "
                f"Java message is inside this text:\n  {exc}\n"
                f"  Most common causes: (1) MALLET_MEMORY={os.environ.get('MALLET_MEMORY')} "
                f"exceeds the SLURM memory allocation, so the JVM cannot reserve "
                f"the heap; (2) --tmp-path {tmp_path} filled up; (3) java not on "
                f"PATH inside the job. See TRAP 1/TRAP 2 in this file's docstring.")
        log.info("LDA finished in %.1f min", (time.time() - t_lda) / 60)

        with open(out_dir / "models.pkl", "wb") as fh:
            pickle.dump(models, fh)

    if not models:
        die("no models to select from")

    # =====================================================================
    # 5. model selection
    # =====================================================================
    # VERIFIED evaluate_models(models, select_model=None, return_model=True,
    #   metrics=['Minmo_2011','loglikelihood','Cao_Juan_2009','Arun_2010'],
    #   min_topics_coh=5, plot=True, figsize=(6.4,4.8), plot_metrics=False,
    #   save=None)
    # NOTE the misspelling 'Minmo_2011' in the metrics default -- the coherence
    # metric key really is spelled that way in the parameter list (the per-topic
    # coherence column is 'Mimno_2011'). Do not "fix" it when passing metrics=.
    # plot=False to avoid plt.show() under SLURM (TRAP 9); save= still writes a PDF.
    log.info("evaluate_models(select_model=%s)", args.select_model)
    model = evaluate_models(
        models,
        select_model=args.select_model,
        return_model=True,
        plot=False,
        save=str(out_dir / "model_selection.pdf"),
    )
    if model is None:
        die("evaluate_models returned None (return_model=True should prevent this)")
    log.info("selected model: %s", model)

    # VERIFIED CistopicObject.add_LDA_model(model) -- in place. It also sets
    # model.region_topic = model.topic_region.loc[self.region_names], which will
    # KeyError if the model was trained on a different object than this one.
    try:
        cistopic_obj.add_LDA_model(model)
    except KeyError as exc:
        die(f"add_LDA_model could not align the model's regions to this object "
            f"({exc}). The models were trained on a different cisTopic object; "
            f"re-run LDA or pass the matching --models-from.")

    # VERIFIED compute_topic_metrics(cistopic_obj, return_metrics=True)
    # -> pd.DataFrame, also stored at selected_model.topic_qc_metrics.
    topic_qc = compute_topic_metrics(cistopic_obj, return_metrics=True)
    topic_qc.to_csv(out_dir / "topic_qc_metrics.tsv", sep="\t")
    log.info("wrote topic_qc_metrics.tsv (%d topics)", topic_qc.shape[0])

    # =====================================================================
    # 6. binarization (before pickling, so topic_ass counts are stored -- TRAP 8)
    # =====================================================================
    # VERIFIED binarize_topics(cistopic_obj, target='region', method='otsu',
    #   smooth_topics=True, ntop=2000, predefined_thr=None, nbins=100,
    #   plot=False, figsize=(6.4,4.8), num_columns=1, save=None)
    # -> dict {'Topic1': DataFrame indexed by region name}
    n_regions_obj = len(cistopic_obj.region_names)
    binarized: dict[str, dict] = {}
    for method in args.binarize_methods:
        kw = {}
        if method == "ntop":
            if args.ntop >= n_regions_obj:
                die(f"--ntop {args.ntop} >= {n_regions_obj} regions on the object; "
                    f"binarize_topics indexes data.iloc[ntop] and will raise "
                    f"IndexError")
            kw["ntop"] = args.ntop
        label = f"Topics_{method}" + (f"_{args.ntop}" if method == "ntop" else "")
        log.info("binarize_topics(target='region', method=%r)%s", method,
                 f" ntop={args.ntop}" if method == "ntop" else "")
        binarized[label] = binarize_topics(
            cistopic_obj, target="region", method=method,
            plot=False, num_columns=5,
            save=str(out_dir / f"binarization_{label}.pdf"),
            **kw,
        )
        folder = rs_dir / label
        folder.mkdir(parents=True, exist_ok=True)
        total = 0
        for topic, df in binarized[label].items():
            total += write_region_set_bed(df.index, folder / f"{topic}.bed",
                                          f"{label}/{topic}")
        log.info("  %s: %d topics, %d regions total -> %s",
                 label, len(binarized[label]), total, folder)

    if args.binarize_cells:
        log.info("binarize_topics(target='cell', method='li')")
        cell_bin = binarize_topics(cistopic_obj, target="cell", method="li",
                                   plot=False, num_columns=5, nbins=100,
                                   save=str(out_dir / "binarization_cells_li.pdf"))
        cdir = out_dir / "binarized_cell_topic"
        cdir.mkdir(exist_ok=True)
        for topic, df in cell_bin.items():
            df.to_csv(cdir / f"{topic}.tsv", sep="\t")
        log.info("  wrote %d binarized cell-topic tables to %s", len(cell_bin), cdir)

    with open(f_obj, "wb") as fh:
        pickle.dump(cistopic_obj, fh)
    log.info("wrote %s (model + topic QC + binarization counts attached)", f_obj)

    # =====================================================================
    # 7. imputed accessibility -> h5ad for the GLUE pairing step
    # =====================================================================
    n_cells_obj = len(cistopic_obj.cell_names)
    itemsize = 4  # int32 when scale_factor is an int != 1, else float32
    dense_gb = n_regions_obj * n_cells_obj * itemsize / 1024 ** 3
    log.info("imputed matrix would be %d regions x %d cells dense = %.1f GB",
             n_regions_obj, n_cells_obj, dense_gb)

    selected_regions = None
    if args.impute_on_variable_regions:
        # First pass on all regions purely to rank them, then impute the subset.
        log.info("--impute-on-variable-regions: first pass to rank regions")
        if dense_gb > args.max_impute_gb:
            die(f"even the ranking pass needs {dense_gb:.1f} GB > "
                f"--max-impute-gb {args.max_impute_gb}. Raise the guard on a "
                f"high-memory partition, or drop this flag and subset regions "
                f"yourself before this script.")
        probe = impute_accessibility(
            cistopic_obj, selected_cells=None, selected_regions=None,
            scale_factor=args.impute_scale_factor,
            chunk_size=args.impute_chunk_size,
        )
        probe_norm = normalize_scores(probe, scale_factor=10 ** 4)
        selected_regions = find_highly_variable_features(
            probe_norm,
            min_disp=args.hvf_min_disp, min_mean=args.hvf_min_mean,
            max_mean=args.hvf_max_mean, max_disp=np.inf,
            n_bins=20, n_top_features=args.hvf_n_top, plot=False,
            save=str(out_dir / "hvf_dispersion_probe.pdf"),
        )
        log.info("  %d variable regions retained", len(selected_regions))
        del probe, probe_norm
        dense_gb = len(selected_regions) * n_cells_obj * itemsize / 1024 ** 3

    if dense_gb > args.max_impute_gb:
        die(f"imputed accessibility would need {dense_gb:.1f} GB of dense memory, "
            f"over --max-impute-gb {args.max_impute_gb}. impute_accessibility "
            f"allocates np.empty((n_regions, n_cells)) up front -- there is no "
            f"sparse path (TRAP 6). Options: --impute-on-variable-regions, "
            f"a subset of cells, or a high-memory partition with a raised guard.")

    # VERIFIED impute_accessibility(cistopic_obj, selected_cells=None,
    #   selected_regions=None, scale_factor=10**6, chunk_size=20000,
    #   project='cisTopic_Impute') -> CistopicImputedFeatures
    #   (.mtx = features x cells dense ndarray, .feature_names, .cell_names)
    # Regions whose whole imputed row is zero are DROPPED, so feature_names is
    # not necessarily the object's region_names.
    log.info("impute_accessibility(scale_factor=%d, chunk_size=%d)",
             args.impute_scale_factor, args.impute_chunk_size)
    imputed = impute_accessibility(
        cistopic_obj,
        selected_cells=None,
        selected_regions=selected_regions,
        scale_factor=args.impute_scale_factor,
        chunk_size=args.impute_chunk_size,
    )
    log.info("%s", imputed)

    # -- write h5ad -------------------------------------------------------
    # pycisTopic has no h5ad writer (anndata appears only in label_transfer.py),
    # so this is ours. AnnData is cells x vars; imputed.mtx is regions x cells.
    try:
        import anndata as ad
    except ImportError:
        die("anndata is not installed, so the .h5ad for the GLUE pairing step "
            "cannot be written. pip install anndata. The cisTopic object and "
            "region sets above are already saved, so re-running with "
            "--models-from will skip LDA.")

    log.info("writing %s (cells x regions, var_names = region names)", h5ad_out)
    X = imputed.mtx
    X = X.T.tocsr() if sp.issparse(X) else np.asarray(X).T
    adata = ad.AnnData(
        X=X,
        obs=pd.DataFrame(index=pd.Index(list(imputed.cell_names), name="cell_id")),
        var=pd.DataFrame(index=pd.Index(list(imputed.feature_names), name="region")),
    )
    if adata.shape != (len(imputed.cell_names), len(imputed.feature_names)):
        die(f"AnnData shape {adata.shape} does not match "
            f"({len(imputed.cell_names)}, {len(imputed.feature_names)}) after "
            f"transpose; refusing to write a mislabelled h5ad")

    # Carry the metadata across so the GLUE step does not have to re-join.
    meta_aligned = cistopic_obj.cell_data.reindex(adata.obs_names)
    for col in meta_aligned.columns:
        s = meta_aligned[col]
        if s.dtype == object:
            adata.obs[col] = s.astype(str).values
        else:
            adata.obs[col] = s.values
    coords = None
    try:
        from pycisTopic.utils import region_names_to_coordinates
        coords = region_names_to_coordinates(list(adata.var_names))
        adata.var["Chromosome"] = coords["Chromosome"].astype(str).values
        adata.var["Start"] = coords["Start"].values
        adata.var["End"] = coords["End"].values
    except Exception as exc:  # noqa: BLE001
        log.warning("could not annotate var with coordinates (%s); var_names "
                    "still carry chr:start-end", exc)
    adata.uns["cistopic"] = {
        "project": args.project,
        "n_topics": int(model.n_topic),
        "scale_factor": int(args.impute_scale_factor),
        "imputed_on_variable_regions": bool(args.impute_on_variable_regions),
        "source_matrix": str(f_mtx),
        "blacklist": str(f_bl),
    }
    adata.write_h5ad(h5ad_out, compression="gzip")
    log.info("wrote %s: %d cells x %d regions", h5ad_out, *adata.shape)
    del adata

    # =====================================================================
    # 8. per-cell-type DARs
    # =====================================================================
    if args.skip_dars:
        log.info("--skip-dars set: no DAR region sets")
    elif args.group_col is None:
        log.info("no --group-col given: skipping per-cell-type DARs (topic region "
                 "sets above are unaffected)")
    else:
        log.info("normalize_scores() then find_highly_variable_features()")
        # VERIFIED normalize_scores(imputed_acc, scale_factor=10**4)
        #   -> same class as input. NOTE it densifies (mtx.toarray()) if sparse.
        norm = normalize_scores(imputed, scale_factor=10 ** 4)
        # VERIFIED find_highly_variable_features(input_mat, min_disp=0.05,
        #   min_mean=0.0125, max_disp=np.inf, max_mean=3, n_bins=20,
        #   n_top_features=None, plot=True, save=None) -> list of feature names.
        #   Argument ORDER is (min_disp, min_mean, max_disp, max_mean) -- the
        #   tutorial passes them as keywords, and so must we: positional use
        #   would swap max_disp and max_mean.
        var_features = find_highly_variable_features(
            norm,
            min_disp=args.hvf_min_disp,
            min_mean=args.hvf_min_mean,
            max_mean=args.hvf_max_mean,
            max_disp=np.inf,
            n_bins=20,
            n_top_features=args.hvf_n_top,
            plot=False,
            save=str(out_dir / "hvf_dispersion.pdf"),
        )
        log.info("  %d variable regions", len(var_features))
        if len(var_features) == 0:
            die("find_highly_variable_features returned nothing; loosen "
                "--hvf-min-disp/--hvf-min-mean or set --hvf-n-top")
        (out_dir / "variable_regions.txt").write_text(
            "\n".join(map(str, var_features)) + "\n")

        log2fc = args.dar_log2fc if args.dar_log2fc is not None else float(np.log2(1.5))
        ray_tmp = Path(args.ray_temp_dir) if args.ray_temp_dir else tmp_path / "ray"
        # VERIFIED find_diff_features(cistopic_obj, imputed_features_obj, variable,
        #   var_features=None, contrasts=None, adjpval_thr=0.05,
        #   log2fc_thr=np.log2(1.5), split_pattern='___', n_cpu=1, **kwargs)
        #   -> dict {group_name: DataFrame[Log2FC, Adjusted_pval, Contrast]}
        # **kwargs goes to ray.init(); _temp_dir is a ray argument (TRAP 7).
        # With contrasts=None it builds one-vs-rest per level of `variable`.
        extra = {}
        if args.dar_n_cpu > 1:
            ray_tmp.mkdir(parents=True, exist_ok=True)
            extra["_temp_dir"] = str(ray_tmp)
            log.info("find_diff_features with n_cpu=%d -> ray.init(_temp_dir=%s)",
                     args.dar_n_cpu, ray_tmp)
        else:
            log.info("find_diff_features with n_cpu=1 -> ray is not initialised")

        log.info("find_diff_features(variable=%r, adjpval_thr=%g, log2fc_thr=%g)",
                 args.group_col, args.dar_adjpval, log2fc)
        try:
            markers_dict = find_diff_features(
                cistopic_obj,
                imputed,
                variable=args.group_col,
                var_features=var_features,
                contrasts=None,
                adjpval_thr=args.dar_adjpval,
                log2fc_thr=log2fc,
                n_cpu=args.dar_n_cpu,
                split_pattern=args.split_pattern,
                **extra,
            )
        except Exception as exc:  # noqa: BLE001
            die(f"find_diff_features failed: {type(exc).__name__}: {exc}\n"
                f"  If this mentions ray, retry with --dar-n-cpu 1 (serial, no "
                f"ray) or point --ray-temp-dir at node-local scratch.")

        dar_dir = rs_dir / f"DARs_{args.group_col}"
        dar_dir.mkdir(parents=True, exist_ok=True)
        summary = []
        for group, df in markers_dict.items():
            safe = "".join(c if (c.isalnum() or c in "._-") else "_" for c in str(group))
            n = write_region_set_bed(df.index, dar_dir / f"{safe}.bed",
                                     f"DARs/{group}")
            df.to_csv(dar_dir / f"{safe}.tsv", sep="\t")
            summary.append({"group": group, "file_stem": safe, "n_dars": n})
            log.info("  %-30s %6d DARs", str(group), n)
        pd.DataFrame(summary).to_csv(dar_dir / "dar_summary.tsv",
                                     sep="\t", index=False)
        log.info("wrote %d DAR region sets to %s", len(summary), dar_dir)
        del norm

    # =====================================================================
    # 9. manifest
    # =====================================================================
    lines = [
        "key\tvalue",
        f"script\t01_cistopic/run_cistopic.py",
        f"pycisTopic_source\tsrc/pycisTopic @ 787ce42 (version 2.0a)",
        f"project\t{args.project}",
        f"n_cells\t{len(cistopic_obj.cell_names)}",
        f"n_regions\t{len(cistopic_obj.region_names)}",
        f"n_topics_scanned\t{','.join(map(str, args.n_topics))}",
        f"n_topics_selected\t{model.n_topic}",
        f"n_iter\t{args.n_iter}",
        f"mallet_binary\t{mallet_bin}",
        f"MALLET_MEMORY\t{os.environ.get('MALLET_MEMORY')}",
        f"mallet_tmp_path\t{tmp_path}",
        f"blacklist\t{f_bl}",
        f"group_col\t{args.group_col}",
        f"impute_scale_factor\t{args.impute_scale_factor}",
        f"h5ad\t{h5ad_out}",
        f"cistopic_obj\t{f_obj}",
        f"region_set_folder\t{rs_dir}",
        f"runtime_min\t{(time.time() - t0) / 60:.1f}",
    ]
    (out_dir / "run_manifest.tsv").write_text("\n".join(lines) + "\n")

    log.info("DONE in %.1f min", (time.time() - t0) / 60)
    log.info("cisTopic object : %s", f_obj)
    log.info("h5ad for GLUE   : %s", h5ad_out)
    log.info("region sets     : %s", rs_dir)


if __name__ == "__main__":
    main()
