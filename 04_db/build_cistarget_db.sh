#!/usr/bin/env bash
#
# build_cistarget_db.sh -- build a custom cisTarget motif database from an ATAC
#                          consensus peak set (hg38, ArchR fixed-width peaks).
#
# ============================================================================
# WHAT THIS DOES
# ============================================================================
#   1. Validates every input and every required binary, up front, loudly.
#   2. Extracts a region FASTA with `bedtools getfasta`, naming each sequence
#      "chr:start-end" so region IDs match your ATAC var_names exactly.
#   3. Builds the motif ID list from the Cluster-Buster motif directory.
#   4. Runs create_cistarget_motif_databases.py, which emits BOTH
#        <prefix>.regions_vs_motifs.rankings.feather
#        <prefix>.regions_vs_motifs.scores.feather
#      (plus the motifs_vs_regions scores intermediate -- see OUTPUTS below).
#
# ============================================================================
# UPSTREAM FACTS THIS SCRIPT IS BUILT ON
# (verified by reading github.com/aertslab/create_cisTarget_databases, not
#  from memory -- if you update the clone, re-verify these)
# ============================================================================
#
# * THE ENTRY POINT IS `create_cistarget_motif_databases.py`.
#   There is NO `create_cistarget_databases.py` in the repo. The repo ships
#   `create_cistarget_motif_databases.py` (motifs, what you want) and
#   `create_cistarget_track_databases.py` (bigWig ChIP tracks, not this).
#
# * THE REPO IS SCRIPTS, NOT A PIP PACKAGE.
#   `pyproject.toml` exists but only configures tooling; there is no install
#   step and no console_scripts entry point. You clone it and call the .py
#   files by path. Hence CTDB_DIR below.
#
# * FLAGS (from create_cistarget_motif_databases.py argparse):
#     -f/--fasta FASTA              (required) region FASTA
#     -M/--motifs_dir DIR           (required) dir of Cluster-Buster .cb motifs
#     -m/--motifs FILE              (required) list of motif IDs, one per line
#     -o/--output PREFIX            (required) feather output prefix
#     -c/--cbust PATH               Cluster-Buster binary   [default: "cbust"]
#     -t/--threads N                threads                 [default: 1]
#     -p/--partial CUR TOTAL        split motif list into TOTAL parts
#     -b/--bgpadding N              bp of bg padding in FASTA [default: 0]
#     -l/--mask                     treat lowercase (repeat-masked) as N
#     -s/--seed N                   seed for breaking ranking ties
#     -F/--fasta-original-species   cross-species only
#     -5/--md5 FILE                 motif-MD5 -> motif-ID map
#     -g/--genes REGEX              collapse regions to genes (NOT for ATAC)
#     --min N / --max N             filter motif files by motif count
#
# * CBUST IS A HARD REQUIREMENT, RESOLVED VIA `shutil.which`.
#   create_cistarget_motif_databases.py does:
#       cluster_buster_path = shutil.which(args.cluster_buster_path)
#       if not cluster_buster_path: <print error>; sys.exit(1)
#   So the value of -c must be findable AND executable: either a bare name on
#   $PATH or a path that `which` resolves. A path that exists but lacks the
#   executable bit fails, as does a relative path not on $PATH. This script
#   checks all of that before submitting hours of compute.
#
#   Every region is scored by invoking cbust once per motif file as:
#       cbust -f 4 -c 0.0 -r 10000 -b <bgpadding> -t 1 [-l] <motif.cb> <fasta>
#   (clusterbuster.py). Note -t 1: cbust itself is single-threaded here; the
#   -t you pass to the Python script sets how many cbust processes run in
#   parallel. Get the binary from
#       https://resources.aertslab.org/cistarget/programs/cbust   (~3.1 MB)
#   then `chmod a+x cbust`.
#
# * MOTIF LIST SEMANTICS (clusterbuster.py::
#   get_motif_id_to_filename_and_nbr_motifs_dict):
#   Each line is a motif ID; a trailing ".cb" is stripped if present; the file
#   read is "${motifs_dir}/${motif_id}.cb" and a missing file is a hard error.
#   So the list must contain motif IDs (or basenames), NOT paths.
#
# * OUTPUT FILENAMES are assembled by DatabaseTypes.create_db_filename() as
#       {prefix}.{column_kind}_vs_{row_kind}.{scores|rankings}.feather
#   and write_db() adds the row kind as a named column holding the row labels.
#   In the regions_vs_motifs files the region IDs are the COLUMN NAMES and
#   there is one extra column literally named "motifs". That is the layout
#   pycistarget reads.
#
#   CAUTION: if you pass --min/--max (or -p), the script silently inserts
#   ".part_0001_of_0001" and/or ".min_N_to_max_..." into the prefix, so the
#   filenames change. This script therefore does not pass --min/--max.
#
# ============================================================================
# OUTPUTS (for a 300k-peak x ~5.9k-motif build)
# ============================================================================
#   <prefix>.motifs_vs_regions.scores.feather    ~3.4 GB   intermediate
#   <prefix>.regions_vs_motifs.scores.feather    ~3.4 GB   keep (DEM/scores)
#   <prefix>.regions_vs_motifs.rankings.feather  ~2-4 GB   keep (SCENIC+ uses)
#
#   The motifs_vs_regions RANKINGS file is deliberately NOT written upstream
#   (a commented-out write_db call notes it can take ~1.5 h for 1M regions
#   because the array is C-ordered and Feather writes column-wise). Do not
#   wait for it.
#
#   DISK: budget ~4x the final rankings file, i.e. 25-30 GB free for a 300k
#   peak build, plus ~1 GB for the region FASTA (300k x 501 bp + headers).
#   Feather is zstd level 6 by default; scores compress better than rankings.
#
# ============================================================================
# WALL-TIME AND CPU SCALE
# ============================================================================
#   Cost model: one cbust invocation per motif file over ALL regions, so total
#   work ~ (n_regions x total_bp) x n_motif_files, parallel over -t.
#
#   For 300k x 501 bp regions (~150 Mbp) and the v10nr_clust collection
#   (10,249 singleton .cb files; the precomputed hg38 SCREEN db contains 5,876
#   scored motif entries after clustering):
#
#     threads |  approximate wall time
#     --------+-----------------------
#        8    |  40-70 h      (not recommended)
#       16    |  20-35 h
#       32    |  10-18 h      <- reasonable target
#       64    |   6-10 h
#
#   Highly machine-dependent (cbust throughput varies ~2x with the math
#   library -- the AMD LibM build is roughly twice as fast as glibc). Treat
#   these as order-of-magnitude and set your SLURM --time with headroom.
#
#   MEMORY is the real constraint. The README gives:
#       scores db   = 4 bytes x n_regions x n_motifs
#       rankings db = 4 bytes x n_regions x n_motifs   (>32768 regions)
#   and warns you need ~3x the database size in RAM while running.
#       300,000 regions x 5,876 motifs x 4 B = 7.1 GB per matrix
#       -> ~21 GB working set; request 64 GB to be safe.
#       400,000 x 10,249 x 4 B = 16.4 GB -> ~49 GB; request 128 GB.
#   If you cannot get that much RAM, use -p/--partial to build the scores db
#   in parts, then combine with
#     combine_partial_regions_or_genes_vs_motifs_scores_cistarget_dbs.py and
#     convert_motifs_or_tracks_vs_regions_or_genes_scores_to_rankings_cistarget_dbs.py
#   (this script does not orchestrate the partial path; see repo README).
#
# ============================================================================
# USAGE
# ============================================================================
#   ./build_cistarget_db.sh \
#       --bed        consensus_peaks.bed \
#       --genome-fa  hg38.fa \
#       --motifs-dir /path/to/v10nr_clust_public/singletons \
#       --out-prefix /path/to/db/hg38_myAML_v10nr_clust \
#       --ctdb-dir   src/ctdb \
#       --threads    32
#
#   Add --dry-run to print the full plan (including the exact command line)
#   without executing anything.
#
set -o errexit
set -o nounset
set -o pipefail

SCRIPT_NAME="$(basename "${0}")"

# --------------------------------------------------------------------------- #
# Defaults
# --------------------------------------------------------------------------- #
BED=""
GENOME_FA=""
MOTIFS_DIR=""
OUT_PREFIX=""
CTDB_DIR="${CTDB_DIR:-}"
MOTIFS_LIST=""
CBUST="${CBUST:-cbust}"
THREADS=1
BG_PADDING=0
MASK=0
SEED=42
FASTA_OUT=""
KEEP_FASTA=1
DRY_RUN=0

# --------------------------------------------------------------------------- #
# Logging / failure helpers
# --------------------------------------------------------------------------- #
log()  { printf '[%s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "${*}" >&2; }
warn() { printf '[%s] WARNING: %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "${*}" >&2; }

die() {
    printf '\n' >&2
    printf '========================================================================\n' >&2
    printf 'ERROR (%s): %s\n' "${SCRIPT_NAME}" "${1}" >&2
    shift
    while [ "${#}" -gt 0 ]; do
        printf '  -> %s\n' "${1}" >&2
        shift
    done
    printf '========================================================================\n' >&2
    exit 1
}

usage() {
    sed -n '2,/^set -o errexit/p' "${0}" | sed 's/^# \{0,1\}//; $d'
    exit "${1:-0}"
}

# --------------------------------------------------------------------------- #
# Argument parsing
# --------------------------------------------------------------------------- #
while [ "${#}" -gt 0 ]; do
    case "${1}" in
        --bed)          BED="${2:?--bed needs a value}"; shift 2 ;;
        --genome-fa)    GENOME_FA="${2:?--genome-fa needs a value}"; shift 2 ;;
        --motifs-dir)   MOTIFS_DIR="${2:?--motifs-dir needs a value}"; shift 2 ;;
        --out-prefix)   OUT_PREFIX="${2:?--out-prefix needs a value}"; shift 2 ;;
        --ctdb-dir)     CTDB_DIR="${2:?--ctdb-dir needs a value}"; shift 2 ;;
        --motifs-list)  MOTIFS_LIST="${2:?--motifs-list needs a value}"; shift 2 ;;
        --cbust)        CBUST="${2:?--cbust needs a value}"; shift 2 ;;
        --threads)      THREADS="${2:?--threads needs a value}"; shift 2 ;;
        --bgpadding)    BG_PADDING="${2:?--bgpadding needs a value}"; shift 2 ;;
        --mask)         MASK=1; shift ;;
        --seed)         SEED="${2:?--seed needs a value}"; shift 2 ;;
        --fasta)        FASTA_OUT="${2:?--fasta needs a value}"; shift 2 ;;
        --no-keep-fasta) KEEP_FASTA=0; shift ;;
        --dry-run)      DRY_RUN=1; shift ;;
        -h|--help)      usage 0 ;;
        *)              die "unknown option: ${1}" \
                            "run '${SCRIPT_NAME} --help' for usage" ;;
    esac
done

# --------------------------------------------------------------------------- #
# Validate required arguments
# --------------------------------------------------------------------------- #
[ -n "${BED}" ]        || die "--bed is required" \
    "path to your consensus peak BED (ArchR fixed-width reproducible peaks)"
[ -n "${GENOME_FA}" ]  || die "--genome-fa is required" \
    "path to the hg38 genome FASTA (must match the build your peaks were called on)"
[ -n "${MOTIFS_DIR}" ] || die "--motifs-dir is required" \
    "directory of Cluster-Buster .cb motif files, e.g. v10nr_clust_public/singletons" \
    "get it from https://resources.aertslab.org/cistarget/motif_collections/v10nr_clust_public/v10nr_clust_public.zip"
[ -n "${OUT_PREFIX}" ] || die "--out-prefix is required" \
    "e.g. /path/to/db/hg38_myAML_v10nr_clust (no .feather suffix)"
[ -n "${CTDB_DIR}" ]   || die "--ctdb-dir is required (or set \$CTDB_DIR)" \
    "this repo is SCRIPTS, not a pip package -- you must point at a clone:" \
    "git clone --depth 1 https://github.com/aertslab/create_cisTarget_databases.git src/ctdb"

# --------------------------------------------------------------------------- #
# Validate inputs exist and are sane
# --------------------------------------------------------------------------- #
[ -f "${BED}" ] || die "peak BED not found: ${BED}"
[ -s "${BED}" ] || die "peak BED is empty: ${BED}"

case "${BED}" in
    *.gz|*.bgz) die "peak BED is compressed: ${BED}" \
        "bedtools getfasta needs a plain BED here; decompress it first:" \
        "  zcat '${BED}' > peaks.bed" ;;
esac

[ -f "${GENOME_FA}" ] || die "genome FASTA not found: ${GENOME_FA}"
case "${GENOME_FA}" in
    *.gz|*.bgz) die "genome FASTA is bgzipped/gzipped: ${GENOME_FA}" \
        "bedtools getfasta requires an uncompressed, .fai-indexed FASTA" ;;
esac

if [ ! -f "${GENOME_FA}.fai" ]; then
    warn "no FASTA index at ${GENOME_FA}.fai -- bedtools will try to create one"
    warn "(this needs write permission in $(dirname "${GENOME_FA}"); if that dir"
    warn " is read-only, run 'samtools faidx ${GENOME_FA}' somewhere writable)"
fi

[ -d "${MOTIFS_DIR}" ] || die "motifs dir not found: ${MOTIFS_DIR}"

CTDB_SCRIPT="${CTDB_DIR}/create_cistarget_motif_databases.py"
[ -d "${CTDB_DIR}" ] || die "create_cisTarget_databases clone not found: ${CTDB_DIR}" \
    "git clone --depth 1 https://github.com/aertslab/create_cisTarget_databases.git '${CTDB_DIR}'"
[ -f "${CTDB_SCRIPT}" ] || die "entry point not found: ${CTDB_SCRIPT}" \
    "expected 'create_cistarget_motif_databases.py' in the clone." \
    "NOTE: there is no 'create_cistarget_databases.py' in this repo -- the" \
    "motif entry point is create_cistarget_motif_databases.py."

# --------------------------------------------------------------------------- #
# Validate BED shape: fixed-width, no track lines, plausible hg38 contigs
# --------------------------------------------------------------------------- #
if head -n 1 "${BED}" | grep -qE '^(track|browser|#)'; then
    die "peak BED starts with a track/browser/# header line: ${BED}" \
        "bedtools getfasta will choke; strip it:" \
        "  grep -vE '^(track|browser|#)' in.bed > peaks.bed"
fi

N_PEAKS="$(grep -cvE '^([[:space:]]*$|track|browser|#)' "${BED}" || true)"
[ "${N_PEAKS}" -gt 0 ] || die "no usable records in ${BED}"

N_COLS="$(awk 'NR==1{print NF}' "${BED}")"
[ "${N_COLS}" -ge 3 ] || die "peak BED has only ${N_COLS} column(s): ${BED}" \
    "need at least chrom/start/end"

# Peak widths. ArchR reproducible peaks are fixed-width (501 bp by default).
WIDTH_INFO="$(awk -F'\t' '
    !/^(track|browser|#)/ && NF>=3 {
        w = $3 - $2
        if (w <= 0) { bad++; next }
        if (n++ == 0) { min = max = w }
        if (w < min) min = w
        if (w > max) max = w
        sum += w
        seen[w] = 1
    }
    END {
        k = 0; for (x in seen) k++
        printf "%d %d %d %.1f %d %d", n, min, max, (n ? sum/n : 0), k, bad+0
    }' "${BED}")"
read -r W_N W_MIN W_MAX W_MEAN W_UNIQ W_BAD <<< "${WIDTH_INFO}"

[ "${W_BAD}" -eq 0 ] || die "${W_BAD} record(s) in ${BED} have end <= start" \
    "fix or drop these before building"

# Contig naming: the DB region IDs inherit whatever is in column 1, and
# pycistarget matches them against your ATAC var_names by string. A 'chr'
# mismatch here produces a database that silently matches nothing.
if ! awk -F'\t' '!/^(track|browser|#)/ && NF>=3 {print $1; exit}' "${BED}" \
        | grep -q '^chr'; then
    warn "first contig in ${BED} does not start with 'chr'."
    warn "hg38 region IDs are conventionally 'chr1:...'. If your ATAC"
    warn "var_names are 'chr1:...' but the BED says '1', the resulting"
    warn "database will match NOTHING. Verify before spending compute."
fi

N_NONSTD="$(awk -F'\t' '!/^(track|browser|#)/ && NF>=3 && $1 ~ /_|^chrEBV|^chrUn/ {n++} END{print n+0}' "${BED}")"
if [ "${N_NONSTD}" -gt 0 ]; then
    warn "${N_NONSTD} peak(s) are on scaffolds/alts/chrEBV."
    warn "These are fine to score but are absent from most precomputed DBs."
fi

# --------------------------------------------------------------------------- #
# Validate binaries
# --------------------------------------------------------------------------- #
command -v bedtools >/dev/null 2>&1 || die "bedtools not found on \$PATH" \
    "conda install -c bioconda bedtools"

PYTHON_BIN="${PYTHON_BIN:-python3}"
command -v "${PYTHON_BIN}" >/dev/null 2>&1 \
    || die "${PYTHON_BIN} not found on \$PATH" \
        "activate the create_cistarget_databases conda env first:" \
        "  conda create -n create_cistarget_databases 'python=3.10' 'numpy=1.21' \\" \
        "      'pandas>=1.4.1' 'pyarrow>=7.0.0' 'numba>=0.55.1' python-flatbuffers" \
        "  conda activate create_cistarget_databases"

MISSING_PY=""
for mod in numpy pandas pyarrow numba flatbuffers; do
    "${PYTHON_BIN}" -c "import ${mod}" >/dev/null 2>&1 || MISSING_PY="${MISSING_PY} ${mod}"
done
[ -z "${MISSING_PY}" ] || die "missing Python module(s):${MISSING_PY}" \
    "create_cistarget_motif_databases.py imports these; install into the env:" \
    "  conda install -c conda-forge 'numpy=1.21' 'pandas>=1.4.1' 'pyarrow>=7.0.0' \\" \
    "      'numba>=0.55.1' python-flatbuffers"

# -- PANDAS VERSION GUARD (this one HANGS, it does not error) -------------- #
# cistarget_db.py::update_scores_for_motif_or_track writes results in place via
#     self.df.to_numpy()[rows, col] = scores
# Under pandas >= 2 copy-on-write semantics, DataFrame.to_numpy() returns a
# READ-ONLY array, so that assignment raises
#     ValueError: assignment destination is read-only
# The assignment happens inside a multiprocessing.Pool *callback*, so the
# exception is raised on the pool's result-handling thread, is never propagated
# to the main thread, and the job HANGS FOREVER at
#     "Initialize dataframe (N regions x M motifs) ..."
# with no error and no output. Verified: pandas 3.0.5 hangs; pandas 1.5.3
# completes normally. On a cluster this burns your entire wall-time allocation
# and produces nothing.
#
# The upstream README pins 'pandas>=1.4.1' with 'python=3.10', which resolved
# to a 1.x pandas when it was written. Pin the upper bound explicitly.
PANDAS_VER="$("${PYTHON_BIN}" -c 'import pandas; print(pandas.__version__)' 2>/dev/null || echo "0")"
PANDAS_MAJOR="${PANDAS_VER%%.*}"
if [ "${PANDAS_MAJOR}" -ge 2 ] 2>/dev/null; then
    # Confirm the actual failure mode rather than trusting the version number:
    # some builds/flags change copy-on-write behaviour.
    WRITEABLE="$("${PYTHON_BIN}" -c '
import pandas as pd, numpy as np
print(pd.DataFrame(np.zeros((2,2), dtype=np.float32)).to_numpy().flags.writeable)
' 2>/dev/null || echo "False")"
    if [ "${WRITEABLE}" != "True" ]; then
        die "pandas ${PANDAS_VER} is INCOMPATIBLE with create_cisTarget_databases" \
            "DataFrame.to_numpy() is read-only under pandas >= 2 copy-on-write." \
            "cistarget_db.py assigns through it inside a multiprocessing callback," \
            "so the failure is swallowed by the pool's result thread and the job" \
            "HANGS INDEFINITELY at 'Initialize dataframe ...' -- no error, no output," \
            "no database, and your whole SLURM wall-time allocation is consumed." \
            "Install a pandas 1.x environment:" \
            "  conda create -n create_cistarget_databases 'python=3.10' 'numpy=1.24' \\" \
            "      'pandas=1.5' 'pyarrow=11' 'numba=0.57' python-flatbuffers" \
            "  conda activate create_cistarget_databases" \
            "(verified working: pandas 1.5.3. Verified hanging: pandas 3.0.5.)"
    fi
    warn "pandas ${PANDAS_VER} is >= 2 but to_numpy() is still writeable here;"
    warn "proceeding, but pandas 1.5.x is the tested configuration."
fi

# Preflight the repo's own sibling-module imports (cistarget_db, clusterbuster,
# orderstatistics, feather_v1_or_v2). These are flat imports of files in the
# clone, which resolve via the script-directory sys.path entry -- an entry that
# PYTHONSAFEPATH=1 / `python -P` removes. Catch it in 1 second instead of after
# hours of scoring.
if ! PYTHONPATH="$(cd "${CTDB_DIR}" && pwd)${PYTHONPATH:+:${PYTHONPATH}}" \
        "${PYTHON_BIN}" -c 'import cistarget_db, clusterbuster' >/dev/null 2>&1; then
    die "cannot import the create_cisTarget_databases sibling modules from ${CTDB_DIR}" \
        "'import cistarget_db, clusterbuster' failed even with PYTHONPATH set." \
        "Check that the clone is complete (it must contain cistarget_db.py," \
        "clusterbuster.py, orderstatistics.py, feather_v1_or_v2.py):" \
        "  ls '${CTDB_DIR}'" \
        "and that numpy/pandas/pyarrow/numba import cleanly in ${PYTHON_BIN}." \
        "Diagnose with:" \
        "  PYTHONPATH='${CTDB_DIR}' ${PYTHON_BIN} -c 'import cistarget_db'"
fi

# -- THE CBUST CHECK ------------------------------------------------------- #
# Mirror exactly what create_cistarget_motif_databases.py does:
#     cluster_buster_path = shutil.which(args.cluster_buster_path)
#     if not cluster_buster_path: error; sys.exit(1)
# shutil.which() requires the target to be resolvable AND executable. Doing
# this check here turns a 6-hour-job failure into a 1-second failure.
CBUST_RESOLVED="$("${PYTHON_BIN}" - "${CBUST}" <<'PYEOF' || true
import shutil, sys
print(shutil.which(sys.argv[1]) or "")
PYEOF
)"
if [ -z "${CBUST_RESOLVED}" ]; then
    die "Cluster-Buster ('${CBUST}') could not be resolved by shutil.which()" \
        "create_cistarget_motif_databases.py resolves -c/--cbust with" \
        "shutil.which(), so the value must be a bare name on \$PATH or a path" \
        "that resolves AND has the executable bit set." \
        "Install the precompiled binary:" \
        "  wget https://resources.aertslab.org/cistarget/programs/cbust  # ~3.1 MB" \
        "  chmod a+x cbust" \
        "  mv cbust \"\${CONDA_PREFIX}/bin/cbust\"   # or pass --cbust /abs/path/cbust" \
        "Note: a path like './cbust' that exists but is NOT executable, and a" \
        "bare relative name not on \$PATH, both fail this check."
fi
[ -x "${CBUST_RESOLVED}" ] || die "cbust resolved to '${CBUST_RESOLVED}' but is not executable" \
    "chmod a+x '${CBUST_RESOLVED}'"

# Confirm it actually runs. cbust has no --version; -h returns non-zero on some
# builds, so treat "produced output" as the success signal.
if ! "${CBUST_RESOLVED}" -h >/dev/null 2>&1; then
    if [ -z "$("${CBUST_RESOLVED}" -h 2>&1 | head -c 1)" ]; then
        die "cbust at '${CBUST_RESOLVED}' did not execute" \
            "wrong architecture, missing libs, or a corrupt download." \
            "Check with: file '${CBUST_RESOLVED}' && ldd '${CBUST_RESOLVED}'"
    fi
fi

# --------------------------------------------------------------------------- #
# Motif list
# --------------------------------------------------------------------------- #
# Each line must be a motif ID whose file is "${MOTIFS_DIR}/${id}.cb"
# (clusterbuster.py strips a trailing ".cb" and hard-errors on a missing file).
N_CB="$(find "${MOTIFS_DIR}" -maxdepth 1 -name '*.cb' | wc -l | tr -d ' ')"
[ "${N_CB}" -gt 0 ] || die "no .cb motif files in ${MOTIFS_DIR}" \
    "expected Cluster-Buster format files, e.g. metacluster_1.1.cb" \
    "download + unzip:" \
    "  wget https://resources.aertslab.org/cistarget/motif_collections/v10nr_clust_public/v10nr_clust_public.zip" \
    "  unzip v10nr_clust_public.zip   # -> v10nr_clust_public/singletons/*.cb (10,249 files)"

GENERATED_LIST=0
if [ -z "${MOTIFS_LIST}" ]; then
    MOTIFS_LIST="${OUT_PREFIX}.motifs.lst"
    GENERATED_LIST=1
else
    [ -f "${MOTIFS_LIST}" ] || die "--motifs-list not found: ${MOTIFS_LIST}"
    [ -s "${MOTIFS_LIST}" ] || die "--motifs-list is empty: ${MOTIFS_LIST}"
fi

# --------------------------------------------------------------------------- #
# Derived paths, sizes, resource estimates
# --------------------------------------------------------------------------- #
OUT_DIR="$(dirname "${OUT_PREFIX}")"
[ -n "${FASTA_OUT}" ] || FASTA_OUT="${OUT_PREFIX}.regions.fa"

N_MOTIFS_EST="${N_CB}"
[ "${GENERATED_LIST}" -eq 1 ] || N_MOTIFS_EST="$(grep -cvE '^([[:space:]]*$|#)' "${MOTIFS_LIST}" || echo "${N_CB}")"

# scores/rankings matrix = 4 bytes x n_regions x n_motifs; ~3x that in RAM.
MATRIX_GB="$(awk -v r="${N_PEAKS}" -v m="${N_MOTIFS_EST}" \
    'BEGIN{printf "%.1f", 4*r*m/1024/1024/1024}')"
RAM_GB="$(awk -v g="${MATRIX_GB}" 'BEGIN{printf "%.0f", g*3}')"
FASTA_GB="$(awk -v r="${N_PEAKS}" -v w="${W_MEAN}" -v b="${BG_PADDING}" \
    'BEGIN{printf "%.2f", r*(w+2*b+30)/1024/1024/1024}')"
DISK_GB="$(awk -v g="${MATRIX_GB}" -v f="${FASTA_GB}" \
    'BEGIN{printf "%.0f", g*3+f+2}')"

OUT_RANKINGS="${OUT_PREFIX}.regions_vs_motifs.rankings.feather"
OUT_SCORES="${OUT_PREFIX}.regions_vs_motifs.scores.feather"
OUT_SCORES_T="${OUT_PREFIX}.motifs_vs_regions.scores.feather"

# --------------------------------------------------------------------------- #
# The command line we will run
# --------------------------------------------------------------------------- #
CTDB_ARGS=(
    -f "${FASTA_OUT}"
    -M "${MOTIFS_DIR}"
    -m "${MOTIFS_LIST}"
    -o "${OUT_PREFIX}"
    -c "${CBUST_RESOLVED}"
    -t "${THREADS}"
    -b "${BG_PADDING}"
    -s "${SEED}"
)
[ "${MASK}" -eq 1 ] && CTDB_ARGS+=(-l)

# --------------------------------------------------------------------------- #
# Plan
# --------------------------------------------------------------------------- #
WIDTH_DESC="variable (min ${W_MIN}, mean ${W_MEAN}, max ${W_MAX})"
[ "${W_UNIQ}" -eq 1 ] && WIDTH_DESC="fixed-width ${W_MIN} bp"

cat >&2 <<PLAN
========================================================================
build_cistarget_db.sh -- PLAN
========================================================================
INPUTS
  peak BED            : ${BED}
                        ${N_PEAKS} peaks, ${N_COLS} column(s), ${WIDTH_DESC}
  genome FASTA        : ${GENOME_FA}
  motifs dir          : ${MOTIFS_DIR}
                        ${N_CB} .cb files present
  motifs list         : ${MOTIFS_LIST}$([ "${GENERATED_LIST}" -eq 1 ] && echo "  (will be generated from motifs dir)")
                        ~${N_MOTIFS_EST} motifs to score
  create_cisTarget dir: ${CTDB_DIR}
  entry point         : ${CTDB_SCRIPT}
  cbust               : ${CBUST_RESOLVED}
  python              : $(command -v "${PYTHON_BIN}")

PARAMETERS
  threads (-t)        : ${THREADS}   (parallel cbust processes; cbust itself
                                      is invoked with -t 1 per process)
  bg padding (-b)     : ${BG_PADDING}
  mask lowercase (-l) : $([ "${MASK}" -eq 1 ] && echo yes || echo no)
  seed (-s)           : ${SEED}   (fixed -> reproducible rankings)

STEP 1  extract region FASTA
  bedtools getfasta -fi <genome> -bed <bed> -nameOnly -fo ${FASTA_OUT}
  Region names are forced to "chr:start-end" (BED col 1-3) so region IDs
  match your ATAC var_names exactly.

STEP 2  build motif list$([ "${GENERATED_LIST}" -eq 1 ] && echo "" || echo " (using supplied list)")

STEP 3  create database
  ${CTDB_SCRIPT} \\
$(printf '      %s\n' "${CTDB_ARGS[@]}" | paste -d' ' - - | sed 's/$/ \\/' | sed '$ s/ \\$//')

EXPECTED OUTPUTS
  ${OUT_RANKINGS}
  ${OUT_SCORES}
  ${OUT_SCORES_T}   (intermediate)
  NOTE: motifs_vs_regions.RANKINGS is intentionally never written upstream.

RESOURCE ESTIMATE  (${N_PEAKS} regions x ~${N_MOTIFS_EST} motifs)
  one matrix          : ~${MATRIX_GB} GB   (4 bytes x regions x motifs)
  peak RAM            : ~${RAM_GB} GB     (README: ~3x database size)
  region FASTA        : ~${FASTA_GB} GB
  free disk needed    : ~${DISK_GB} GB
  wall time           : ~10-18 h at 32 threads; ~20-35 h at 16.
                        Scales ~linearly in n_regions x n_motifs / threads.
========================================================================
PLAN

# --------------------------------------------------------------------------- #
# Disk space check
# --------------------------------------------------------------------------- #
if [ "${DRY_RUN}" -eq 0 ]; then
    mkdir -p "${OUT_DIR}" || die "cannot create output dir: ${OUT_DIR}"
    [ -w "${OUT_DIR}" ] || die "output dir is not writable: ${OUT_DIR}"

    AVAIL_GB="$(df -Pk "${OUT_DIR}" | awk 'NR==2{printf "%.0f", $4/1024/1024}')"
    log "free disk in ${OUT_DIR}: ${AVAIL_GB} GB (need ~${DISK_GB} GB)"
    if [ "${AVAIL_GB}" -lt "${DISK_GB}" ]; then
        die "not enough free disk in ${OUT_DIR}: ${AVAIL_GB} GB free, ~${DISK_GB} GB needed" \
            "the scores + rankings feathers plus the intermediate transpose are large" \
            "free space, or point --out-prefix at a bigger filesystem"
    fi
fi

if [ "${DRY_RUN}" -eq 1 ]; then
    log "--dry-run: all checks passed, nothing executed."
    log "Re-run without --dry-run to build."
    exit 0
fi

# --------------------------------------------------------------------------- #
# STEP 1: region FASTA
# --------------------------------------------------------------------------- #
# -nameOnly uses BED column 4 as the sequence name, so we synthesise column 4
# as "chr:start-end" from columns 1-3. We do NOT use bedtools' own `-name+`
# or default naming:
#   * default naming gives ">chr:start-end" but appends strand for BED6+ and
#     varies across bedtools versions;
#   * -name+ gives "name::chr:start-end".
# Forcing column 4 ourselves makes the region ID exactly reproducible and
# identical to the ArchR/AnnData var_names convention "chr:start-end".
if [ -s "${FASTA_OUT}" ] && [ "${KEEP_FASTA}" -eq 1 ]; then
    N_FA="$(grep -c '^>' "${FASTA_OUT}" || echo 0)"
    if [ "${N_FA}" -eq "${N_PEAKS}" ]; then
        log "STEP 1: reusing existing FASTA (${N_FA} sequences): ${FASTA_OUT}"
    else
        warn "existing ${FASTA_OUT} has ${N_FA} sequences but BED has ${N_PEAKS}; regenerating"
        rm -f "${FASTA_OUT}"
    fi
fi

if [ ! -s "${FASTA_OUT}" ]; then
    log "STEP 1: extracting region FASTA -> ${FASTA_OUT}"
    TMP_BED="$(mktemp "${OUT_DIR}/.ctdb_bed.XXXXXX")"
    TMP_FA="$(mktemp "${OUT_DIR}/.ctdb_fa.XXXXXX")"
    trap 'rm -f "${TMP_BED}" "${TMP_FA}"' EXIT

    awk -F'\t' -v OFS='\t' '
        !/^(track|browser|#)/ && NF >= 3 && ($3 - $2) > 0 {
            print $1, $2, $3, $1 ":" $2 "-" $3
        }' "${BED}" > "${TMP_BED}"

    N_TMP="$(wc -l < "${TMP_BED}" | tr -d ' ')"
    [ "${N_TMP}" -gt 0 ] || die "produced an empty BED for getfasta -- check ${BED}"

    # Region IDs must be unique or scores collide silently.
    N_UNIQ="$(cut -f4 "${TMP_BED}" | sort -u | wc -l | tr -d ' ')"
    [ "${N_UNIQ}" -eq "${N_TMP}" ] || die \
        "peak BED contains duplicate intervals (${N_TMP} records, ${N_UNIQ} unique)" \
        "duplicate region IDs collapse in the database and corrupt scores" \
        "de-duplicate first: sort -k1,1 -k2,2n -u"

    if ! bedtools getfasta \
            -fi "${GENOME_FA}" \
            -bed "${TMP_BED}" \
            -nameOnly \
            -fo "${TMP_FA}" 2> "${TMP_FA}.err"; then
        sed 's/^/  bedtools: /' "${TMP_FA}.err" >&2 || true
        rm -f "${TMP_FA}.err"
        die "bedtools getfasta failed" \
            "common causes:" \
            "  * contig names in the BED are absent from the genome FASTA" \
            "    (chr1 vs 1, or hg19 peaks against an hg38 FASTA)" \
            "  * peak coordinates run past the end of a contig" \
            "  * ${GENOME_FA}.fai missing and its directory not writable"
    fi
    rm -f "${TMP_FA}.err"

    N_SEQ="$(grep -c '^>' "${TMP_FA}" || echo 0)"
    [ "${N_SEQ}" -eq "${N_TMP}" ] || die \
        "FASTA has ${N_SEQ} sequences but BED had ${N_TMP} regions" \
        "bedtools skipped regions -- almost always contig-name or bounds problems"

    # A region that is all-N scores nothing and just wastes cbust time.
    N_ALLN="$(awk '/^>/{next} {if ($0 ~ /^[Nn]+$/) n++} END{print n+0}' "${TMP_FA}")"
    [ "${N_ALLN}" -eq 0 ] && log "  no all-N sequences" \
        || warn "  ${N_ALLN} sequence(s) are entirely N (assembly gaps); they will score 0"

    mv "${TMP_FA}" "${FASTA_OUT}"
    rm -f "${TMP_BED}"
    trap - EXIT
    log "STEP 1: done -- ${N_SEQ} sequences, $(du -h "${FASTA_OUT}" | cut -f1)"
    log "  first region ID: $(head -n 1 "${FASTA_OUT}")"
fi

# --------------------------------------------------------------------------- #
# STEP 2: motif list
# --------------------------------------------------------------------------- #
if [ "${GENERATED_LIST}" -eq 1 ]; then
    log "STEP 2: generating motif list -> ${MOTIFS_LIST}"
    # Motif IDs = .cb basenames with the extension stripped.
    find "${MOTIFS_DIR}" -maxdepth 1 -name '*.cb' -exec basename {} .cb \; \
        | LC_ALL=C sort > "${MOTIFS_LIST}"
    N_LIST="$(wc -l < "${MOTIFS_LIST}" | tr -d ' ')"
    [ "${N_LIST}" -gt 0 ] || die "generated an empty motif list from ${MOTIFS_DIR}"
    log "STEP 2: done -- ${N_LIST} motif IDs"
else
    log "STEP 2: using supplied motif list ${MOTIFS_LIST}"
    # Fail now, not 3 hours in, if a listed motif file is missing.
    MISSING="$(awk -v d="${MOTIFS_DIR}" '
        /^[[:space:]]*$/ || /^#/ { next }
        { id = $0; sub(/\.cb$/, "", id)
          f = d "/" id ".cb"
          if ((getline line < f) < 0) { print id; n++ }
          close(f)
          if (n >= 5) exit }' "${MOTIFS_LIST}")"
    [ -z "${MISSING}" ] || die \
        "motif file(s) listed in ${MOTIFS_LIST} are missing from ${MOTIFS_DIR}" \
        "first missing: $(echo "${MISSING}" | tr '\n' ' ')" \
        "clusterbuster.py raises OSError on the first missing '<id>.cb'"
fi

# --------------------------------------------------------------------------- #
# STEP 3: build the database
# --------------------------------------------------------------------------- #
log "STEP 3: running create_cistarget_motif_databases.py"
log "  ${CTDB_SCRIPT} $(printf '%s ' "${CTDB_ARGS[@]}")"
log "  this is the long step -- see the wall-time table in this script's header"

START_TS="$(date +%s)"
# PYTHONPATH is load-bearing. create_cistarget_motif_databases.py does
#     from cistarget_db import ...
#     from clusterbuster import ...
# i.e. flat imports of sibling modules in the clone (the repo is scripts, not
# an installed package). That normally works because CPython prepends the
# script's own directory to sys.path -- but NOT when PYTHONSAFEPATH=1 or
# `python -P` is in effect, which some cluster module files and hardened
# images set globally. The failure is an opaque
#     ModuleNotFoundError: No module named 'cistarget_db'
# hours into a job. Setting PYTHONPATH explicitly makes it work either way.
CTDB_DIR_ABS="$(cd "${CTDB_DIR}" && pwd)"
if [ -x "${CTDB_SCRIPT}" ] && [ -z "${PYTHONSAFEPATH:-}" ]; then
    PYTHONPATH="${CTDB_DIR_ABS}${PYTHONPATH:+:${PYTHONPATH}}" \
        "${CTDB_SCRIPT}" "${CTDB_ARGS[@]}"
else
    # Clone may have lost the executable bit (e.g. copied via a zip), or
    # PYTHONSAFEPATH is set and we must not rely on the shebang's plain python3.
    PYTHONPATH="${CTDB_DIR_ABS}${PYTHONPATH:+:${PYTHONPATH}}" \
        "${PYTHON_BIN}" "${CTDB_SCRIPT}" "${CTDB_ARGS[@]}"
fi
END_TS="$(date +%s)"
ELAPSED=$(( END_TS - START_TS ))

# --------------------------------------------------------------------------- #
# Verify the outputs we care about actually exist
# --------------------------------------------------------------------------- #
FAILED=0
for f in "${OUT_RANKINGS}" "${OUT_SCORES}"; do
    if [ -s "${f}" ]; then
        log "OK  $(du -h "${f}" | cut -f1)  ${f}"
    else
        warn "MISSING or empty: ${f}"
        FAILED=1
    fi
done

[ "${FAILED}" -eq 0 ] || die \
    "the build finished but an expected feather is missing" \
    "if you passed --min/--max or -p/--partial, the prefix gains a" \
    "'.part_0001_of_0001' / '.min_N_to_max_motifs' infix and the filenames differ" \
    "check: ls -la $(dirname "${OUT_PREFIX}")/$(basename "${OUT_PREFIX}")*"

printf '\n' >&2
log "BUILD COMPLETE in $(( ELAPSED / 3600 ))h $(( (ELAPSED % 3600) / 60 ))m"
log "  rankings (SCENIC+ / pycistarget input): ${OUT_RANKINGS}"
log "  scores   (DEM / raw CRM scores)       : ${OUT_SCORES}"
log ""
log "NEXT: confirm your peaks are representable in what you just built:"
log "  python 04_db/peak_overlap_audit.py \\"
log "      --peaks '${BED}' --db '${OUT_RANKINGS}' \\"
log "      --out-prefix '${OUT_PREFIX}.audit'"
log "(A custom DB built from these same peaks should report ~0% dropped."
log " Anything else means the region IDs do not line up.)"
